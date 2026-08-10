#!/usr/bin/env python3
"""
OSDR Metadata Enrichment: Query factor values for all Arabidopsis spaceflight
studies and merge with the parsed sample metadata.

This supplements mine_osdr.py by fetching structured factor values from the
OSDR query API (study.factor value.spaceflight, organism part, light, etc.)
which are not captured in the initial broad query.
"""
import urllib.request, urllib.parse, csv, io, os, re, time
import pandas as pd

API_BASE = "https://visualization.osdr.nasa.gov/biodata/api"
OUTPUT_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"


def fetch_csv(url, timeout=60):
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "text/csv"})
            resp = urllib.request.urlopen(req, timeout=timeout)
            raw = resp.read().decode("utf-8")
            reader = csv.DictReader(io.StringIO(raw))
            return list(reader)
        except Exception as e:
            if attempt < 2:
                time.sleep(2)
            else:
                print(f"  ERROR: {e}")
                return []


def enrich_metadata():
    """Fetch factor values for all spaceflight Col-0/WS studies and merge."""
    df = pd.read_csv(os.path.join(OUTPUT_DIR, "sample_metadata.tsv"), sep="\t")
    print(f"Loaded {len(df)} samples from sample_metadata.tsv")

    # All factor value fields to request
    factor_fields = [
        "study.factor value.spaceflight",
        "study.factor value.organism part",
        "study.factor value.genotype",
        "study.factor value.light",
        "study.factor value.altered gravity",
        "study.factor value.ecotype",
        "study.factor value.growth time",
        "study.factor value.age",
        "study.factor value.cultivar",
        "study.factor value.treatment",
        "study.characteristics.ecotype",
        "study.characteristics.genotype",
        "study.characteristics.cultivar",
        "study.characteristics.condition",
        "study.characteristics.spaceflight",
    ]

    all_studies = sorted(df["id.accession"].unique())
    enriched_rows = []

    for osd in all_studies:
        params = [("id.accession", osd), ("id.sample name", "")] + [(f, "") for f in factor_fields]
        qs = urllib.parse.urlencode([(k, v) for k, v in params])
        url = f"{API_BASE}/v2/query/metadata/?{qs}"
        rows = fetch_csv(url, timeout=60)
        if rows:
            enriched_rows.extend(rows)
            print(f"  {osd}: {len(rows)} samples enriched")

    enriched_df = pd.DataFrame(enriched_rows)
    print(f"\nEnriched {len(enriched_df)} samples with factor values")

    # Clean up factor values
    def clean_factor(v):
        if not v or str(v).lower() == "nan":
            return None
        return str(v).strip()

    # Map spaceflight factor to FLT/GC
    def map_flight(v):
        if not v:
            return None
        vl = str(v).lower()
        if "space" in vl or "flight" in vl or "flt" in vl:
            return "FLT"
        elif "ground" in vl or "control" in vl or "gc" in vl:
            return "GC"
        return None

    # Map organism part to normalized organ
    def map_organ(v):
        if not v:
            return None
        vl = str(v).lower()
        organ_map = {
            "plant roots": "root", "plant root": "root", "root": "root",
            "plant shoots": "shoot", "plant shoot": "shoot", "shoot": "shoot",
            "plant leaves": "leaf", "leaf": "leaf", "leaves": "leaf",
            "hypocotyl": "hypocotyl", "hypocotyl cell culture": "hypocotyl",
            "seedlings": "whole_seedling", "seedling": "whole_seedling",
            "whole organism": "whole_seedling", "whole plant": "whole_seedling",
            "cotyledon": "cotyledon", "rosette": "rosette",
        }
        return organ_map.get(vl, vl)

    # Map light factor
    def map_light(v):
        if not v:
            return None
        vl = str(v).lower()
        if "red" in vl and "light" in vl:
            return "red_light"
        elif "dark" in vl:
            return "dark"
        elif "light" in vl:
            return "light"
        return None

    # Map ecotype from genotype/characteristics
    def map_ecotype(row):
        for col in ["study.factor value.ecotype", "study.characteristics.ecotype",
                     "study.characteristics.cultivar", "study.factor value.cultivar"]:
            v = clean_factor(row.get(col, ""))
            if v:
                vl = v.lower()
                if "col" in vl and "0" in vl:
                    return "Col-0"
                if "ws" in vl or "wassilewskij" in vl or "wasselewskij" in vl:
                    return "WS"
                if "col-0" in vl and "ws" in vl:
                    return "Col-0+WS"
        # Check genotype for wild type
        gt = clean_factor(row.get("study.factor value.genotype", "") or row.get("study.characteristics.genotype", ""))
        if gt and "wild type" in gt.lower():
            return "Col-0"  # Most OSDR WT is Col-0
        return None

    # Map altered gravity
    def map_gravity(v):
        if not v:
            return None
        vl = str(v).lower()
        if "micro" in vl or "ug" in vl or "0g" in vl or "0 g" in vl:
            return "microgravity"
        elif "0.34" in vl or "mars" in vl:
            return "Mars_g"
        elif "0.16" in vl or "moon" in vl or "lunar" in vl:
            return "Moon_g"
        elif "1g" in vl or "1 g" in vl:
            return "1g"
        return v

    enriched_df["flight_factor"] = enriched_df["study.factor value.spaceflight"].apply(map_flight)
    enriched_df["organ_factor"] = enriched_df["study.factor value.organism part"].apply(map_organ)
    enriched_df["light_factor"] = enriched_df["study.factor value.light"].apply(map_light)
    enriched_df["gravity_factor"] = enriched_df["study.factor value.altered gravity"].apply(map_gravity)
    enriched_df["ecotype_factor"] = enriched_df.apply(map_ecotype, axis=1)
    enriched_df["genotype_factor"] = enriched_df.apply(
        lambda r: clean_factor(r.get("study.factor value.genotype", "") or r.get("study.characteristics.genotype", "")), axis=1)
    enriched_df["growth_time_factor"] = enriched_df["study.factor value.growth time"].apply(clean_factor)

    # Merge enriched factor values back into the main metadata
    merge_cols = ["id.accession", "id.sample name", "flight_factor", "organ_factor",
                  "light_factor", "gravity_factor", "ecotype_factor", "genotype_factor", "growth_time_factor"]
    enriched_subset = enriched_df[merge_cols].drop_duplicates(subset=["id.accession", "id.sample name"])

    # Merge
    merged = df.merge(enriched_subset, on=["id.accession", "id.sample name"], how="left")

    # Use factor values as primary, parsed as fallback
    merged["flight_final"] = merged["flight_factor"].fillna(merged["flight_parsed"])
    merged["organ_final2"] = merged["organ_factor"].fillna(merged["organ_final"])
    merged["light_final"] = merged["light_factor"].fillna(merged["light_parsed"])
    merged["ecotype_final2"] = merged["ecotype_factor"].fillna(merged["ecotype_final"])

    # Normalize organ_final2
    merged["organ_final2"] = merged["organ_final2"].str.lower().map(
        lambda x: {"plant roots": "root", "plant root": "root", "root": "root",
                   "plant shoots": "shoot", "shoot": "shoot",
                   "plant leaves": "leaf", "leaf": "leaf", "leaves": "leaf",
                   "hypocotyl": "hypocotyl", "seedlings": "whole_seedling",
                   "seedling": "whole_seedling", "whole organism": "whole_seedling",
                   "whole_seedling": "whole_seedling", "cotyledon": "cotyledon",
                   "rosette": "rosette", "unknown": "unknown"}.get(str(x).lower(), x)
    )

    # Drop old columns and rename
    merged = merged.drop(columns=["organ_final", "ecotype_final", "flight_parsed", "light_parsed"])
    merged = merged.rename(columns={"organ_final2": "organ_final", "ecotype_final2": "ecotype_final",
                                     "flight_final": "flight_final", "light_final": "light_final"})

    # Save enriched metadata
    merged.to_csv(os.path.join(OUTPUT_DIR, "sample_metadata_enriched.tsv"), sep="\t", index=False)
    print(f"\nSaved sample_metadata_enriched.tsv ({len(merged)} rows)")

    # Print final factorial coverage
    print("\n=== Final factorial coverage ===")
    for factor in ["hardware", "flight_final", "organ_final", "ecotype_final", "omics_type", "light_final", "gravity_factor"]:
        vc = merged[factor].value_counts(dropna=False)
        print(f"  {factor}: {vc.to_dict()}")

    # Cross-tabulation: hardware × flight × organ
    print("\n=== Hardware × Flight × Organ coverage ===")
    ct = merged.groupby(["hardware", "flight_final", "organ_final", "omics_type"]).size().reset_index(name="n")
    ct = ct.sort_values("n", ascending=False)
    print(ct.to_string(index=False))
    ct.to_csv(os.path.join(OUTPUT_DIR, "factorial_coverage_enriched.tsv"), sep="\t", index=False)

    return merged


if __name__ == "__main__":
    enrich_metadata()
