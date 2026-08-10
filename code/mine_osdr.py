#!/usr/bin/env python3
"""
OSDR Miner: Query NASA OSDR Biological Data API for Arabidopsis spaceflight
multiomics data, parse sample names into factorial factors, and build a
harmonized metadata table.

API base: https://visualization.osdr.nasa.gov/biodata/api/v2/
"""
import urllib.request, urllib.parse, json, csv, io, os, re, time
import pandas as pd

API_BASE = "https://visualization.osdr.nasa.gov/biodata/api"
OUTPUT_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def fetch_json(url, timeout=30):
    """Fetch JSON from URL with retry."""
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            resp = urllib.request.urlopen(req, timeout=timeout)
            return json.loads(resp.read())
        except Exception as e:
            if attempt < 2:
                time.sleep(2)
            else:
                print(f"  ERROR fetching {url}: {e}")
                return None


def fetch_csv(url, timeout=90):
    """Fetch CSV from query endpoint."""
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
                print(f"  ERROR fetching {url}: {e}")
                return []


def get_all_datasets():
    """Get list of all OSDR dataset accessions."""
    data = fetch_json(f"{API_BASE}/v2/datasets/")
    if data:
        return sorted(data.keys(), key=lambda x: int(x.split("-")[1]) if x.split("-")[1].isdigit() else 9999)
    return []


def query_arabidopsis_samples():
    """Query all Arabidopsis thaliana samples with key metadata fields."""
    params = [
        ("id.accession", ""),
        ("id.assay name", ""),
        ("id.sample name", ""),
        ("study.characteristics.organism", "arabidopsis thaliana"),
        ("study.characteristics.ecotype", ""),
        ("study.characteristics.tissue", ""),
        ("study.characteristics.organism part", ""),
        ("study.factors", ""),
    ]
    qs = urllib.parse.urlencode(params)
    url = f"{API_BASE}/v2/query/metadata/?{qs}"
    rows = fetch_csv(url, timeout=120)
    print(f"  Retrieved {len(rows)} Arabidopsis samples")
    return rows


def get_study_metadata(osd_id):
    """Fetch study-level metadata (platform, title, factors, project type)."""
    data = fetch_json(f"{API_BASE}/v2/dataset/{osd_id}/")
    if not data:
        return {}
    meta = data.get(osd_id, {}).get("metadata", {})

    def clean(v):
        if isinstance(v, list):
            return " | ".join(str(x) for x in v)
        return str(v) if v else ""

    return {
        "osd_id": osd_id,
        "study_title": clean(meta.get("study title", "")),
        "project_title": clean(meta.get("project title", "")),
        "experiment_platform": clean(meta.get("experiment platform", "")),
        "study_factor_name": clean(meta.get("study factor name", "")),
        "study_factor_type": clean(meta.get("study factor type", "")),
        "project_type": clean(meta.get("project type", "")),
        "flight_program": clean(meta.get("flight program", "")),
        "mission": clean(meta.get("mission", "")),
        "study_assay_technology_type": clean(meta.get("study assay technology type", "")),
        "study_assay_measurement_type": clean(meta.get("study assay measurement type", "")),
        "study_description": clean(meta.get("study description", ""))[:500],
    }


def get_study_assays(osd_id):
    """Get assay names for a study."""
    data = fetch_json(f"{API_BASE}/v2/dataset/{osd_id}/assays/")
    if not data:
        return []
    assays = data.get(osd_id, {}).get("assays", {})
    return list(assays.keys())


def classify_assay(assay_name):
    """Classify assay name into omics type."""
    al = assay_name.lower()
    if "rna-seq" in al or "rna_seq" in al or "rna sequencing" in al:
        return "RNAseq"
    elif "microarray" in al or "dna-microarray" in al:
        return "microarray"
    elif "bisulfite" in al or "methylation" in al:
        return "methylation"
    elif "mirna" in al or "micro-rna" in al or "micro_rna" in al:
        return "microRNA"
    elif "mass-spectrometry" in al or "proteom" in al:
        return "proteomics"
    elif "image" in al or "photography" in al:
        return "imaging"
    elif "glycom" in al or "cell-wall" in al or "carbohydr" in al:
        return "glycomics"
    elif "genome" in al:
        return "genome"
    elif "protein-dna" in al or "chip" in al:
        return "protein-DNA"
    else:
        return assay_name.split("_")[1] if "_" in assay_name else "other"


def classify_hardware(platform_str):
    """Classify experiment platform string into hardware category."""
    pl = platform_str.lower()
    if "bric-led" in pl or "bric led" in pl:
        return "BRIC-LED"
    elif "bric" in pl:
        return "BRIC"
    elif "veggie" in pl or "vps" in pl or "vegetable production" in pl:
        return "Veggie"
    elif "apex" in pl or "tages" in pl:
        return "APEX-TAGES"
    elif "emcs" in pl or "european modular" in pl:
        return "EMCS"
    elif "cara" in pl:
        return "CARA"
    elif "abrs" in pl or "advanced biological research" in pl:
        return "ABRS"
    elif "kft" in pl or "ksc fixation" in pl:
        return "KFT"
    elif "lpgc" in pl or "low pressure" in pl:
        return "LPGC"
    elif "spaceshiptwo" in pl or "new shepard" in pl or "virgin galactic" in pl:
        return "Suborbital"
    elif "f-104" in pl or "c-9" in pl or "parabolic" in pl:
        return "Parabolic"
    elif "stratospher" in pl or "boron" in pl:
        return "Balloon"
    else:
        return "other"


def parse_sample_name(sample_name):
    """Parse sample name to extract factorial factors using regex.

    Common patterns:
    - Atha_Col-0-PhyD_root_FLT_Alight_Rep1_GSM2493783_Day03
    - Atha_WS-0_Col-0_Hypocotyl_FLT_Rep1
    - Atha_Col-0_sl-pool_FLT_Rep1_R1-FL-A1
    - Atha_WT-Col-0_sl_FLT_Rep1_G1S1_membrane
    """
    s = sample_name
    sl = s.lower()

    # Ecotype
    ecotype = "unknown"
    if "col-0" in sl or "col0" in sl or "wt-col" in sl or "col_0" in sl:
        ecotype = "Col-0"
    if "ws-0" in sl or "ws0" in sl or "ws_0" in sl or "wassilewskij" in sl or "wasselewskij" in sl:
        ecotype = "WS"
    # If both present, mark as mixed
    if "col-0" in sl and ("ws-0" in sl or "ws0" in sl):
        ecotype = "Col-0+WS"

    # Organ/tissue
    organ = "unknown"
    if "root" in sl:
        organ = "root"
    elif "hypocotyl" in sl:
        organ = "hypocotyl"
    elif "shoot" in sl:
        organ = "shoot"
    elif "leaf" in sl or "leaves" in sl:
        organ = "leaf"
    elif "cotyledon" in sl:
        organ = "cotyledon"
    elif "seedling" in sl or "sl-pool" in sl or "sl_pool" in sl or "sl_" in sl or "whole" in sl:
        organ = "whole_seedling"
    elif "rosette" in sl:
        organ = "rosette"

    # Flight condition
    flight = "unknown"
    if "_flt" in sl or "flt_" in sl or "flt-" in sl or "flight" in sl:
        flight = "FLT"
    elif "_gc" in sl or "gc_" in sl or "gc-" in sl or "ground" in sl:
        flight = "GC"

    # Light condition
    light = "unspecified"
    if "alight" in sl or "a_light" in sl or "_light" in sl or "light_" in sl:
        light = "light"
    if "blind" in sl or "dark" in sl or "blight" in sl:
        light = "dark"
    if "red" in sl and "light" in sl:
        light = "red_light"

    # Replicate
    rep_match = re.search(r"rep(\d+)", sl, re.IGNORECASE)
    replicate = int(rep_match.group(1)) if rep_match else None

    # Timepoint/age
    day_match = re.search(r"day(\d+)", sl, re.IGNORECASE)
    day = int(day_match.group(1)) if day_match else None

    return {
        "ecotype_parsed": ecotype,
        "organ_parsed": organ,
        "flight_parsed": flight,
        "light_parsed": light,
        "replicate_parsed": replicate,
        "day_parsed": day,
    }


def main():
    print("=" * 70)
    print("OSDR Miner: Arabidopsis spaceflight multiomics data retrieval")
    print("=" * 70)

    # Step 1: Query all Arabidopsis samples
    print("\n[1/5] Querying all Arabidopsis samples from OSDR...")
    samples = query_arabidopsis_samples()
    if not samples:
        print("ERROR: No samples retrieved. Exiting.")
        return

    samples_df = pd.DataFrame(samples)
    print(f"  Columns: {list(samples_df.columns)}")

    # Step 2: Get study-level metadata for all studies
    print("\n[2/5] Fetching study-level metadata...")
    all_studies = sorted(samples_df["id.accession"].unique(),
                         key=lambda x: int(x.split("-")[1]) if x.split("-")[1].isdigit() else 9999)
    print(f"  {len(all_studies)} unique studies")

    study_meta_list = []
    for i, osd in enumerate(all_studies):
        meta = get_study_metadata(osd)
        meta["assays"] = get_study_assays(osd)
        meta["omics_types"] = sorted(set(classify_assay(a) for a in meta["assays"]))
        meta["hardware"] = classify_hardware(meta["experiment_platform"])
        study_meta_list.append(meta)
        if (i + 1) % 10 == 0:
            print(f"  ... {i+1}/{len(all_studies)} studies processed")

    study_meta_df = pd.DataFrame(study_meta_list)
    print(f"  Study metadata retrieved for {len(study_meta_df)} studies")

    # Step 3: Parse sample names for factorial factors
    print("\n[3/5] Parsing sample names for factorial factors...")
    parsed = samples_df["id.sample name"].apply(parse_sample_name)
    parsed_df = pd.DataFrame(parsed.tolist(), index=samples_df.index)
    samples_df = pd.concat([samples_df, parsed_df], axis=1)

    # Merge study-level metadata into samples
    samples_df = samples_df.merge(
        study_meta_df[["osd_id", "experiment_platform", "hardware", "project_type",
                       "study_factor_name", "flight_program", "mission"]],
        left_on="id.accession", right_on="osd_id", how="left"
    )

    # Use structured metadata where available, fall back to parsed
    samples_df["ecotype_final"] = samples_df["study.characteristics.ecotype"].apply(
        lambda x: x if x and x.lower() != "nan" else None
    ).fillna(samples_df["ecotype_parsed"])

    samples_df["organ_final"] = samples_df.apply(
        lambda r: r["study.characteristics.organism part"] if r["study.characteristics.organism part"] and r["study.characteristics.organism part"].lower() != "nan"
        else (r["study.characteristics.tissue"] if r["study.characteristics.tissue"] and r["study.characteristics.tissue"].lower() != "nan"
              else r["organ_parsed"]), axis=1
    )

    # Normalize organ names
    organ_map = {
        "plant roots": "root", "plant root": "root", "root": "root",
        "plant shoots": "shoot", "plant shoot": "shoot", "shoot": "shoot",
        "plant leaves": "leaf", "leaf": "leaf", "leaves": "leaf",
        "hypocotyl": "hypocotyl", "hypocotyl cell culture": "hypocotyl",
        "seedlings": "whole_seedling", "whole organism": "whole_seedling",
        "whole_seedling": "whole_seedling", "cotyledon": "cotyledon",
        "rosette": "rosette", "unknown": "unknown",
    }
    samples_df["organ_final"] = samples_df["organ_final"].str.lower().map(lambda x: organ_map.get(x, x))

    # Classify omics type per sample
    samples_df["omics_type"] = samples_df["id.assay name"].apply(classify_assay)

    # Step 4: Filter to spaceflight studies with Col-0 or WS
    print("\n[4/5] Filtering to spaceflight studies with Col-0 or WS...")
    spaceflight_mask = samples_df["project_type"].str.contains("paceflight", case=False, na=False)
    ecotype_mask = samples_df["ecotype_final"].isin(["Col-0", "WS", "Col-0+WS", "Wasselewskija", "Wassilewskija", "Wassilewskija ecotype", "Wasselewskija and Col-0"])

    # Normalize WS ecotype names
    ws_names = {"Wasselewskija", "Wassilewskija", "Wassilewskija ecotype", "Wasselewskija and Col-0"}
    samples_df.loc[samples_df["ecotype_final"].isin(ws_names), "ecotype_final"] = "WS"
    samples_df.loc[samples_df["ecotype_final"] == "Col-0+WS", "ecotype_final"] = "Col-0+WS"

    filtered = samples_df[spaceflight_mask & ecotype_mask].copy()
    print(f"  Spaceflight + Col-0/WS samples: {len(filtered)}")
    print(f"  Studies: {filtered['id.accession'].nunique()}")
    print(f"  Ecotypes: {filtered['ecotype_final'].value_counts().to_dict()}")
    print(f"  Hardware: {filtered['hardware'].value_counts().to_dict()}")
    print(f"  Organs: {filtered['organ_final'].value_counts().to_dict()}")
    print(f"  Flight: {filtered['flight_parsed'].value_counts().to_dict()}")
    print(f"  Omics: {filtered['omics_type'].value_counts().to_dict()}")

    # Step 5: Save outputs
    print("\n[5/5] Saving outputs...")

    # Full sample metadata (all Arabidopsis)
    samples_df.to_csv(os.path.join(OUTPUT_DIR, "all_arabidopsis_samples.tsv"), sep="\t", index=False)
    print(f"  Saved all_arabidopsis_samples.tsv ({len(samples_df)} rows)")

    # Filtered sample metadata (spaceflight Col-0/WS)
    filtered.to_csv(os.path.join(OUTPUT_DIR, "sample_metadata.tsv"), sep="\t", index=False)
    print(f"  Saved sample_metadata.tsv ({len(filtered)} rows)")

    # Study metadata
    study_meta_df.to_csv(os.path.join(OUTPUT_DIR, "study_metadata.tsv"), sep="\t", index=False)
    print(f"  Saved study_metadata.tsv ({len(study_meta_df)} rows)")

    # Factorial coverage matrix
    print("\n  Building factorial coverage matrix...")
    if len(filtered) > 0:
        coverage = filtered.groupby(["hardware", "flight_parsed", "organ_final", "ecotype_final", "omics_type"]).size().reset_index(name="n_samples")
        coverage.to_csv(os.path.join(OUTPUT_DIR, "factorial_coverage.tsv"), sep="\t", index=False)
        print(f"  Saved factorial_coverage.tsv ({len(coverage)} cells)")

        # Summary table
        print("\n  Factorial coverage summary:")
        for factor in ["hardware", "flight_parsed", "organ_final", "ecotype_final", "omics_type", "light_parsed"]:
            print(f"    {factor}: {filtered[factor].value_counts().to_dict()}")

    print("\n" + "=" * 70)
    print("OSDR mining complete.")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
