#!/usr/bin/env python3
"""
Final metadata cleanup: fix remaining organ/hardware gaps using study-specific
knowledge and sample-name patterns.
"""
import pandas as pd, os, re

OUTPUT_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"

df = pd.read_csv(os.path.join(OUTPUT_DIR, "sample_metadata_enriched.tsv"), sep="\t")
print(f"Loaded {len(df)} samples")

# Fix hardware: OSD-120 is CARA, OSD-678 is Veggie
df.loc[df["id.accession"] == "OSD-120", "hardware"] = "CARA"
df.loc[df["id.accession"] == "OSD-678", "hardware"] = "Veggie"

# Fix organ from sample names for remaining nan/unknown
def fix_organ(row):
    organ = row.get("organ_final", "")
    if pd.notna(organ) and str(organ).lower() not in ("nan", "unknown", ""):
        return organ
    sname = str(row.get("id.sample name", "")).lower()
    osd = row.get("id.accession", "")
    # Study-specific
    if osd == "OSD-120":
        return "root"  # CARA study is all roots
    if osd == "OSD-678":
        return "leaf"  # Light study is leaves
    # Pattern-based
    if "root" in sname:
        return "root"
    if "hypocotyl" in sname or "hypocotylcc" in sname:
        return "hypocotyl"
    if "shoot" in sname:
        return "shoot"
    if "leaf" in sname or "leaves" in sname:
        return "leaf"
    if "sl-pool" in sname or "sl_pool" in sname or "_sl_" in sname or "seedling" in sname:
        return "whole_seedling"
    if "_wo_" in sname or "whole" in sname:
        return "whole_seedling"
    if "undifferentiated" in sname or "cell culture" in str(row.get("organ_factor","")).lower():
        return "undifferentiated_cell_culture"
    # OSD-321 (BRIC, bzip28 mutant) - seedlings
    if osd == "OSD-321":
        return "whole_seedling"
    # OSD-37, OSD-38 - seedlings
    if osd in ("OSD-37", "OSD-38"):
        return "whole_seedling"
    return "unknown"

df["organ_final"] = df.apply(fix_organ, axis=1)

# Fix ecotype: OSD-37 has 4 ecotypes including WS
# Check sample names for WS in OSD-37
def fix_ecotype(row):
    eco = row.get("ecotype_final", "")
    if pd.notna(eco) and str(eco).lower() not in ("nan", "unknown", ""):
        return eco
    sname = str(row.get("id.sample name", "")).lower()
    if "ws-0" in sname or "ws0" in sname or "ws_" in sname:
        return "WS"
    if "col-0" in sname or "col0" in sname or "wt-col" in sname or "col_0" in sname:
        return "Col-0"
    if "wassilewskij" in sname or "wasselewskij" in sname:
        return "WS"
    # Default for studies known to be Col-0
    osd = row.get("id.accession", "")
    if osd in ("OSD-120", "OSD-314", "OSD-321", "OSD-416", "OSD-427", "OSD-44", "OSD-147", "OSD-205", "OSD-522", "OSD-615", "OSD-625", "OSD-678"):
        return "Col-0"
    return eco if pd.notna(eco) else "unknown"

df["ecotype_final"] = df.apply(fix_ecotype, axis=1)

# Fix light: default to "light" for spaceflight studies unless dark specified
def fix_light(row):
    light = row.get("light_final", "")
    if pd.notna(light) and str(light).lower() not in ("nan", "unspecified", ""):
        return light
    sname = str(row.get("id.sample name", "")).lower()
    if "dark" in sname or "blind" in sname:
        return "dark"
    if "alight" in sname or "a_light" in sname or "_light" in sname:
        return "light"
    # Most spaceflight studies grow in light
    return "light"

df["light_final"] = df.apply(fix_light, axis=1)

# Save final metadata
df.to_csv(os.path.join(OUTPUT_DIR, "sample_metadata_final.tsv"), sep="\t", index=False)
print(f"\nSaved sample_metadata_final.tsv ({len(df)} rows)")

# Print final coverage
print("\n=== FINAL factorial coverage ===")
for factor in ["hardware", "flight_final", "organ_final", "ecotype_final", "omics_type", "light_final"]:
    vc = df[factor].value_counts(dropna=False)
    print(f"  {factor}: {vc.to_dict()}")

# Key cross-tabs
print("\n=== Hardware × Flight × Omics ===")
ct = df.groupby(["hardware", "flight_final", "omics_type"]).size().reset_index(name="n")
ct = ct.sort_values(["hardware", "flight_final", "omics_type"])
print(ct.to_string(index=False))

print("\n=== Hardware × Organ ===")
ct2 = df.groupby(["hardware", "organ_final"]).size().reset_index(name="n")
ct2 = ct2.sort_values(["hardware", "n"], ascending=[True, False])
print(ct2.to_string(index=False))

print("\n=== Ecotype × Hardware ===")
ct3 = df.groupby(["ecotype_final", "hardware"]).size().reset_index(name="n")
ct3 = ct3.sort_values(["ecotype_final", "n"], ascending=[True, False])
print(ct3.to_string(index=False))
