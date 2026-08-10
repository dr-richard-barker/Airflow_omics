#!/usr/bin/env python3
"""
Download processed data files from OSDR for all included spaceflight studies.
Focuses on normalized counts, unnormalized counts, differential expression,
and sample tables for each omics type.
"""
import urllib.request, json, os, time, re
import pandas as pd

API_BASE = "https://visualization.osdr.nasa.gov/biodata/api"
DATA_DIR = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
RAW_DIR = os.path.join(DATA_DIR, "raw_osdr")
os.makedirs(RAW_DIR, exist_ok=True)

# Studies to download (spaceflight Col-0/WS)
STUDIES = {
    "OSD-7": ["array_normalized_expression_probeset", "array_normalized_intensities_probe", "array_microarray_v0_runsheet"],
    "OSD-37": ["rna_seq_Normalized_Counts_GLbulkRNAseq", "rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_differential_expression_GLbulkRNAseq", "rna_seq_SampleTable_GLbulkRNAseq"],
    "OSD-38": ["rna_seq_Normalized_Counts", "rna_seq_Unnormalized_Counts", "rna_seq_differential_expression", "rna_seq_SampleTable", "membrane_proteome_WyattLab_membrane", "soluble_proteome_WyattLab_soluble"],
    "OSD-120": ["rna_seq_Normalized_Counts_GLbulkRNAseq", "rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_differential_expression_rRNArm_GLbulkRNAseq", "rna_seq_SampleTable_GLbulkRNAseq"],
    "OSD-147": ["array_normalized_expression_probeset", "array_microarray_v0_runsheet"],
    "OSD-17": ["array_normalized_expression_probeset", "array_microarray_v0_runsheet"],
    "OSD-205": ["array_normalized_expression_probeset", "array_microarray_v0_runsheet"],
    "OSD-217": ["rna_seq_Normalized_Counts", "rna_seq_Unnormalized_Counts", "rna_seq_differential_expression", "rna_seq_SampleTable", "rna_seq_contrasts"],
    "OSD-321": ["rna_seq_Normalized_Counts", "rna_seq_RSEM_Unnormalized_Counts", "rna_seq_STAR_Unnormalized_Counts", "rna_seq_differential_expression", "rna_seq_SampleTable"],
    "OSD-416": [],  # No processed files found
    "OSD-522": ["rna_seq_Normalized_Counts_GLbulkRNAseq", "rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_differential_expression_GLbulkRNAseq", "rna_seq_SampleTable_GLbulkRNAseq", "proteomics_GO_Shoot_MEM_Report", "proteomics_GO_Root_MEM_Report", "proteomics_GO_Shoot_SOL_Report", "proteomics_GO_Root_SOL_Report"],
    "OSD-615": ["glycomics_Raw-APEX-03-01_SUBMITTED", "glycomics_APEX-03-01_TRANSFORMED"],
    "OSD-625": [],
    "OSD-314": ["rna_seq_Normalized_Counts", "rna_seq_RSEM_Unnormalized_Counts", "rna_seq_STAR_Unnormalized_Counts", "rna_seq_differential_expression_GLbulkRNAseq", "rna_seq_SampleTable_GLbulkRNAseq"],
    "OSD-678": ["rna_seq_Normalized_Counts_GLbulkRNAseq", "rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq", "rna_seq_differential_expression_GLbulkRNAseq", "rna_seq_SampleTable_GLbulkRNAseq"],
}


def get_file_list(osd):
    """Get file listing for a study."""
    url = f"{API_BASE}/v2/dataset/{osd}/files/"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        resp = urllib.request.urlopen(req, timeout=30)
        data = json.loads(resp.read())
        files_info = data.get(osd, {}).get("files", {})
        if isinstance(files_info, dict):
            return files_info
        return {}
    except Exception as e:
        print(f"  ERROR getting file list for {osd}: {e}")
        return {}


def download_file(filename, file_info, dest_dir):
    """Download a single file using the URL field from file info."""
    # The file info dict has a 'URL' field with the direct download link
    if isinstance(file_info, dict):
        dl_url = file_info.get("URL", "") or file_info.get("url", "")
    else:
        dl_url = ""

    if not dl_url:
        return False

    try:
        dest_path = os.path.join(dest_dir, filename)
        req = urllib.request.Request(dl_url)
        resp = urllib.request.urlopen(req, timeout=180)
        with open(dest_path, 'wb') as f:
            f.write(resp.read())
        size = os.path.getsize(dest_path)
        return size > 0
    except Exception as e:
        print(f"    Download error for {filename}: {e}")
        return False


def main():
    print("=" * 70)
    print("Downloading processed data files from OSDR")
    print("=" * 70)

    total_downloaded = 0
    download_log = []

    for osd, patterns in STUDIES.items():
        if not patterns:
            print(f"\n{osd}: No processed files to download (skipping)")
            continue

        print(f"\n{osd}: Looking for {len(patterns)} file patterns...")
        files = get_file_list(osd)
        if not files:
            print(f"  No files found")
            continue

        osd_dir = os.path.join(RAW_DIR, osd)
        os.makedirs(osd_dir, exist_ok=True)

        downloaded = 0
        for fname, finfo in files.items():
            # Check if filename matches any of our patterns
            fl = fname.lower()
            if any(p.lower() in fl for p in patterns):
                dest = os.path.join(osd_dir, fname)
                if os.path.exists(dest) and os.path.getsize(dest) > 0:
                    print(f"  EXISTS: {fname} ({os.path.getsize(dest)} bytes)")
                    downloaded += 1
                    download_log.append({"osd": osd, "filename": fname, "status": "exists", "size": os.path.getsize(dest)})
                    continue
                print(f"  Downloading: {fname}...", end=" ")
                success = download_file(fname, finfo, osd_dir)
                if success:
                    size = os.path.getsize(dest)
                    print(f"OK ({size} bytes)")
                    downloaded += 1
                    download_log.append({"osd": osd, "filename": fname, "status": "downloaded", "size": size})
                else:
                    print("FAILED")
                    download_log.append({"osd": osd, "filename": fname, "status": "failed", "size": 0})
                time.sleep(0.5)  # Rate limit

        print(f"  {osd}: {downloaded}/{len(patterns)} files downloaded")
        total_downloaded += downloaded

    # Save download log
    log_df = pd.DataFrame(download_log)
    log_path = os.path.join(DATA_DIR, "download_log.tsv")
    log_df.to_csv(log_path, sep="\t", index=False)

    print(f"\n{'='*70}")
    print(f"Download complete: {total_downloaded} files total")
    print(f"Download log: {log_path}")
    print(f"Raw data: {RAW_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
