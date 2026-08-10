"""
Extract CFD covariates from LunarLeaf-CFD results tables and map them
to each OSDR sample based on hardware x flight condition.

CFD model scenarios (from LunarLeaf-CFD):
  BRIC   = hermetically sealed (zero-flux walls, enclosure drifts)
  CARA   = gas-permeable tape (semi-permeable lid, enclosure ~ambient)
  VEGGIE = vented + forced airflow (both enclosure & boundary layer controlled)

Hardware mapping (OSDR hardware -> CFD scenario):
  BRIC       -> BRIC  (sealed)
  BRIC-LED   -> BRIC  (sealed, LED lighting - same enclosure physics)
  APEX-TAGES -> BRIC  (sealed canister on Space Shuttle)
  CARA       -> CARA  (semi-permeable)
  Veggie     -> VEGGIE (vented + fan)
  EMCS       -> VEGGIE (controlled atmosphere with gas exchange)

Gravity mapping:
  FLT (Space Flight) -> microgravity (g=0)
  GC  (Ground Control) -> 1g (g=9.81)

Scale mapping (organ -> CFD spatial scale):
  leaf, shoot, hypocotyl, undifferentiated cell culture -> leaf
  whole_seedling -> rosette
  root -> leaf (boundary layer physics applies; no photosynthesis term)
  unknown -> rosette (conservative default for seedlings)

Outputs:
  cfd_covariates_by_scenario.tsv  -- one row per (hardware_cfd, gravity, scale)
  cfd_covariates_by_sample.tsv    -- one row per OSDR sample with CFD scalars
"""
import os
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CFD_TAB = "/workspace/LunarLeaf-CFD/results/tables"
DATA = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/data"
OUT = "/mnt/shared-workspace/microgravity_atmospheric_adaptation/cfd_covariates"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load all CFD tables
# ---------------------------------------------------------------------------
t2  = pd.read_csv(f"{CFD_TAB}/T2_model_sweep.csv")          # gravity sweep: dC_H2O/CO2/O2 mean & peak
t6  = pd.read_csv(f"{CFD_TAB}/T6_forced_airflow.csv")        # forced airflow: fan speed vs dC
t8  = pd.read_csv(f"{CFD_TAB}/T8_enclosure_timescales.csv")  # BRIC enclosure timescales
t9  = pd.read_csv(f"{CFD_TAB}/T9_fan_by_scale.csv")          # fan by scale
t10 = pd.read_csv(f"{CFD_TAB}/T10_hardware_by_scale.csv")    # hardware x scale: surf dC
t11 = pd.read_csv(f"{CFD_TAB}/T11_photosynthesis_feedback.csv")  # 12h carbon % Earth
t12 = pd.read_csv(f"{CFD_TAB}/T12_feedback_spatial.csv")     # net assimilation % by scenario
t13 = pd.read_csv(f"{CFD_TAB}/T13_boundary_layer.csv")       # g_bl, delta, Sherwood by gravity

# ---------------------------------------------------------------------------
# 2. Build per-scenario covariate table
#    Rows = (cfd_hardware, gravity, scale)
# ---------------------------------------------------------------------------
# Hardware scenarios modeled in CFD
HW_CFD = ["BRIC", "CARA", "VEGGIE"]
GRAVITIES = {"ug": 0.0, "earth": 9.81}
SCALES = ["leaf", "rosette", "canopy"]

# -- From T13: boundary layer conductance by gravity x scale --
t13_clean = t13.copy()
t13_clean["scale"] = t13_clean["scenario"].str.split("-").str[0]
t13_clean["gravity_key"] = t13_clean["scenario"].str.split("-").str[1]
t13_clean["gravity_g"] = t13_clean["gravity_g"].astype(float)

# -- From T2: full gas-exchange sweep (dC mean/peak for H2O, CO2, O2) --
t2_clean = t2.copy()
t2_clean["scale"] = t2_clean["scenario"].str.split("-").str[0]
t2_clean["gravity_key"] = t2_clean["scenario"].str.split("-").str[1]

# -- From T10: hardware-specific surface gradients (microgravity) --
# T10 has surf_dC_H2O, surf_dC_CO2, dishmean_CO2 per (scale, hardware)
t10_ug = t10.copy()
t10_ug["gravity_key"] = "ug"  # T10 is microgravity-specific

# -- From T11: 12h carbon fixation (% Earth) by enclosure --
t11_map = {"Earth": ("earth", None), "VEGGIE": ("ug", "VEGGIE"),
           "CARA": ("ug", "CARA"), "BRIC": ("ug", "BRIC")}

# -- From T12: net assimilation % by scenario --
t12_clean = t12.copy()

# -- From T5: physical predictions (CO2/O2 drawdown in ppm) --
t5 = pd.read_csv(f"{CFD_TAB}/T5_physical_prediction.csv")
# Parse scenario column: "Leaf · Earth" -> scale=leaf, gravity=earth
t5_clean = t5.copy()
t5_clean["scale"] = t5_clean["scenario"].str.split("·").str[0].str.strip().str.lower()
t5_clean["gravity_key"] = t5_clean["scenario"].str.split("·").str[1].str.strip().str.lower().str.replace("µg", "ug")

# ---------------------------------------------------------------------------
# 3. Assemble the scenario covariate table
# ---------------------------------------------------------------------------
rows = []

for hw in HW_CFD:
    for gkey, gval in GRAVITIES.items():
        for scale in SCALES:
            row = {
                "cfd_hardware": hw,
                "gravity_key": gkey,
                "gravity_g": gval,
                "scale": scale,
            }

            # --- Boundary layer (T13) ---
            bl = t13_clean[(t13_clean["scale"] == scale) & (t13_clean["gravity_key"] == gkey)]
            if len(bl) > 0:
                bl = bl.iloc[0]
                row["g_bl_mol_m2_s"] = bl["g_bl_mol_m2_s"]
                row["delta_mm"] = bl["delta_mm"]
                row["sherwood"] = bl["Sherwood"]
                row["dC_CO2_mean_bl"] = bl["dC_CO2_mean"]
                row["o2_excess_ppm_bl"] = bl["o2_excess_ppm"]
            else:
                row["g_bl_mol_m2_s"] = np.nan
                row["delta_mm"] = np.nan
                row["sherwood"] = np.nan
                row["dC_CO2_mean_bl"] = np.nan
                row["o2_excess_ppm_bl"] = np.nan

            # --- Full gas-exchange sweep (T2) ---
            ge = t2_clean[(t2_clean["scale"] == scale) & (t2_clean["gravity_key"] == gkey)]
            if len(ge) > 0:
                ge = ge.iloc[0]
                row["u_max"] = ge["u_max"]
                row["Ra_H2O"] = ge["Ra_H2O"]
                row["dC_H2O_mean"] = ge["dC_H2O_mean"]
                row["dC_H2O_peak"] = ge["dC_H2O_peak"]
                row["dC_CO2_mean"] = ge["dC_CO2_mean"]
                row["dC_CO2_peak"] = ge["dC_CO2_peak"]
                row["dC_O2_mean"] = ge["dC_O2_mean"]
                row["dC_O2_peak"] = ge["dC_O2_peak"]
            else:
                for c in ["u_max", "Ra_H2O", "dC_H2O_mean", "dC_H2O_peak",
                          "dC_CO2_mean", "dC_CO2_peak", "dC_O2_mean", "dC_O2_peak"]:
                    row[c] = np.nan

            # --- Hardware-specific surface gradients (T10, microgravity only) ---
            if gkey == "ug":
                hw10 = t10_ug[(t10_ug["scale"] == scale) & (t10_ug["hardware"] == hw)]
                if len(hw10) > 0:
                    hw10 = hw10.iloc[0]
                    row["surf_dC_H2O_hw"] = hw10["surf_dC_H2O"]
                    row["surf_dC_CO2_hw"] = hw10["surf_dC_CO2"]
                    row["dishmean_CO2_hw"] = hw10["dishmean_CO2"]
                else:
                    row["surf_dC_H2O_hw"] = np.nan
                    row["surf_dC_CO2_hw"] = np.nan
                    row["dishmean_CO2_hw"] = np.nan
            else:
                # Ground control: hardware doesn't matter for boundary layer
                # (free convection dominates); use the generic earth values
                row["surf_dC_H2O_hw"] = row.get("dC_H2O_mean", np.nan)
                row["surf_dC_CO2_hw"] = row.get("dC_CO2_mean", np.nan)
                row["dishmean_CO2_hw"] = 0.0  # vented/steady at ambient on Earth

            # --- 12h carbon fixation % (T11) ---
            if gkey == "earth":
                row["carbon_12h_pct_earth"] = 100.0
            else:
                t11_row = t11[t11["enclosure"] == hw]
                if len(t11_row) > 0:
                    row["carbon_12h_pct_earth"] = t11_row.iloc[0]["12h carbon (% Earth)"]
                    row["A_start_umol_m2_s"] = t11_row.iloc[0]["A start (µmol/m2/s)"]
                    row["A_end_umol_m2_s"] = t11_row.iloc[0]["A end (µmol/m2/s)"]
                else:
                    row["carbon_12h_pct_earth"] = np.nan
                    row["A_start_umol_m2_s"] = np.nan
                    row["A_end_umol_m2_s"] = np.nan

            # --- Net assimilation % (T12) ---
            if gkey == "earth":
                t12_row = t12_clean[t12_clean["scenario"].str.lower() == f"{scale} earth"]
            else:
                t12_row = t12_clean[t12_clean["scenario"].str.lower() == f"{scale} ug {hw.lower().replace('veggie','veggie')}"]
            if len(t12_row) > 0:
                row["net_assimilation_pct"] = t12_row.iloc[0]["net_assimilation_pct"]
            else:
                # Fallback: try generic ug scenario
                if gkey == "ug":
                    t12_fb = t12_clean[t12_clean["scenario"].str.lower() == f"{scale} ug"]
                    if len(t12_fb) > 0:
                        row["net_assimilation_pct"] = t12_fb.iloc[0]["net_assimilation_pct"]
                    else:
                        row["net_assimilation_pct"] = np.nan
                else:
                    row["net_assimilation_pct"] = np.nan

            # --- Physical predictions: CO2/O2 drawdown (T5) ---
            pp = t5_clean[(t5_clean["scale"] == scale) & (t5_clean["gravity_key"] == gkey)]
            if len(pp) > 0:
                pp = pp.iloc[0]
                row["co2_drawdown_mean_ppm"] = pp["CO2 drawdown mean (ppm)"]
                row["co2_drawdown_peak_ppm"] = pp["CO2 drawdown peak (ppm)"]
                row["o2_buildup_mean_ppm"] = pp["O2 build-up mean (ppm)"]
            else:
                row["co2_drawdown_mean_ppm"] = np.nan
                row["co2_drawdown_peak_ppm"] = np.nan
                row["o2_buildup_mean_ppm"] = np.nan

            # --- Enclosure classification ---
            enclosure_map = {"BRIC": "sealed", "CARA": "semi-permeable", "VEGGIE": "vented+fan"}
            row["enclosure_type"] = enclosure_map[hw]

            # --- Enclosure CO2 depletion timescale (T8, BRIC only) ---
            if hw == "BRIC" and gkey == "ug":
                row["co2_depletion_min"] = 7.0       # minutes (T8)
                row["co2_stress_hr"] = 9.4           # hours (T8)
                row["o2_hypoxia_days"] = 6.5         # days (T8)
            elif hw == "BRIC" and gkey == "earth":
                # On Earth, free convection mixes the enclosure; but BRIC is sealed
                # The depletion still occurs but slower due to convection
                row["co2_depletion_min"] = np.nan
                row["co2_stress_hr"] = np.nan
                row["o2_hypoxia_days"] = np.nan
            else:
                # CARA/VEGGIE vent to cabin - bounded
                row["co2_depletion_min"] = np.nan
                row["co2_stress_hr"] = np.nan
                row["o2_hypoxia_days"] = np.nan

            # --- Forced airflow equivalent (T6/T9) ---
            # Earth-equivalent fan speed that nulls the microgravity penalty
            fan_equiv = {"leaf": 2.8, "rosette": 11.0, "canopy": np.nan}  # cm/s
            row["earth_equiv_fan_cm_s"] = fan_equiv.get(scale, np.nan)
            # VEGGIE has forced airflow (~2.8 cm/s or more)
            if hw == "VEGGIE" and gkey == "ug":
                row["has_forced_airflow"] = 1
                row["fan_cm_s"] = 2.8  # approximate VEGGIE fan speed
            else:
                row["has_forced_airflow"] = 0
                row["fan_cm_s"] = 0.0

            rows.append(row)

scenario_df = pd.DataFrame(rows)
scenario_df.to_csv(f"{OUT}/cfd_covariates_by_scenario.tsv", sep="\t", index=False)
print(f"Scenario table: {scenario_df.shape}")
print(scenario_df[["cfd_hardware", "gravity_key", "scale", "g_bl_mol_m2_s",
                    "delta_mm", "carbon_12h_pct_earth", "enclosure_type"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 4. Map CFD covariates to each OSDR sample
# ---------------------------------------------------------------------------
meta = pd.read_csv(f"{DATA}/sample_metadata_final.tsv", sep="\t")

# Hardware -> CFD scenario mapping
HW_MAP = {
    "BRIC": "BRIC",
    "BRIC-LED": "BRIC",
    "APEX-TAGES": "BRIC",
    "CARA": "CARA",
    "Veggie": "VEGGIE",
    "EMCS": "VEGGIE",
}

# Organ -> CFD scale mapping
SCALE_MAP = {
    "leaf": "leaf",
    "shoot": "leaf",
    "hypocotyl": "leaf",
    "undifferentiated cell culture": "leaf",
    "whole_seedling": "rosette",
    "root": "leaf",
    "unknown": "rosette",
}

# Flight -> gravity mapping
def gravity_key(flight):
    if pd.isna(flight):
        return "earth"
    return "ug" if flight == "FLT" else "earth"

# Build lookup from scenario_df
def get_cfd_row(hw, flight, organ):
    cfd_hw = HW_MAP.get(hw, "BRIC")
    gkey = gravity_key(flight)
    scale = SCALE_MAP.get(organ, "rosette")
    match = scenario_df[(scenario_df["cfd_hardware"] == cfd_hw) &
                        (scenario_df["gravity_key"] == gkey) &
                        (scenario_df["scale"] == scale)]
    if len(match) == 0:
        # Fallback to rosette scale
        match = scenario_df[(scenario_df["cfd_hardware"] == cfd_hw) &
                            (scenario_df["gravity_key"] == gkey) &
                            (scenario_df["scale"] == "rosette")]
    if len(match) == 0:
        return pd.Series(dtype=float)
    return match.iloc[0]

# Apply mapping
cfd_cols = [c for c in scenario_df.columns if c not in ["cfd_hardware", "gravity_key", "scale"]]
sample_cfd = []
for _, srow in meta.iterrows():
    cfd = get_cfd_row(srow["hardware"], srow["flight_final"], srow["organ_final"])
    out = {
        "sample_name": srow["id.sample name"],
        "osd_id": srow["osd_id"],
        "hardware": srow["hardware"],
        "cfd_hardware": HW_MAP.get(srow["hardware"], "BRIC"),
        "flight_final": srow["flight_final"],
        "gravity_key": gravity_key(srow["flight_final"]),
        "organ_final": srow["organ_final"],
        "cfd_scale": SCALE_MAP.get(srow["organ_final"], "rosette"),
    }
    for c in cfd_cols:
        out[c] = cfd.get(c, np.nan)
    sample_cfd.append(out)

sample_cfd_df = pd.DataFrame(sample_cfd)
sample_cfd_df.to_csv(f"{OUT}/cfd_covariates_by_sample.tsv", sep="\t", index=False)
print(f"\nSample-level CFD covariates: {sample_cfd_df.shape}")
print(sample_cfd_df[["sample_name", "cfd_hardware", "gravity_key", "cfd_scale",
                      "g_bl_mol_m2_s", "delta_mm", "carbon_12h_pct_earth"]].head(20).to_string(index=False))

# ---------------------------------------------------------------------------
# 5. Summary: covariate coverage
# ---------------------------------------------------------------------------
print("\n=== CFD covariate coverage by hardware x flight ===")
print(sample_cfd_df.groupby(["cfd_hardware", "gravity_key"])[["g_bl_mol_m2_s",
    "delta_mm", "carbon_12h_pct_earth", "net_assimilation_pct"]].mean().to_string())

print("\n=== Key CFD scalars by scenario (rosette scale) ===")
rosette = scenario_df[scenario_df["scale"] == "rosette"]
print(rosette[["cfd_hardware", "gravity_key", "g_bl_mol_m2_s", "delta_mm",
               "dC_CO2_mean", "carbon_12h_pct_earth", "enclosure_type"]].to_string(index=False))

print(f"\nSaved: {OUT}/cfd_covariates_by_scenario.tsv")
print(f"Saved: {OUT}/cfd_covariates_by_sample.tsv")
