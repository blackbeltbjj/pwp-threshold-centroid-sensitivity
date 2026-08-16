#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
OSAF PROGRAM 35 — PWP LCC-ONLY FULL-RECORD ROBUSTNESS
VERSION 1.0.0
===============================================================================
PURPOSE
    Close the Paper-1 LCC-only science-freeze experiment WITHOUT rereading or
    redefining daily OISST. This program consumes the canonical full-record
    Program-31 daily connectivity table and treats the largest-connected-
    component (LCC) centroid as the alternative observable.

CANONICAL INPUT
    outputs/tables/threshold_comparison/pwp_long_term_connectivity/
        pwp_daily_connectivity_diagnostics.csv

EXPECTED INPUT COLUMNS
    date, threshold_c, largest_component_area_km2,
    largest_centroid_lon_360, largest_centroid_lat,
    plus Program-31 full-domain diagnostics.

OUTPUTS
    outputs/tables/threshold_comparison/pwp_lcc_only_robustness/
        pwp_lcc_only_threshold_summary.csv
        pwp_lcc_only_interthreshold_separation_daily.csv
        pwp_lcc_only_interthreshold_separation_summary.csv
        pwp_lcc_only_28_to_29_change_summary.csv

    outputs/reports/threshold_comparison/pwp_lcc_only_robustness/
        PROGRAM35_PWP_LCC_ONLY_FULL_RECORD_ROBUSTNESS.txt

SCIENTIFIC RULES
    - thresholds: 28.0, 28.5, 29.0 °C
    - common dates only
    - equal-day spherical mean centroid (3-D unit-vector mean)
    - continuous/unwrapped Pacific longitude branch for longitudinal percentiles
    - radial excursion = great-circle distance from the threshold-specific
      long-term spherical mean LCC centroid
    - same-day inter-threshold separation = great-circle distance
    - no p-values are computed across the three thresholds
===============================================================================
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = (28.0, 28.5, 29.0)
EARTH_RADIUS_KM = 6371.0088

INPUT = ROOT / "outputs" / "tables" / "threshold_comparison" / "pwp_long_term_connectivity" / "pwp_daily_connectivity_diagnostics.csv"
TABLE_DIR = ROOT / "outputs" / "tables" / "threshold_comparison" / "pwp_lcc_only_robustness"
REPORT_DIR = ROOT / "outputs" / "reports" / "threshold_comparison" / "pwp_lcc_only_robustness"

def spherical_mean(lon_deg, lat_deg):
    lon = np.deg2rad(np.asarray(lon_deg, float))
    lat = np.deg2rad(np.asarray(lat_deg, float))
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    x = np.mean(np.cos(lat) * np.cos(lon))
    y = np.mean(np.cos(lat) * np.sin(lon))
    z = np.mean(np.sin(lat))
    norm = np.sqrt(x*x + y*y + z*z)
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("Degenerate spherical mean vector.")
    x, y, z = x/norm, y/norm, z/norm
    mlon = np.rad2deg(np.arctan2(y, x)) % 360.0
    mlat = np.rad2deg(np.arctan2(z, np.sqrt(x*x+y*y)))
    return float(mlon), float(mlat)

def gc_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    dlon = lon2-lon1
    dlat = lat2-lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return EARTH_RADIUS_KM * 2*np.arcsin(np.sqrt(np.clip(a,0,1)))

def unwrap_about(lon_deg, center_deg):
    lon = np.asarray(lon_deg, float)
    return center_deg + ((lon-center_deg+180.0) % 360.0 - 180.0)

def main():
    print("="*78)
    print("OSAF PROGRAM 35 — PWP LCC-ONLY FULL-RECORD ROBUSTNESS")
    print("="*78)
    if not INPUT.is_file():
        raise FileNotFoundError(f"Canonical Program-31 table not found:\n{INPUT}")
    df = pd.read_csv(INPUT)
    required = {
        "date","threshold_c","largest_component_area_km2",
        "largest_centroid_lon_360","largest_centroid_lat"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    df["date"] = pd.to_datetime(df["date"]).dt.normalize()
    df["threshold_c"] = df["threshold_c"].astype(float)

    # exact common dates
    sets = [set(df.loc[np.isclose(df.threshold_c,t),"date"]) for t in THRESHOLDS]
    common = sorted(set.intersection(*sets))
    if not common:
        raise ValueError("No common dates across thresholds.")
    df = df[df.date.isin(common)].copy()
    expected = len(common) * len(THRESHOLDS)
    if len(df) != expected:
        raise RuntimeError(f"Expected {expected:,} rows after common-date filtering; found {len(df):,}.")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    summaries=[]
    daily_parts=[]
    for t in THRESHOLDS:
        g=df[np.isclose(df.threshold_c,t)].sort_values("date").copy()
        lon=g["largest_centroid_lon_360"].to_numpy(float)
        lat=g["largest_centroid_lat"].to_numpy(float)
        area=g["largest_component_area_km2"].to_numpy(float)
        mlon,mlat=spherical_mean(lon,lat)
        lon_u=unwrap_about(lon,mlon)
        radial=gc_km(mlon,mlat,lon,lat)
        g["lcc_lon_continuous_degE"]=lon_u
        g["lcc_radial_excursion_km"]=radial
        daily_parts.append(g)
        p05,p95=np.quantile(lon_u,[.05,.95])
        summaries.append({
            "threshold_c":t, "n_days":len(g),
            "start_date":g.date.min().date().isoformat(),
            "end_date":g.date.max().date().isoformat(),
            "mean_spherical_lcc_lon_360":mlon,
            "mean_spherical_lcc_lat":mlat,
            "median_lcc_lon_continuous_degE":float(np.median(lon_u)),
            "lcc_lon_p05_degE":float(p05),
            "lcc_lon_p95_degE":float(p95),
            "lcc_lon_p95_minus_p05_deg":float(p95-p05),
            "lcc_lon_sd_deg":float(np.std(lon_u,ddof=1)),
            "lcc_radial_p95_km":float(np.quantile(radial,.95)),
            "lcc_radial_max_km":float(np.max(radial)),
            "mean_lcc_area_million_km2":float(np.mean(area)/1e6),
            "median_lcc_area_million_km2":float(np.median(area)/1e6),
        })
    summary=pd.DataFrame(summaries)
    summary.to_csv(TABLE_DIR/"pwp_lcc_only_threshold_summary.csv",index=False,float_format="%.8f")

    # Same-day inter-threshold LCC separations
    wide={}
    for t in THRESHOLDS:
        tag=str(t).replace(".","p")
        g=df[np.isclose(df.threshold_c,t)].set_index("date")
        wide[t]=g
    sep_rows=[]
    pairs=((28.0,28.5),(28.0,29.0),(28.5,29.0))
    for date in pd.DatetimeIndex(common):
        row={"date":date.date().isoformat()}
        for a,b in pairs:
            ga,gb=wide[a].loc[date],wide[b].loc[date]
            dist=float(gc_km(
                ga["largest_centroid_lon_360"],ga["largest_centroid_lat"],
                gb["largest_centroid_lon_360"],gb["largest_centroid_lat"]))
            dlon=float((gb["largest_centroid_lon_360"]-ga["largest_centroid_lon_360"]+180)%360-180)
            dlat=float(gb["largest_centroid_lat"]-ga["largest_centroid_lat"])
            key=f"{a:g}_vs_{b:g}".replace(".","p")
            row[f"{key}_distance_km"]=dist
            row[f"{key}_delta_lon_b_minus_a_deg"]=dlon
            row[f"{key}_delta_lat_b_minus_a_deg"]=dlat
        sep_rows.append(row)
    sep=pd.DataFrame(sep_rows)
    sep.to_csv(TABLE_DIR/"pwp_lcc_only_interthreshold_separation_daily.csv",index=False,float_format="%.8f")

    ss=[]
    for a,b in pairs:
        key=f"{a:g}_vs_{b:g}".replace(".","p")
        d=sep[f"{key}_distance_km"]
        dlo=sep[f"{key}_delta_lon_b_minus_a_deg"]
        dla=sep[f"{key}_delta_lat_b_minus_a_deg"]
        ss.append({
            "threshold_a_c":a,"threshold_b_c":b,"n_days":len(sep),
            "distance_median_km":float(d.median()),
            "distance_p95_km":float(d.quantile(.95)),
            "distance_max_km":float(d.max()),
            "median_delta_lon_b_minus_a_deg":float(dlo.median()),
            "median_delta_lat_b_minus_a_deg":float(dla.median()),
        })
    sep_summary=pd.DataFrame(ss)
    sep_summary.to_csv(TABLE_DIR/"pwp_lcc_only_interthreshold_separation_summary.csv",index=False,float_format="%.8f")

    s28=summary.loc[np.isclose(summary.threshold_c,28.0)].iloc[0]
    s29=summary.loc[np.isclose(summary.threshold_c,29.0)].iloc[0]
    changes=[]
    for metric in ["mean_lcc_area_million_km2","lcc_lon_p95_minus_p05_deg","lcc_lon_sd_deg","lcc_radial_p95_km"]:
        a=float(s28[metric]); b=float(s29[metric])
        changes.append({"metric":metric,"value_28C":a,"value_29C":b,
                        "absolute_change_29_minus_28":b-a,
                        "percent_change_28_to_29":100.0*(b-a)/a})
    change=pd.DataFrame(changes)
    change.to_csv(TABLE_DIR/"pwp_lcc_only_28_to_29_change_summary.csv",index=False,float_format="%.8f")

    report=[]
    report += ["OSAF PROGRAM 35 — PWP LCC-ONLY FULL-RECORD ROBUSTNESS","",
               f"Canonical input: {INPUT}",f"Common dates: {len(common):,} | {common[0]:%Y-%m-%d} to {common[-1]:%Y-%m-%d}","",
               "THRESHOLD-SPECIFIC LCC OBSERVABLES"]
    for _,r in summary.iterrows():
        report.append(
            f"{r.threshold_c:4.1f} °C | N={int(r.n_days):,} | "
            f"mean spherical LCC centroid=({r.mean_spherical_lcc_lon_360:.3f}°E, {r.mean_spherical_lcc_lat:+.3f}°) | "
            f"lon P05-P95={r.lcc_lon_p05_degE:.3f}–{r.lcc_lon_p95_degE:.3f}°E "
            f"(range={r.lcc_lon_p95_minus_p05_deg:.3f}°) | radial P95={r.lcc_radial_p95_km:.1f} km"
        )
    report += ["","SAME-DAY INTER-THRESHOLD LCC SEPARATIONS"]
    for _,r in sep_summary.iterrows():
        report.append(
            f"{r.threshold_a_c:.1f} vs {r.threshold_b_c:.1f} °C | N={int(r.n_days):,} | "
            f"median={r.distance_median_km:.1f} km | P95={r.distance_p95_km:.1f} km | max={r.distance_max_km:.1f} km | "
            f"median Δlon(b-a)={r.median_delta_lon_b_minus_a_deg:+.2f}° | median Δlat(b-a)={r.median_delta_lat_b_minus_a_deg:+.2f}°"
        )
    report += ["","28→29 °C CHANGES"]
    for _,r in change.iterrows():
        report.append(f"{r.metric}: {r.value_28C:.6f} -> {r.value_29C:.6f} | {r.percent_change_28_to_29:+.2f}%")
    report += ["","INTERPRETATION RULE:",
               "If the monotonic threshold effects persist for the LCC-only observable, detached warm patches cannot be the sole explanation.",
               "If they weaken materially, fragmentation contributes to the full-domain threshold sensitivity.",
               "This remains a geometric/statistical sensitivity test, not evidence of ocean-dynamical inertia."]
    (REPORT_DIR/"PROGRAM35_PWP_LCC_ONLY_FULL_RECORD_ROBUSTNESS.txt").write_text("\n".join(report),encoding="utf-8")
    print("\n".join(report))
    print("\nPROGRAM 35 COMPLETED SUCCESSFULLY.")

if __name__ == "__main__":
    main()
