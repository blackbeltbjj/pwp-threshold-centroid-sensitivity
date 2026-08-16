#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
OSAF PROGRAM 36 — DAILY PWP OCCURRENCE / PERSISTENCE CLIMATOLOGY
VERSION 1.0.0
===============================================================================
PURPOSE
    For each OISST grid cell and threshold T, calculate the fraction of VALID
    daily observations for which the cell belongs to the threshold-defined PWP:

        F_i(T) = sum_t V_i(t) I[SST_i(t) >= T] / sum_t V_i(t)

    where V_i(t)=1 only when the daily SST value is finite and the cell belongs
    to the fixed Pacific mask.

IMPORTANT
    This program thresholds EACH DAILY SST FIELD. It does NOT threshold the
    temporally averaged SST field.

CANONICAL UPSTREAM METHODS
    src_a/17_methodological_domain_figure.py
    data/processed/pacific_mask_oisst.npy
    data/processed/grid_lat.npy
    data/processed/grid_lon.npy
    daily NOAA OISST under data/raw/

TEMPORAL DOMAIN
    Exact common Program-05 centroid dates for 28.0, 28.5, 29.0 °C.
    Expected frozen Paper-1 record: 1981-09-01 through 2026-07-29, N=16,403.

OUTPUTS
    data/processed/pwp_occurrence_persistence/
        pwp_daily_occurrence_persistence_1981-09-01_2026-07-29.npz

    outputs/tables/threshold_comparison/pwp_occurrence_persistence/
        pwp_occurrence_persistence_summary.csv

    outputs/figures/threshold_comparison/pwp_occurrence_persistence/
        pwp_daily_occurrence_persistence_thresholds.png
        pwp_daily_occurrence_persistence_thresholds.pdf

    outputs/reports/threshold_comparison/pwp_occurrence_persistence/
        PROGRAM36_PWP_DAILY_OCCURRENCE_PERSISTENCE.txt

CHECKPOINT
    Counts are checkpointed periodically as compressed NPZ so a long run can
    resume safely.
===============================================================================
"""
from pathlib import Path
import importlib.util, sys, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
P17_FILE=ROOT/"src_a"/"17_methodological_domain_figure.py"
THRESHOLDS=(28.0,28.5,29.0)
CHECKPOINT_EVERY=100

PROCESSED=ROOT/"data"/"processed"/"pwp_occurrence_persistence"
TABLE_DIR=ROOT/"outputs"/"tables"/"threshold_comparison"/"pwp_occurrence_persistence"
FIG_DIR=ROOT/"outputs"/"figures"/"threshold_comparison"/"pwp_occurrence_persistence"
REPORT_DIR=ROOT/"outputs"/"reports"/"threshold_comparison"/"pwp_occurrence_persistence"
CHECKPOINT=PROCESSED/"pwp_occurrence_checkpoint.npz"

def load_p17():
    if not P17_FILE.is_file():
        raise FileNotFoundError(P17_FILE)
    spec=importlib.util.spec_from_file_location("pwp_program17",P17_FILE)
    mod=importlib.util.module_from_spec(spec); sys.modules[spec.name]=mod
    spec.loader.exec_module(mod)
    return mod
P17=load_p17()

def centroid_path(t):
    folder=str(int(t)) if float(t).is_integer() else f"{t:g}"
    candidates=[
        ROOT/"data"/"processed"/folder/"centroid"/"pwp_centroid_series.csv",
        ROOT/"data"/"processed"/folder/"pwp_centroid_series.csv",
    ]
    for p in candidates:
        if p.is_file(): return p
    raise FileNotFoundError("\n".join(map(str,candidates)))

def common_dates():
    sets=[]
    for t in THRESHOLDS:
        f=pd.read_csv(centroid_path(t))
        col=next((c for c in ("date","Date","DATE","time","Time","TIME","datetime","timestamp") if c in f.columns),None)
        if col is None: raise ValueError(f"No date column for {t} °C.")
        d=pd.to_datetime(f[col],errors="coerce").dropna().dt.normalize()
        sets.append(set(d))
    return pd.DatetimeIndex(sorted(set.intersection(*sets)))

def load_checkpoint(shape):
    if not CHECKPOINT.is_file():
        return 0, np.zeros(shape,np.uint32), {t:np.zeros(shape,np.uint32) for t in THRESHOLDS}
    z=np.load(CHECKPOINT)
    next_index=int(z["next_index"])
    valid=z["valid_count"].astype(np.uint32)
    hits={t:z[f"hits_{str(t).replace('.','p')}"].astype(np.uint32) for t in THRESHOLDS}
    if valid.shape != shape: raise ValueError("Checkpoint grid shape mismatch.")
    return next_index,valid,hits

def save_checkpoint(next_index,valid,hits):
    PROCESSED.mkdir(parents=True,exist_ok=True)
    payload={"next_index":np.array(next_index,dtype=np.int64),"valid_count":valid}
    for t in THRESHOLDS: payload[f"hits_{str(t).replace('.','p')}"]=hits[t]
    np.savez_compressed(CHECKPOINT,**payload)

def main():
    print("="*78)
    print("OSAF PROGRAM 36 — DAILY PWP OCCURRENCE / PERSISTENCE CLIMATOLOGY")
    print("="*78)
    dates=common_dates()
    if len(dates)==0: raise ValueError("No common dates.")
    print(f"Common record: {len(dates):,} days | {dates.min():%Y-%m-%d} to {dates.max():%Y-%m-%d}")

    first=P17.load_daily_field(dates[0])
    shape=first.sst_c.shape
    mask=first.pacific_mask.astype(bool)
    lat=first.latitude.astype(float); lon=first.longitude.astype(float)
    start,valid_count,hits=load_checkpoint(shape)
    print(f"Resume index: {start:,}")

    for i in range(start,len(dates)):
        field=first if i==0 else P17.load_daily_field(dates[i])
        if field.sst_c.shape != shape: raise ValueError(f"Grid shape changed on {dates[i]:%Y-%m-%d}")
        valid=mask & np.isfinite(field.sst_c)
        valid_count[valid]+=1
        for t in THRESHOLDS:
            hits[t][valid & (field.sst_c>=t)]+=1
        if (i+1)%CHECKPOINT_EVERY==0 or i+1==len(dates):
            save_checkpoint(i+1,valid_count,hits)
            print(f"{i+1:6,d}/{len(dates):,} | {dates[i]:%Y-%m-%d}")

    frac={}
    for t in THRESHOLDS:
        a=np.full(shape,np.nan,float)
        ok=mask & (valid_count>0)
        a[ok]=hits[t][ok]/valid_count[ok]
        frac[t]=a

    PROCESSED.mkdir(parents=True,exist_ok=True)
    TABLE_DIR.mkdir(parents=True,exist_ok=True)
    FIG_DIR.mkdir(parents=True,exist_ok=True)
    REPORT_DIR.mkdir(parents=True,exist_ok=True)

    tag=f"{dates.min():%Y-%m-%d}_{dates.max():%Y-%m-%d}"
    npz=PROCESSED/f"pwp_daily_occurrence_persistence_{tag}.npz"
    np.savez_compressed(
        npz, latitude=lat, longitude=lon, pacific_mask=mask.astype(np.uint8),
        valid_count=valid_count,
        occurrence_28p0=frac[28.0], occurrence_28p5=frac[28.5], occurrence_29p0=frac[29.0],
        thresholds_c=np.array(THRESHOLDS), start_date=str(dates.min().date()),
        end_date=str(dates.max().date()), n_common_days=np.array(len(dates),dtype=np.int64)
    )

    rows=[]
    for t in THRESHOLDS:
        x=frac[t][mask & np.isfinite(frac[t])]
        rows.append({
            "threshold_c":t,"n_common_days":len(dates),
            "start_date":str(dates.min().date()),"end_date":str(dates.max().date()),
            "pacific_cells_with_valid_denominator":int(x.size),
            "mean_gridcell_occurrence_fraction":float(np.mean(x)),
            "median_gridcell_occurrence_fraction":float(np.median(x)),
            "gridcells_occurrence_ge_0p10":int(np.sum(x>=.10)),
            "gridcells_occurrence_ge_0p50":int(np.sum(x>=.50)),
            "gridcells_occurrence_ge_0p90":int(np.sum(x>=.90)),
        })
    summary=pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR/"pwp_occurrence_persistence_summary.csv",index=False,float_format="%.8f")

    # Publication figure: identical axes and color normalization.
    fig,axes=plt.subplots(3,1,figsize=(11,10),sharex=True,sharey=True,constrained_layout=True)
    LON,LAT=np.meshgrid(lon,lat)
    mappable=None
    for ax,t,label in zip(axes,THRESHOLDS,["A)","B)","C)"]):
        mappable=ax.pcolormesh(LON,LAT,100*frac[t],shading="auto",vmin=0,vmax=100,cmap="viridis")
        ax.contour(LON,LAT,mask.astype(float),levels=[0.5],colors="0.35",linewidths=0.45)
        ax.set_ylabel("Latitude (°N)")
        ax.set_title(f"{label} Daily PWP occurrence — SST ≥ {t:.1f} °C")
        ax.grid(False)
    axes[-1].set_xlabel("Longitude (°E)")
    cbar=fig.colorbar(mappable,ax=axes,orientation="vertical",fraction=0.025,pad=0.02)
    cbar.set_label("Fraction of valid days satisfying threshold (%)")
    fig.suptitle(f"Pacific Warm Pool daily occurrence/persistence | {dates.min():%Y-%m-%d}–{dates.max():%Y-%m-%d} | N={len(dates):,}",fontsize=13)
    png=FIG_DIR/"pwp_daily_occurrence_persistence_thresholds.png"
    pdf=FIG_DIR/"pwp_daily_occurrence_persistence_thresholds.pdf"
    fig.savefig(png,dpi=300,bbox_inches="tight")
    fig.savefig(pdf,bbox_inches="tight")
    plt.close(fig)

    report=[
        "OSAF PROGRAM 36 — DAILY PWP OCCURRENCE / PERSISTENCE CLIMATOLOGY","",
        f"Program 17: {P17_FILE}",
        f"Common record: {len(dates):,} days | {dates.min():%Y-%m-%d} to {dates.max():%Y-%m-%d}",
        "Definition: per-cell fraction of VALID daily observations satisfying SST >= threshold.",
        "IMPORTANT: daily fields were thresholded individually; long-term mean SST was NOT thresholded.",
        f"Canonical processed output: {npz}",f"Figure: {png}",f"Figure: {pdf}","",
        "SUMMARY"
    ]
    for _,r in summary.iterrows():
        report.append(
            f"{r.threshold_c:.1f} °C | valid Pacific cells={int(r.pacific_cells_with_valid_denominator):,} | "
            f"mean cell occurrence={100*r.mean_gridcell_occurrence_fraction:.2f}% | "
            f"median={100*r.median_gridcell_occurrence_fraction:.2f}%"
        )
    (REPORT_DIR/"PROGRAM36_PWP_DAILY_OCCURRENCE_PERSISTENCE.txt").write_text("\n".join(report),encoding="utf-8")
    if CHECKPOINT.is_file(): CHECKPOINT.unlink()
    print("\n".join(report))
    print("\nPROGRAM 36 COMPLETED SUCCESSFULLY.")

if __name__=="__main__":
    main()
