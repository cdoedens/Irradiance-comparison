#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import datetime
import xesmf as xe
import re

from dask.distributed import Client, wait
import os, sys
import psutil

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset

LOG = logger.get_logger(__name__)

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

    file_path = Path(f"/scratch/er8/cd3022/Irradiance-comparisons/himawari_barrar2_diffct.zarr")
    ds = xr.open_zarr(file_path)

    # BARRA Cloud Types
    BARRA = Path('/g/data/ob53/BARRA2/output/reanalysis/')
    BARRA_R2 = BARRA / "AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
    BARRA_R2_DIR = BARRA_R2 / "1hr"

    vars_to_open = [
        'clh',
        'clm',
        'cll'
    ]
    def _preprocess(ds):
        return ds.sel(lat=slice(-44.5, -10), lon=slice(112, 156.26), time=slice('2016', '2024'))

    ds_vars = []
    for var in vars_to_open:

        files = sorted([f for f in BARRA_R2_DIR.glob(f'{var}/latest/*1hr_20*.nc')])
        ds_b = xr.open_mfdataset(files,
                            chunks={'time':24, 'lat':157, 'lon':31},
                            concat_dim='time',
                            combine='nested',
                            data_vars='minimal',
                            coords='minimal',
                            compat='override',
                            parallel=True,
                            preprocess=_preprocess
                            )
        ds_vars.append(ds_b)
    barra_r2 = xr.merge(ds_vars)
    LOG.info('BARRA-R2 opened')

    # assign cloud_type based on fractions at different levels
    barra_ct = (
        (barra_r2.cll >= 50).astype(int) * 1 +
        (barra_r2.clm >= 50).astype(int) * 2 +
        (barra_r2.clh >= 50).astype(int) * 4
    )
    barra_ct.attrs['description'] = '0: all < 50%; 1: low >= 50; 2: mid >= 50%; 3: low and mid >= 50%; 4: high >= 50%; 5: low and high >= 50%; 6: mid and high >= 50%; 7: all >= 50%'

    him_ct_18 = ds.ct.sel(time='2018-01')
    diff_18 = ds.ghi_diff.sel(time=him_ct_18.time, method='nearest')
    barra_ct_18 = barra_ct.sel(time=him_ct_18.time, method='nearest')

    cloud_aggs = {
        'clear_sky':[1,2,3,4],                  # 0
        'low':[5,6],                            # 1
        'med':[7],                              # 2
        'high_opaque':[8,9],                    # 3
        'fractional': [10],                     # 4
        'high_semitransparent': [11,12,13,14]   # 5
    }

    for i, agg in enumerate(cloud_aggs):
        ct_codes = cloud_aggs[agg]
        him_ct_18 = xr.where(him_ct_18.isin(ct_codes), i, him_ct_18)

    h_categories = {int(num): desc for num, desc in enumerate(cloud_aggs.keys())}

    b_matches = re.findall(r'(\d+):\s+([^;]+)', barra_ct.attrs['description'])
    b_categories = {int(num): desc.strip() for num, desc in b_matches}

    # him_ct_18, barra_ct_18, diff_18 = xr.align(him_ct_18, barra_ct_18, diff_18, join="inner")

    h_ct = him_ct_18.values.ravel()
    b_ct = barra_ct_18.values.ravel()
    diff_vals = diff_18.values.ravel()

    df = pd.DataFrame({
        "h_ct": h_ct,
        "b_ct": b_ct,
        "diff": diff_vals
    }).dropna(subset=["h_ct", "b_ct", "diff"])

    # Mean bias error (signed mean)
    mbe = df.groupby(["h_ct", "b_ct"])["diff"].mean().unstack(fill_value=np.nan)

    # Reindex to fill all possible categories
    mbe = mbe.reindex(
        index=range(len(h_categories)),
        columns=range(len(b_categories)),
        fill_value=np.nan
    )

    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(
        mbe,
        cmap="RdBu_r",          # Red = positive bias, Blue = negative bias
        center=0,               # zero-centered colormap
        annot=True,
        fmt=".1f",
        cbar_kws={"label": "Mean Bias Error (W/m²)"},
        square=True
    )

    # Center ticks on boxes
    ax.set_xticks(np.arange(mbe.shape[1]) + 0.5)
    ax.set_yticks(np.arange(mbe.shape[0]) + 0.5)

    ax.set_xticklabels(list(b_categories.values()), rotation=30, ha='right')
    ax.set_yticklabels(list(h_categories.values()))

    plt.xlabel("BARRA cloud type")
    plt.ylabel("Himawari cloud type")
    plt.title("Mean Bias Error of Surface SW Irradiance by Cloud Type (Jan 2018)")
    plt.tight_layout()
    plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/ct_matrix_mbe.png')
    plt.close()




    # cloud optical depth and cloud top pressure vs errors
