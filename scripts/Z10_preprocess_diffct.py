#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

from dask.distributed import Client, wait
import os, sys
import psutil

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset

LOG = logger.get_logger(__name__)
year = sys.argv[1]
LOG.info(f'Running script for year: {year}')

if __name__ == '__main__':
    client = Client(
        n_workers=48,
        threads_per_worker=1
    )
    LOG.info('Started client')

    # OPEN BARRA-R2
    BARRA = Path('/g/data/ob53/BARRA2/output/reanalysis/')
    BARRA_R2 = BARRA / "AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
    BARRA_R2_DIR = BARRA_R2 / "1hr"
    var = 'rsds'
    files = sorted([f for f in BARRA_R2_DIR.glob(f'{var}/latest/*1hr_20*.nc')])
    def _preprocess(ds):
        return ds.sel(lat=slice(-44.5, -10), lon=slice(112, 156.26))
    barra_r2 = xr.open_mfdataset(
        files,
        chunks={'time':24, 'lat':256, 'lon':256},
        concat_dim='time',
        combine='nested',
        data_vars='minimal',
        coords='minimal',
        compat='override',
        parallel=True,
        preprocess=_preprocess
    )
    LOG.info('BARRA-R2 opened')

    # OPEN HIMAWARI CT
    ds_ct = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
    lat=slice(-10, -44.5)
    lon=slice(112, 156.26)
    ds_ct = ds_ct.sel(
        lat=lat,
        lon=lon,
    )
    LOG.info('Himawari CT opened')
    
    # OPEN HIMAWARI GHI
    himawari_list = []
    for month in range (1, 13):
        date = f'{year}-{month:02d}' 
        himawari = read_dataset(
                dataset='himawari',
                resolution='hourly',
                date=date
            )
        himawari_list.append(himawari)
    him_ds = xr.concat(himawari_list, dim='time', data_vars='minimal')
    LOG.info('Himawari GHI opened')

    # PROCESS DATA
    start = him_ds.isel(time=0).time
    end = him_ds.isel(time=-1).time
    barra_date = barra_r2.sel(
        time=slice(start, end)
    )
    LOG.info('BARR-R2 time slice')
    
    him_ds = him_ds.interp(
        lat=barra_date.lat,
        lon=barra_date.lon,
        method='linear'
    )
    LOG.info('Himawari GHI regridded to BARRA-R2')
    diff = him_ds.ghi - barra_date.rsds

    ds_ct_date = ds_ct.sel(
        time=diff.time,
        method='nearest'
    )
    ds_ct_date = ds_ct_date.interp(
        lat=diff.lat,
        lon=diff.lon,
        method='nearest'
    )
    LOG.info('Himawari CT regridded to diff')
    
    if diff["time"].to_index().duplicated().any():
        diff = diff.sel(time=~diff.get_index("time").duplicated())
    
    if ds_ct_date["time"].to_index().duplicated().any():
        ds_ct_date = ds_ct_date.sel(time=~ds_ct_date.get_index("time").duplicated())

    # CONSTRUCT DATASET TO BE SAVED
    final_ds = xr.Dataset(
        {
            "ghi_diff": diff,
            "ct": ds_ct_date.ct
        }
    )
    LOG.info('Final DS constructed')
    
    # Remove any pre-existing chunk encodings
    for v in final_ds:
        if "chunks" in final_ds[v].encoding:
            del final_ds[v].encoding["chunks"]
    
    final_ds = final_ds.chunk({"time": 24, "lat": 157, "lon": 31})
    LOG.info('Final DS chunked')
    # ensure all variables match these chunks
    
    
    file_path = Path("/scratch/er8/cd3022/Irradiance-comparisons/")
    os.makedirs(file_path, exist_ok=True)
    file_name = f'himawari_barrar2_diffct_{year}'
    LOG.info('Saving DS as zarr')
    final_ds.to_zarr(f"{file_path}/{file_name}.zarr", mode="w", consolidated=True, zarr_format=2)
    LOG.info('DONE')