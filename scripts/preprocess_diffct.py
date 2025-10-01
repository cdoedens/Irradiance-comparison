#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import xesmf as xe

from dask.distributed import Client, wait
import os, sys
import psutil

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset

LOG = logger.get_logger(__name__)
# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28

if __name__ == '__main__':
    client = Client(
        n_workers=48,
        threads_per_worker=1
    )
    LOG.info('Dask client started')

    # OPEN BARRA-R2 DATA
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
    LOG.info('BARRA-R2 data opened')

    # OPEN HIMAWARI CT DATA
    ds_ct = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
    
    lat=slice(-10, -44.5)
    lon=slice(112, 156.26)
    ds_ct = ds_ct.sel(
        lat=lat,
        lon=lon,
    )
    LOG.info('Himawari CT data opened')
    
    year = 2017
    LOG.info(f'Beginning loop for year: {year}')
    for month in range (1, 13):
        LOG.info(f'starting month: {month}')
        date = f'{year}-{month:02d}' 
        himawari = read_dataset(
                dataset='himawari',
                resolution='hourly',
                date=date
            )
        LOG.info('Himawari irradiance opened')
        
        # time slice based on himawari
        start = himawari.isel(time=0).time
        end = himawari.isel(time=-1).time
        barra_date = barra_r2.sel(
            time=slice(start, end)
        )
    
        
        regridder =  xe.Regridder(himawari, barra_date, "bilinear", reuse_weights=False)
        himawari = regridder(himawari)
        LOG.info('Himawari irradiance regridded to BARRA-R2')
        diff = himawari.ghi - barra_date.rsds
        LOG.info('Calculated diff between Himawari and BARRA-R2 irradiance')
    
        ds_ct_date = ds_ct.sel(
            time=diff.time,
            method='nearest'
        )
        
        regridder_ct = xe.Regridder(ds_ct_date, diff, method='bilinear')
        ds_ct_date = regridder_ct(ds_ct_date)
        LOG.info('Himawari CT regridded to diff')
    
        if diff["time"].to_index().duplicated().any():
            diff = diff.sel(time=~diff.get_index("time").duplicated())
            LOG.info('diff duplicates removed')
        if ds_ct_date["time"].to_index().duplicated().any():
            ds_ct_date = ds_ct_date.sel(time=~ds_ct_date.get_index("time").duplicated())
            LOG.info('CT duplicates removed')
        
        final_ds = xr.Dataset(
            {
                "ghi_diff": diff,
                "ct": ds_ct_date.ct
            }
        )
        LOG.info('final_ds ')
        
        # Remove any pre-existing chunk encodings
        for v in final_ds:
            if "chunks" in final_ds[v].encoding:
                del final_ds[v].encoding["chunks"]
        
        final_ds = final_ds.chunk({"time": -1, "lat": 157, "lon": 31})
        # ensure all variables match these chunks
    
        
        file_path = Path("/scratch/er8/cd3022/Irradiance-comparisons/")
        os.makedirs(file_path, exist_ok=True)
        file_name = f'himawari_barrar2_diffct'
        
        if not os.path.exists(f"{file_path}/{file_name}.zarr"):
            final_ds.to_zarr(f"{file_path}/{file_name}.zarr", mode="w", consolidated=True)
            print(f'writing {date}')
        else:
            final_ds.to_zarr(f"{file_path}/{file_name}.zarr", append_dim="time", consolidated=True)
            print(f'appending {date}')