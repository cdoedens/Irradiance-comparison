#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
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
# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28

if __name__ == '__main__':
    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

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
        ds = xr.open_mfdataset(files,
                            chunks={'time':24, 'lat':256, 'lon':256},
                            concat_dim='time',
                            combine='nested',
                            data_vars='minimal',
                            coords='minimal',
                            compat='override',
                            parallel=True,
                            preprocess=_preprocess
                            )
        ds_vars.append(ds)
    barra_r2 = xr.merge(ds_vars)
    LOG.info('BARRA-R2 opened')

    # assign cloud_type based on fractions at different levels
    barra_ct = (
        (barra_r2.cll >= 50).astype(int) * 1 +
        (barra_r2.clm >= 50).astype(int) * 2 +
        (barra_r2.clh >= 50).astype(int) * 4
    )
    LOG.info('Cloud type data array created')
    barra_ct.attrs['data_description'] = '0: all cloud fractall < 50%. 1: cll >= 50%. 2: clm >= 50%. 3: cll and clm >= 50%. 4: clh >= 50%. 5: cll and clh >= 50%. 6: clm and clh >= 50%. 7: all >= 50%'

    # match chunks from barra with ds
    if "chunks" in barra_ct.encoding:
        del barra_ct.encoding["chunks"]
    barra_ct = barra_ct.chunk({"time": 24, "lat": 157, "lon": 31})
    LOG.info('BARRA cloud type data chunked')

    file_path = Path(f"/scratch/er8/cd3022/Irradiance-comparisons/himawari_barrar2_diffct.zarr")
    ds = xr.open_zarr(file_path)
    LOG.info('Opened zarr data store')

    ds['barra_ct'] = barra_ct
    LOG.info('Added barra_ct to data store')

    ds.to_zarr(file_path, mode="a", consolidated=False, zarr_format=2)
    LOG.info('Data store saved')