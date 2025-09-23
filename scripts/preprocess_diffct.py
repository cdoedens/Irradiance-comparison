#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import xesmf as xe

import os, sys
import psutil

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset
from dask_setup import setup_dask_client

LOG = logger.get_logger(__name__)
# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28

client = setup_dask_client(workload_type="io")
LOG.info('Dask client started')


date='2020-01'


# LOAD HIMAWARI
himawari = read_dataset(
        dataset='himawari',
        resolution='hourly',
        date=date
    )
LOG.info('Himawari data opened')

# LOAD BARRA-R2

BARRA = Path('/g/data/ob53/BARRA2/output/reanalysis/')
BARRA_R2 = BARRA / "AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
BARRA_R2_DIR = BARRA_R2 / "1hr"

file = BARRA_R2_DIR.glob(f'rsds/latest/*.nc')

renaming_dir = {
    'rsds': 'ghi'
}

# time slice based on himawari
start = ds1.isel(time=0).time
end = ds1.isel(time=-1).time

barra_r2 = xr.open_mfdataset(file, chunks='auto', parallel=False).rename(
    renaming_dir
).sel(
    lat=slice(-44.5, -10),
    lon=slice(112, 156.26),
    time=slice(start, end)
)
LOG.info('BARRA-R2 data opened')


# REGRID HIMAWARI TO BARRA-R2
regridder =  xe.Regridder(ds1, barra_r2, "bilinear", reuse_weights=False)
ds1 = regridder(ds1)
LOG.info('Himawari regridded to BARRA-R2')
# GET DIFF
diff = ds1.ghi - barra_r2.ghi
LOG.info('Diff calculated')


# LOAD CLOUD TYPES
ds_ct = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
LOG.info('Cloud type data opened')

lat=slice(-10, -44.5)
lon=slice(112, 156.26)
ds_ct = ds_ct.sel(
    lat=lat,
    lon=lon,
    time=diff.time
)




# REGRID DIFFS TO CT
regridder_ct = xe.Regridder(diff, ds_ct, method='bilinear')
diff = regridder_ct(diff)
LOG.info('Diffs regridded to CT')


# COMBINE INTO SINGLE DATASET
final_ds = xr.Dataset(
    {
        "ghi_diff": diff,
        "ct": ds_ct.ct
    }
)

# Remove any pre-existing chunk encodings
for v in final_ds.variables:
    if "chunks" in final_ds[v].encoding:
        del final_ds[v].encoding["chunks"]

final_ds = final_ds.chunk({"time": 24, "lat": 128, "lon": 128})
LOG.info('Final_ds rechunked')

# SAVE AS ZARR
file_path = Path("/scratch/er8/cd3022/Irradiance-comparisons/")
os.makedirs(file_path, exist_ok=True)
file_name = f'himawari_barrar2_diffct_{date}.zarr'
final_ds["ghi_diff"].to_zarr(f"{file_path}/{file_name}", mode="w")
final_ds["ct"].to_zarr(f"{file_path}/{file_name}", mode="a")