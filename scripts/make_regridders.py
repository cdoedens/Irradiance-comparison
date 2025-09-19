import xarray as xr
from pathlib import Path
import xesmf as xe

from read_datasets import read_dataset

# date is arbitrary as grid is always the same
date  ='2020-01'
himawari = read_dataset(
            dataset='himawari',
            resolution='daily',
            date=date
        )

barra_r2 = read_dataset(
            dataset='barra-r2',
            resolution='daily',
            date=date
        )

him_ct = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
him_ct = him_ct.sel(lat=slice(-10, -44.5), lon=slice(112, 156.26))

file_path = Path("/g/data/er8/users/cd3022/regridder_weights")
file_path.mkdir(parents=True, exist_ok=True)
regridder = xe.Regridder(himawari, him_ct, "bilinear",
                         filename= f'{file_path}/himawari_to_himct_weights.nc',
                         reuse_weights=False)