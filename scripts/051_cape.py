import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
from dask.distributed import Client
import sys, os

if __name__ == '__main__':
    client = Client(
                n_workers=24,
                threads_per_worker=1
            )

    year = sys.argv[1]

    file_path = Path('/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/CAPE/latest/')
    files = [f for f in file_path.glob(f'*{year}??.nc')]

    ds = xr.open_mfdataset(files)
    ds = ds.assign_coords({'year': year})

    ds_hours = ds.groupby('time.hour').mean(dim='time')
    ds_seasons = ds.groupby('time.season').mean(dim='time')

    ds_hours.attrs['history'] = datetime.date.today().strftime('%D')
    ds_hours.attrs['source_script'] = 'data produced by the script "051_cape.py"'
    ds_seasons.attrs['history'] = datetime.date.today().strftime('%D')
    ds_seasons.attrs['source_script'] = 'data produced by the script "051_cape.py"'

    save_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/cape')
    os.makedirs(save_path, exit_ok=True)
    ds_hours.to_netcdf(f'{save_path}/hourly_barra-r2_{year}.nc')
    ds_seasons.to_netcdf(f'{save_path}/seasons_barra-r2_{year}.nc')