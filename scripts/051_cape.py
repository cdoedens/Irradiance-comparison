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

    file_path = Path('BARRA/PATH/TO/CAPE/HERE/')
    files = [f for f in file_path.glob(f'*{year}*.nc')]

    ds = xr.open_mfdataset(files)

    ds_hours = ds.groupby('time.hour').mean(dim='time')
    ds_seasons = ds.groupby('time.season').mean(dim='time')

    ds_hours.to_netcdf('.../DATA/PATH/file.nc')
    ds_seasons.to_netcdf('...')
