import xarray as xr
from pathlib import Path
import logger
from config import (
    HIMAWARI,
    ERA5,
    BARRA_R2,
    BARRA_C2,
    CERES
)

LOG = logger.get_logger(__name__)

######################################
# TO DO:
# 1. Generalise to more variables
# 2. Add more datasets
# 3. Flexible time range
######################################

def read_dataset(dataset, resolution, date):
    dataset = dataset.lower()
    resolution = resolution.lower()

    if resolution not in ('hourly', 'daily', 'monthly'):
        raise ValueError(f'Invalid resolution! Resolution must be: "hourly", "daily", or "monthly"')
    
    DATASET_READER = {
        'himawari': read_himawari,
        'era5': read_era5,
        'barra-r2': read_barra_r2,
        'barra-c2': read_barra_c2,
        'ceres': read_ceres
    }

    try:
        reader = DATASET_READER[dataset]
    except KeyError:
        raise ValueError(f'Unsupported dataset "{dataset}"! Supported datasets are: "himawari", "era5", "barra-r2", "barra-c2", "ceres"')
    return reader(resolution, date)

def read_himawari(resolution, date):
    if resolution == 'hourly':
        HIMAWARI_DIR = HIMAWARI / 'p1h'
    elif resolution == 'daily':
        HIMAWARI_DIR = HIMAWARI / 'p1d'

    from datetime import datetime
    date_dt = datetime.strptime(date, '%Y-%m')
    if date_dt <= datetime.strptime('2019-03-31', '%Y-%m-%d'):
        version = 'v1.0'
    else:
        version = 'v1.1'

    year, month = date.split('-')

    directory = HIMAWARI_DIR / version / year/ month
    files = [f for f in directory.rglob('*.nc')]

    renaming_dir = {
        'latitude': 'lat',
        'longitude': 'lon',
        f'{resolution}_integral_of_surface_global_irradiance': 'ghi'
        }
    ds = xr.open_mfdataset(files).rename(
        renaming_dir
    )
    # convert units from MJ to W
    ds['ghi'] = ds['ghi'] * (1e6 / 3600)
    if resolution == 'daily':
        ds['ghi'] = ds['ghi'] / 24
    ds['ghi'].attrs['units'] = 'W m-2'
    return ds

def read_era5(resolution, date):
    if resolution == 'hourly':
        ERA5_DIR = ERA5 / "reanalysis"
    elif resolution == 'monthly':
        ERA5_DIR = ERA5 / "monthly-averaged"

    year, month = date.split('-')

    file = ERA5_DIR.glob(f"msdwswrf/{year}/*{year}{month}*.nc")

    renaming_dir = {
        'longitude': 'lon',
        'latitude': 'lat',
        'msdwswrf': 'ghi'
    }

    return xr.open_mfdataset(file).rename(
        renaming_dir
        ).sel(
            lat=slice(-10, -44.5), lon=slice(112, 156.26)
            )

def read_barra_r2(resolution, date):
    if resolution == 'hourly':
        BARRA_R2_DIR = BARRA_R2 / "1hr"
    elif resolution == 'daily':
        BARRA_R2_DIR = BARRA_R2 / 'day'
    elif resolution == 'monthly':
        BARRA_R2_DIR = BARRA_R2 / 'mon'

    year, month = date.split('-')

    file = BARRA_R2_DIR.glob(f'rsds/latest/*{year}{month}*.nc')

    renaming_dir = {
        'rsds': 'ghi'
    }

    return xr.open_mfdataset(file).rename(
        renaming_dir
    ).sel(
        lat=slice(-44.5, -10,), lon = slice(112, 156.26)
    )

def read_barra_c2(resolution, date):
    if resolution == 'hourly':
        BARRA_C2_DIR = BARRA_C2 / "1hr"
    elif resolution == 'daily':
        BARRA_C2_DIR = BARRA_C2 / 'day'
    elif resolution == 'monthly':
        BARRA_C2_DIR = BARRA_C2 / 'mon'

    year, month = date.split('-')

    file = BARRA_C2_DIR.glob(f'rsds/latest/*{year}{month}*.nc')

    renaming_dir = {
        'rsds': 'ghi'
    }

    return xr.open_mfdataset(file).rename(
        renaming_dir
    ).sel(
        lat=slice(-44.5, -10,), lon = slice(112, 156.26)
    )

def read_ceres(resolution, date):
    LOG.info('CERES only has monthly resolution')
    year, month = date.split('-')
    renaming_dir = {
        'sfc_sw_down_all_mon': 'ghi'
    }
    return xr.open_dataset(CERES).rename(
        renaming_dir
    ).sel(
        lat=slice(-44.5, -10), lon=slice(112, 156.26), time=date
    )
    
