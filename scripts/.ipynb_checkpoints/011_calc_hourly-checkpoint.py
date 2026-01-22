import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import datetime

import os, sys
sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
from read_datasets import read_dataset
import logger
LOG = logger.get_logger(__name__)

from dask.distributed import Client


if __name__ == '__main__':

    client = Client(
        n_workers=24,
        threads_per_worker=1
    )

    dir_path = "/g/data/er8/users/cd3022/Irradiance-comparisons/"
    dataset = sys.argv[1]
    year = sys.argv[2]
    # month = int(sys.argv[3])
    resolution = 'instant'

    month_list = []
    for month in range(1, 13):
        date = f'{year}-{month:02d}'
        ds_month = read_dataset(
                dataset=dataset,
                resolution=resolution,
                date=date
            )

        
        LOG.info(f'Reading date {month}-{year} for dataset {dataset}')

        ghi_month = ds_month.ghi
        if dataset == "himawari":
            # fill missing time steps
            full_time = pd.date_range(
                start=ghi_month.time.min().item(),
                end=ghi_month.time.max().item(),
                freq="10min"
            )
            ghi_month = ghi_month.reindex(time=full_time, fill_value=0)
            ghi_month = ghi_month.resample(time="1h", label="left", closed="right").mean()
            LOG.info(f'Resampled to 1h')
        else:
            LOG.info('Only set up for Himawari atm')

        ghi_month = ghi_month.assign_coords({'month': month})
        month_list.append(ghi_month)

    year_ds = xr.concat(month_list, dim='time')
    year_hours = year_ds.groupby("time.hour").mean("time")
    LOG.info(f'Calculated hour means')

    year_hours.to_netcdf(f"{dir_path}hourly_{dataset}_{year}.nc")
    LOG.info(f'DONE')