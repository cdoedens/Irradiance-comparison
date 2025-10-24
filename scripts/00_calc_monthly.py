#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logger
import datetime

from read_datasets import read_dataset

from dask.distributed import Client, as_completed
import os, sys

dataset = sys.argv[1]
resolution = sys.argv[2]

LOG = logger.get_logger(__name__)
def single_month_mean(date):
    '''
    LOAD DATA FOR A MONTH,
    RETURN THE MEAN FOR THAT MONTH
    '''
    ds = read_dataset(
        dataset=dataset,
        resolution=resolution,
        date=date
    )
    LOG.info(f'data opened: {dataset}, {resolution}, {date}')
    da = ds['ghi']
    return da.mean(dim='time')

def multi_month_mean(ds):
    return ds.mean(dim='date', skipna=True)




if __name__ == '__main__':
    client = Client(
        n_workers = 24,
        memory_limit = int(os.environ['PBS_VMEM']) / int(os.environ['PBS_NCPUS']),
    )
    LOG.info('Dask client started')
    LOG.info(f'Dataset = {dataset}')
    LOG.info(f'Resolution = {resolution}')

    '''
    Inner Loop: calculate mean for an individual month, looping over years
    Outer Loop: calculate mean for a month across all the years, looping over months
    '''

    futures = {}
    for month in range(1,13):

        # submit inner jobs to client
        futures_2 = {}
        for year in range(2016, 2025):
            date = f'{year}-{month:02d}'
            future_2 = client.submit(single_month_mean, date)
            LOG.info(f'inner loop future_2 submitted to client with date = {date}')
            futures_2[future_2] = date


        # collect results for the month
        results_2 = []
        labels_2 = []
        for future_2 in as_completed(futures_2):
           result_2 = future_2.result()
           LOG.info('single month mean calculated')
           results_2.append(result_2)
           labels_2.append(futures_2[future_2])
        combined = xr.concat(results_2, dim='date')
        combined = combined.assign_coords(date=labels_2)

        # Submit outer jobs to client
        future = client.submit(multi_month_mean, combined)
        futures[future] = f'{month:02d}'
        LOG.info('outer loop future submitted to client')
    
    # collect results for ALL months
    results = []
    labels = []
    for future in as_completed(futures):
        result = future.result()
        LOG.info('multi month mean calculated')
        results.append(result)
        labels.append(futures[future])

    all_months = xr.concat(results, dim='month')
    all_months = all_months.assign_coords(month=labels)

    # assign attributes
    all_months.attrs["units"] = "Wm-2"
    script_name = Path(__file__).name
    all_months.attrs["source"] = script_name
    current_time = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    all_months.attrs["history"] = f"Data created on: {current_time}"

    file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/{dataset}_ghi_means/')
    file_name = f'{dataset}_monthly.nc'
    os.makedirs(file_path, exist_ok=True)
    all_months.to_netcdf(f'{file_path}/{file_name}')
