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

def setup_dask_client(
    workload_type="io",        # "cpu", "io", or "mixed"
    max_workers=None,           # Limit max workers (default = all cores)
    reserve_mem_gb=50,          # Memory reserved for system and overhead
    max_mem_gb=None,            # Cap usable memory if needed
    dashboard=True
):
    """
    Setup Dask client with auto-scaling for workload type.

    Parameters:
        workload_type (str): Type of workload - "cpu", "io", or "mixed"
        max_workers (int): Max logical cores to use (default = system max)
        reserve_mem_gb (int): System memory to reserve
        max_mem_gb (int): Cap memory usage (defaults to system memory)
        dashboard (bool): Whether to print the dashboard link

    Returns:
        dask.distributed.Client
    """
    assert workload_type in ("cpu", "io", "mixed"), "Invalid workload_type"

    logical_cores = psutil.cpu_count(logical=True)
    total_memory_gb = psutil.virtual_memory().total / 1e9

    if max_workers is None:
        max_workers = logical_cores

    if max_mem_gb is None:
        max_mem_gb = total_memory_gb

    usable_mem_gb = max_mem_gb - reserve_mem_gb

    # Recommended presets based on workload type
    if workload_type == "cpu":
        threads_per_worker = 1
        n_workers = min(max_workers, logical_cores)
    elif workload_type == "io":
        threads_per_worker = 8
        n_workers = max(1, logical_cores // threads_per_worker)
    else:  # "mixed"
        threads_per_worker = 4
        n_workers = max(1, logical_cores // threads_per_worker)

    memory_per_worker = usable_mem_gb // n_workers

    client = Client(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_limit=f"{int(memory_per_worker)}GB"
    )

    if dashboard:
        print(f"Dask dashboard: {client.dashboard_link}")

    return client

if __name__ == '__main__':
    client = setup_dask_client(workload_type="io")
    client


    date = '2020-05'
    himawari = read_dataset(
            dataset='himawari',
            resolution='hourly',
            date=date
        )
    
    BARRA = Path('/g/data/ob53/BARRA2/output/reanalysis/')
    BARRA_R2 = BARRA / "AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
    BARRA_R2_DIR = BARRA_R2 / "1hr"
    
    file = BARRA_R2_DIR.glob(f'rsds/latest/*.nc')
    
    renaming_dir = {
        'rsds': 'ghi'
    }
    
    # time slice based on himawari
    start = himawari.isel(time=0).time
    end = himawari.isel(time=-1).time
    
    barra_r2 = xr.open_mfdataset(file, chunks='auto', parallel=False).rename(
        renaming_dir
    ).sel(
        lat=slice(-44.5, -10),
        lon=slice(112, 156.26),
        time=slice(start, end)
    )
    
    regridder =  xe.Regridder(himawari, barra_r2, "bilinear", reuse_weights=False)
    himawari = regridder(himawari)
    diff = himawari.ghi - barra_r2.ghi
    ds_ct = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
    
    lat=slice(-10, -44.5)
    lon=slice(112, 156.26)
    ds_ct = ds_ct.sel(
        lat=lat,
        lon=lon,
    )
    ds_ct = ds_ct.sel(
        time=diff.time,
        method='nearest'
    )
    
    regridder_ct = xe.Regridder(ds_ct, diff, method='bilinear')
    ds_ct = regridder_ct(ds_ct)
    
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
    
    file_path = Path("/scratch/er8/cd3022/Irradiance-comparisons/")
    os.makedirs(file_path, exist_ok=True)
    file_name = f'himawari_barrar2_diffct_{date}'
    # final_ds.to_netcdf(f'{file_path}/{file_name}.nc')
    final_ds.to_zarr(f"{file_path}/{file_name}.zarr", mode="w")