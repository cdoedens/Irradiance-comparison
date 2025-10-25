import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74

dataset_1 = sys.argv[1]
dataset_2 = sys.argv[2]

# Load datasets
file_path_1 = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/{dataset_1}_ghi_means')
ds1 = xr.open_dataset(f'{file_path_1}/hourly_{dataset_1}.nc')

file_path_2 = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/{dataset_2}_ghi_means')
ds2 = xr.open_dataset(f'{file_path_2}/hourly_{dataset_2}.nc')

# Prepare data

# regrid fine resolution to coarse resolution so datasets match
if ds1.ghi.shape[0] > ds2.ghi.shape[0]:
    ds1 = ds1.interp(
        lat=ds2.lat,
        lon=ds2.lon,
        method='nearest' # 'linear', 'nearest', 'cubic'
    )
else:
    ds2 = ds2.interp(
        lat=ds1.lat,
        lon=ds1.lon,
        method='linear'
    )
diff = ds2 - ds1


# PLOT

fig, ax = plt.subplots(ncols=6, nrows=4, figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})
ax = ax.flatten()

for i, hour in enumerate(ds1.hour):
    diff_data = diff.sel(hour=hour)


    mesh = ax[i].pcolormesh(
        diff_data.lon,
        diff_data.lat,
        diff_data.ghi,
        cmap='RdBu_r', shading='auto',
        vmin=-200, vmax=200,
        transform=ccrs.PlateCarree()
    )
    ax[i].coastlines()
    ax[i].set_title(f'{hour:02d}:00 UTC')
fig.suptitle(f'{dataset_2.upper()} - {dataset_1.upper()}', fontsize=18)
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("Surface Downward SW bias (W/m²)")
plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/{dataset_2}-{dataset_1}_hourly.png')