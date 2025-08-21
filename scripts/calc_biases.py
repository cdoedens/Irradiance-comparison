import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74

dataset = sys.argv[1]

# Load Himawari data
himawari = xr.open_dataset('/g/data/er8/users/cd3022/Irradiance-comparisons/himawari_ghi_means/himawari_monthly.nc')

# Load comparison dataset
file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/{dataset}_ghi_means')
ds = xr.open_dataset(f'{file_path}/{dataset}_monthly.nc')

# PLOT
month_names = {
    '01': 'jan',
    '02': 'feb',
    '03': 'mar',
    '04': 'apr',
    '05': 'may',
    '06': 'jun',
    '07': 'jul',
    '08': 'aug',
    '09': 'sep',
    '10': 'oct',
    '11': 'nov',
    '12': 'dec'
}

fig, ax = plt.subplots(ncols=4, nrows=3, figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})
ax = ax.flatten()

for i, month in enumerate(himawari.month):
    him_data = himawari.sel(month=month)
    ds_data = ds.sel(month=month)

    # regrid himawari to BARRA-C2 resolution
    him_on_ds_grid = him_data.interp(
        lat=ds_data.lat,
        lon=ds_data.lon,
        method="linear"   # or "nearest", "cubic" etc.
    )
    diff = ds_data - him_on_ds_grid
    # diff = diff * 277.78 / 24 # MJ to Wh (and hence avg W)
    lon2d, lat2d = np.meshgrid(diff['lon'].values, diff['lat'].values) 
    mesh = ax[i].pcolormesh(
        lon2d,
        lat2d,
        diff.ghi.values,
        cmap='RdBu_r', shading='auto',
        vmin=-40, vmax=40,
        transform=ccrs.PlateCarree()
    )
    ax[i].coastlines()
    ax[i].set_title(month_names[str(month.data)])
fig.suptitle(f'{dataset.upper()} - Himawari', fontsize=18)
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("Surface Downward SW bias (MJ/m²)")
plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/{dataset}-Himwari_monthly.png')