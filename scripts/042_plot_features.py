import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74

wf = 'cycmask'

# Load datasets
file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features/{wf}')
him_files = [f for f in file_path.glob('himawari*.nc')]
barra_files = [f for f in file_path.glob('barra-r2*.nc')]
him = xr.open_mfdataset(him_files, combine='nested', concat_dim='year')
bar = xr.open_mfdataset(barra_files, combine='nested', concat_dim='year')

him = him.annual_mean_ghi.mean(dim='year')
bar = bar.annual_mean_ghi.mean(dim='year')

him = him.interp(
    lat=bar.lat,
    lon=bar.lon,
    method='linear'
)

diff = bar - him


fig, ax = plt.subplots(figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})

mesh = ax.pcolormesh(
    diff.lon,
    diff.lat,
    diff,
    cmap='RdBu_r', shading='auto',
    vmin=-250, vmax=250,
    transform=ccrs.PlateCarree()
)
ax.coastlines()
fig.suptitle(f'BARRA-R2 - Himwari {wf}', fontsize=18)
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("Surface Downward SW bias (W/m²)")
plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/weather-features/barra-r2_himawari_{wf}.png')
