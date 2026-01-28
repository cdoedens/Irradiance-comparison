import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74

# Load datasets
file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/cape')
files = [f for f in file_path.glob('hourly*.nc')]
ds = xr.open_mfdataset(files, engine='h5netcdf', combine='nested', concat_dim='year')


ds = ds.sel(
    lat=slice(-44.5, -10),
    lon=slice(112, 156.126)
)
ds = ds.mean(dim='year')

ds["dCAPE"] = ds["CAPE"] - ds["CAPE"].roll(hour=1, roll_coords=False)



# PLOT

fig, ax = plt.subplots(ncols=5, nrows=3, figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})
ax = ax.flatten()

for i, hour in enumerate([19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]):
    data = ds.sel(hour=hour)


    mesh = ax[i].pcolormesh(
        data.lon,
        data.lat,
        data.dCAPE,
        cmap='RdBu_r', shading='auto',
        vmin=-200, vmax=200,
        transform=ccrs.PlateCarree()
    )
    ax[i].coastlines()
    ax[i].set_title(f'{hour:02d}:00 UTC')
fig.suptitle(f'BARRA-R2 dCAPE/dt', fontsize=18)
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("CAPE (J kg-1 hr-1)")
plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/barra-r2_cape_hourly.png')
plt.close()