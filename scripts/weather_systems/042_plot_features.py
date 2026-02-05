import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74

wf = sys.argv[1]

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

# Weather feature frequency
def _preprocess(_ds):
    return _ds.sel(
        latitude=slice(him.lat.max(), him.lat.min()),
        longitude=slice(him.lon.min(), him.lon.max()),
    ).rename(
        {'latitude':'lat',
        'longitude':'lon'}
    )

file_path = Path('/scratch/er8/cd3022/weather-features/cyc-anticyc-front-cao')
wf_files = [f for f in file_path.glob(f'ea.ans.20*{wf}.nc')]
wf_ds = xr.open_mfdataset(wf_files, preprocess=_preprocess)
wf_ds = wf_ds.sel(time=slice('2016', '2022'))

# remove the level dimension if looking at frontal volume
if wf == 'frovo_id':
    wf_ds = wf_ds.isel(lev=0)

wf_ds = wf_ds.interp(
    lat=bar.lat,
    lon=bar.lon,
    method='nearest'
)

# make sure data is binary
wf_da = xr.where(wf_ds[wf] > 0, 1, wf_ds[wf])

wf_count = wf_da.sum(dim='time')
wf_tot = wf_da.sizes['time']
wf_freq = wf_count / wf_tot

# combine difference and freq to represent total effect of the weather feature
comb = diff * wf_freq

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20,5), subplot_kw={'projection': ccrs.PlateCarree()})

mesh_freq = ax[0].pcolormesh(
    wf_freq.lon,
    wf_freq.lat,
    wf_freq,
    cmap='magma', shading='auto',
    vmin=0, vmax=0.1,
    transform=ccrs.PlateCarree()
)
mesh_diff = ax[1].pcolormesh(
    diff.lon,
    diff.lat,
    diff,
    cmap='RdBu_r', shading='auto',
    vmin=-200, vmax=200,
    transform=ccrs.PlateCarree()
)
mesh_comb = ax[2].pcolormesh(
    comb.lon,
    comb.lat,
    comb,
    cmap='RdBu_r', shading='auto',
    vmin=-20, vmax=20,
    transform=ccrs.PlateCarree()
)
ax[0].coastlines()
ax[1].coastlines()
ax[2].coastlines()
fig.suptitle(f'BARRA-R2 - Himwari {wf}', fontsize=18)

cbar_freq= fig.colorbar(mesh_freq, ax=ax[0], orientation='horizontal', fraction=0.05, pad=0.04)
cbar_freq.set_label(f"Frequency of {wf} [-]")
cbar_diff = fig.colorbar(mesh_diff, ax=ax[1], orientation='horizontal', fraction=0.05, pad=0.04)
cbar_diff.set_label("Surface Downward SW bias (W/m²)")
cbar_comb = fig.colorbar(mesh_comb, ax=ax[2], orientation='horizontal', fraction=0.05, pad=0.04)
cbar_comb.set_label("Frequency x Bias (W/m²)")

plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/weather-features/barra-r2_himawari_{wf}.png')
