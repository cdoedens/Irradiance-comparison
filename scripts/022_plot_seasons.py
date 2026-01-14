import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs
import sys

# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74



##############################################################################

# NOT COMPLETE YET, STILL FIXING

#####################################################################################################################








datasets = [
    'barra-r2',
    'barra-c2',
    'era5'
]
file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/')

himawari = xr.open_dataset(
    f'{file_path}himawari_ghi_means/himawari_monthly.nc',
    engine='h5netcdf'
)

re_ds_list = []
for re in datasets:
    ds = xr.open_dataset(
        f'{file_path}{re}_ghi_means/{re}_monthly.nc',
        engine='h5netcdf'
        )
    ds = ds.assign_coords({'dataset': re})
    re_ds_list.append(ds)
re_ds = xr.concat(re_ds_list)


# PLOT
season_mapping = {
    '01': 'DJF',
    '02': 'DJF',
    '03': 'MAM',
    '04': 'MAM',
    '05': 'MAM',
    '06': 'JJA',
    '07': 'JJA',
    '08': 'JJA',
    '09': 'SON',
    '10': 'SON',
    '11': 'SON',
    '12': 'DJF'
}

fig, ax = plt.subplots(ncols=4, nrows=3, figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})
ax = ax.flatten()

for i, month in enumerate(ds1.month):
    ds1_data = ds1.sel(month=month)
    ds2_data = ds2.sel(month=month)

    # regrid fine resolution to coarse resolution so datasets match
    if ds1_data.ghi.shape[0] > ds2_data.ghi.shape[0]:
        ds1_data = ds1_data.interp(
            lat=ds2_data.lat,
            lon=ds2_data.lon,
            method='nearest' # 'linear', 'nearest', 'cubic'
        )
    else:
        ds2_data = ds2_data.interp(
            lat=ds1_data.lat,
            lon=ds1_data.lon,
            method='linear'
        )
    diff = ds2_data - ds1_data
    mesh = ax[i].pcolormesh(
        diff.lon,
        diff.lat,
        diff.ghi,
        cmap='RdBu_r', shading='auto',
        vmin=-40, vmax=40,
        transform=ccrs.PlateCarree()
    )
    ax[i].coastlines()
    ax[i].set_title(month_names[str(month.data)])
fig.suptitle(f'{dataset_2.upper()} - {dataset_1.upper()}', fontsize=18)
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("Surface Downward SW bias (W/m²)")
plt.savefig(f'/home/548/cd3022/figures/Irradiance_comparisons/{dataset_2}-{dataset_1}_monthly.png')