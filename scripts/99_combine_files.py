import xarray as xr
from pathlib import Path

ds_list = []
for year in range(2016, 2025):
    file_path = Path(f"/g/data/er8/users/cd3022/Irradiance-comparisons/hourly_himawari_{year}.nc")
    if not file_path.exists():
        print(f'Skipping year {year}, file not found')
        continue
    ds_year = xr.open_dataset(file_path, engine='h5netcdf')
    ds_year = ds_year.assign_coords({'year':year})
    ds_list.append(ds_year)
ds = xr.concat(ds_list, dim='year')
ds = ds.groupby('hour').mean('year')

ds.to_netcdf(f"/g/data/er8/users/cd3022/Irradiance-comparisons/himawari_ghi_means/hourly_himawari.nc")