import xarray as xr
from pathlib import Path

ds_list = []
for year in range(2016, 2025):
    file_path = Path(f"/scratch/er8/cd3022/Irradiance-comparisons/himawari_barrar2_diffct_{year}.zarr")
    ds_year = xr.open_zarr(file_path, consolidated=False)
    ds_list.append(ds_year)
ds = xr.concat(ds_list, dim='time')

# Remove any pre-existing chunk encodings
for v in ds:
    if "chunks" in ds[v].encoding:
        del ds[v].encoding["chunks"]

ds = ds.chunk({"time": 24, "lat": 157, "lon": 31})

ds.to_zarr(f"/scratch/er8/cd3022/Irradiance-comparisons/himawari_barrar2_diffct.zarr", mode="w", consolidated=True, zarr_format=2)