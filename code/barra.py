import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

barra_path = Path('/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/tas/latest/tas_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_202501-202501.nc')

ds = xr.open_dataset(barra_path)