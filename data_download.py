# %% Imports
import logging
import subprocess

import ee
import xarray as xr
from icoscp_core import auth

from conf import GEE_PROJECT_ID, era5_land_path, fluxcom_x_base_path
from functions import data_logging

# Authentication GEE
ee.Authenticate()
ee.Initialize(
    project=GEE_PROJECT_ID, opt_url="https://earthengine-highvolume.googleapis.com"
)

# Authentication ICOS
auth.init_config_file()  # Only do this once per computer!

era5_land_path.mkdir(exist_ok=True, parents=True)
fluxcom_x_base_path.mkdir(exist_ok=True, parents=True)
script_path = __file__

# %% Download ERA5_Land total evaporation
data_logging(era5_land_path, script_path)
logging.info("Downloading monthly (summed) total evaporation data from ERA5-Land")
logging.info(
    "Data source: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY_AGGR#bands"
)
logging.info("Script Author: Olivier Bonte")

ds = xr.open_dataset("ee://ECMWF/ERA5_LAND/MONTHLY_AGGR", engine="ee", scale=0.1)
ds_E = ds["total_evaporation_sum"]
ds_E.attrs["units"] = "m"
ds_E.attrs["standard_name"] = "lwe_thickness_of_water_evaporation_amount"
ds_E.attrs[
    "description"
] = """Accumulated amount of water that has evaporated from the Earth's surface, 
    including a simplified representation of transpiration (from vegetation), into vapor in the air above. 
    This variable is accumulated from the beginning of the forecast to the end of the forecast step"""
# Download from 1980 until 2025
for year in range(1980, 2026):
    logging.info(f"Downloading ERA5-Land total evaporation data for {year}")
    ds_E_ = (
        ds_E.sel(time=str(year)) * -1
    )  # Convention IFS = negative values if outgoing -> change this to positive
    ds_E_.sel(time=str(year)).to_netcdf(
        era5_land_path / f"ERA5_Land_total_evaporation_monthly_{year}.nc"
    )

# %% FLUXCOM X X-BASE download
data_logging(fluxcom_x_base_path, script_path)
logging.info("Downloading FLUXCOM-X X-BASE ET data")
logging.info(
    "Data sources: https://gitlab.gwdg.de/fluxcom/fluxcomxdata/-/blob/main/docs/01-aggregation.md"
)
logging.info("Script Author: Olivier Bonte")
subprocess.run(
    f"python external_code/download_xbase_from_icos.py ET 005_monthly -o {fluxcom_x_base_path}",
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
)
