# State of the Climate reports

The generate the figures of land evaporation for the State of the Climate (SotC) reports (see e.g. the 2023 report [here](https://doi.org/10.1175/2024BAMSStateoftheClimate.1))

- [`data_download_and_process_sotc.py`](data_download_and_process_sotc.py) script is provided to download Southern Oscillation Index (SOI) data (see [here](https://crudata.uea.ac.uk/cru/data/soi/soi_3dp.dat))and process GLEAM4 data (see the `# GLEAM 4 Paths` in [`conf.py`](../conf.py)). Also see the `# SotC Paths` in [`conf.py`](../conf.py) to set the output paths.
- [`plots_sotc.py`](plots_sotc.py) script generates the figures. Note that in [`conf.py`](../conf.py) you can set figure settings under `# %% SotC options`:
  - Fig x: Yearly anomalies time series
  - Fig y: Hovmöller diagram monthly anomalies
  - Plate x: Spatial anomaly map
