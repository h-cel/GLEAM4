# GLEAM4 [![DOI](https://zenodo.org/badge/872264216.svg)](https://doi.org/10.5281/zenodo.14056593)

Repository for reproducing the figures of the GLEAM4 paper ([DOI](https://doi.org/10.1038/s41597-025-04610-y)): *GLEAM4: global land evaporation and soil moisture dataset at 0.1° resolution from 1980 to near present (Miralles et al., 2025)*

The required dependencies to run the python code in this repository is given in `environment.yml`. Install via [(Mini)conda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/) (in this case, replace `conda` by `mamba` in the command below) in the command line interface:
```
conda env create -f environment.yml
conda activate gleam4_paper
```

For running the MATLAB scripts, a MATLAB license is required. 

## Data

In situ data (from eddy covariance sites) combined with model data predictions at these sites can be found on [zenodo](https://doi.org/10.5281/zenodo.14054257). The data is automatically downloaded when running the scripts `figures_timeseries.py` or `figures_validation.py`. If problems were te occur, considering downloading them manually and save as `data/sites`. 

For comparison purposes, the following datasets are also used (see [downloading](data_download.py) and [processing](data_processing.py) scripts):
- GLEAM v3.8: Data access provided upon request (info@gleam.eu)
- ERA5-Land monthly: Downloaded from [Google Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY_AGGR#description)
- FLUXCOM-X X-BASE: Downloaded from ICOS with the script (see [here](external_code/download_xbase_from_icos.py)) as provided in [this gitlab repo](https://gitlab.gwdg.de/fluxcom/fluxcomxdata/-/tree/main). 

The GLEAM4 data can be accessed via https://www.gleam.eu/

NOTE: You need both a GEE and an ICOS account to be able to perform the downloads! 

## Scripts
The main scripts are:
- [`figures_global_patterns_2_3_4.py`](figures_global_patterns_2_3_4.py): Scripts with the plots assessing the global behaviour of the GLEAM4 dataset (Figure 2,3 and 4)
    - Note: [`figures_global_patterns_2_3_4.m`](figures_global_patterns_2_3_4.m) is the original MATLAB script used for the publication. The Python script is a replicate of this.
- [`figures_intercomparison_5_6.m`](figures_intercomparison_5_6.m): Scripts to generate the comparison between GLEAM4, GLEAM v3.8, ERA5-Land and FLUXCOM (Figure 5 and 6)
    - Note: For Figure 5 (global intercomparison), there is also a Python script available: [`figures_intercomparison_5.py`](figures_intercomparison_5.py) which uses the more recent [FLUXCOM-X X-BASE](https://doi.org/10.5194/bg-21-5079-2024
    ).
- [`figures_validation_7.py`](figures_validation_7.py): Script with the plots to generate the comparison to in situ eddy-covariance measurements (Figure 7).
- [`figures_timeseries_8.py`](figures_timeseries_8.py): Script with the plots to generate the 4 sites (more details below) their timeseries (Figure 8).

The following are auxiliary scripts:
- [`conf.py`](conf.py): configuration file. **As user, you have to modify some paths according to your local setup.**
- [`functions.py`](functions.py): contains functions used throughout the main scripts

### Selection of years for timeseries

- US-Ne3: 2012 is chosen because of drought conditions, see [here](https://www.drought.gov/states/nebraska#historical-conditions) and [here](https://www.dallasfed.org/research/economy/~/media/documents/research/swe/2012/swe1204c.pdf)
- DE-Tha: 2003 is chosen because of exceptionally hot summer, which had a larger influence on ecosystem function (see for example [Ciais et al., 2005](https://doi.org/10.1038/nature03972))
- AU-How: 2010 as an exceptionally wet year, see [here](https://en.wikipedia.org/wiki/2000s_Australian_drought#2010_and_2011:_La_Ni%C3%B1a_finally_breaks_the_drought) and [here](http://www.bom.gov.au/climate/annual_sum/2010/index.shtml#:~:text=The%20report%20notes%20that%202010,its%20driest%20year%20on%20record.)
- FR-Pue: 2016 as a drought year, see [Garcia-Herrera et al. (2019)](https://doi.org/10.1175/JCLI-D-18-0331.1)