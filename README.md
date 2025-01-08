# GLEAM4 [![DOI](https://zenodo.org/badge/872264216.svg)](https://doi.org/10.5281/zenodo.14056593)

Repository for reproducing the figures of the GLEAM4 paper ([preprint DOI](https://doi.org/10.21203/rs.3.rs-5488631/v1)): *GLEAM4: global land evaporation dataset at 0.1° resolution from 1980 to near present (Miralles et al., 2025)*

The required dependencies to run the python code in this repository is given in `environment.yml`. Install via [(Mini)conda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/) (in this case, replace `conda` by `mamba` in the command below) in the command line interface:
```
conda env create -f environment.yml
conda activate gleam4_paper
```

For running the MATLAB scripts, a MATLAB license is required. 

## Data

In situ data (from eddy covariance sites) combined with model data predictions at these sites can be found on [zenodo](https://doi.org/10.5281/zenodo.14054257). The data is automatically downloaded when running the scripts `figures_timeseries.py` or `figures_validation.py`. If problems were te occur, considering downloading them manually and save as `data/sites`. 

## Scripts
The main scripts are:
- `figures_global_patterns.m`: Scripts with the plots assessing the global behaviour of the GLEAM4 dataset (Figure 2,3 and 4)
- `figures_intercomparison.m`: Scripts to generate the comparison between GLEAM4, GLEAM v3.8, ERA5-Land and FLUXcoM (Figure 5)
- `figures_validation.py`: Script with the plots to generate the comparison to in situ eddy-covariance measurements (Figure 6).
- `figures_timeseries.py`: Script with the plots to generate the 4 sites (more details below) their timeseries (Figure 7). 

The following are auxiliary scripts:
- `conf.py`: configuration file 
- `functions.py`: contains functions used throughout the main scrips

### Selection of years for timeseries

- US-Ne3: 2012 is chosen because of drought conditions, see [here](https://www.drought.gov/states/nebraska#historical-conditions) and [here](https://www.dallasfed.org/research/economy/~/media/documents/research/swe/2012/swe1204c.pdf)
- DE-Tha: 2003 is chosen because of exceptionally hot summer, which had a larger influence on ecosystem function (see for example [Ciais et al., 2005](https://doi.org/10.1038/nature03972))
- AU-How: 2010 as an exceptionally wet year, see [here](https://en.wikipedia.org/wiki/2000s_Australian_drought#2010_and_2011:_La_Ni%C3%B1a_finally_breaks_the_drought) and [here](http://www.bom.gov.au/climate/annual_sum/2010/index.shtml#:~:text=The%20report%20notes%20that%202010,its%20driest%20year%20on%20record.)
- FR-Pue: 2016 as a drought year, see [Garcia-Herrera et al. (2019)](https://doi.org/10.1175/JCLI-D-18-0331.1)