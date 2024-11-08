# GLEAM4

Repository for reproducing the figures of the GLEAM4 paper: *GLEAM4: global land evaporation dataset at 0.1° resolution from 1980 to near present (Miralles et al., in review)*

The required dependencies to run the python code in this repository is given in `environment.yml`. Install via [(Mini)conda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/) (in this case, replace `conda` by `mamba` in the command below) in the command line interface:
```
conda env create -f environment.yml
conda activate gleam4_paper
```
## Data

Data with in situ data (from eddy covariance sites) combined with model predictions at these sites can be found on [zenodo](https://doi.org/10.5281/zenodo.14054258). The data is automatically downloaded when running the scripts `figures_timeseries.py` or `figures_validation.py`. If problems were te occur, considering downloading them manually and save as `data/sites`. 

## Scripts
The main scripts are:

- `figures_validation.py`: Script with the plots to generate the comparison to in situ eddy-covariance measurements (Figure 6).
- `figures_timeseries.py`: Script with the plots to generate the 4 sites (more details below) their timeseries. 

The following are auxiliary scripts:
- `conf.py`: configuration file 
- `functions.py`: contains functions used throughout the main scrips

## Selection of years for timeseries

- US-Ne3: 2012 is chosen because of drought conditions, see [here](https://www.drought.gov/states/nebraska#historical-conditions) and [here](https://www.dallasfed.org/research/economy/~/media/documents/research/swe/2012/swe1204c.pdf)
- DE-Tha: 2003 is chosen because of exceptionally hot summer, which had a larger influence on ecosystem function (see for example [Ciais et al., 2005](https://doi.org/10.1038/nature03972))
- AU-How: 2010 as an exceptionally wet year, see [here](https://en.wikipedia.org/wiki/2000s_Australian_drought#2010_and_2011:_La_Ni%C3%B1a_finally_breaks_the_drought) and [here](http://www.bom.gov.au/climate/annual_sum/2010/index.shtml#:~:text=The%20report%20notes%20that%202010,its%20driest%20year%20on%20record.)
- FR-Pue: 2016 as a drought year, see [Garcia-Herrera et al. (2019)](https://doi.org/10.1175/JCLI-D-18-0331.1)