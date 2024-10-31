# GLEAM4

Repository for reproducing the figures of the GLEAM4 paper.

The required dependencies to run the python code in this repository is given in `environment.yml`. Install via [(Mini)conda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/) (in this case, replace `conda` by `mamba` in the command below) in the command line interface:
```
conda env create -f environment.yml
conda activate gleam4_paper
```

## Selection of years for timeseries

- US-Ne3: 2012 is chosen because of drought conditions (see [here](https://www.drought.gov/states/nebraska#historical-conditions))
- DE-Tha: 2003 is chosen because of exceptionally hot summer, which had a larger influence on ecosystem function (see for example [Ciais et al., 2005](https://doi.org/10.1038/nature03972))
- AU-How: 2010 as an exceptionally wet year, see [here](https://en.wikipedia.org/wiki/2000s_Australian_drought#2010_and_2011:_La_Ni%C3%B1a_finally_breaks_the_drought)
- FR-Pue: 2016 as a drought year (see [here](https://doi.org/10.1175/JCLI-D-18-0331.1))