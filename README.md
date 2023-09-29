# Python EchoPro

Python EchoPro ("EchoPro") combines acoustic data analysis results with biological information from trawls (such as length, age, etc.) to produce estimates of biomass, biomass density, abundance, and other characteristics for Pacific hake. It is based on the Matlab EchoPro software originally developed for the [Joint U.S.-Canada Integrated Ecosystem and Pacific Hake Acoustic Trawl Survey](https://www.fisheries.noaa.gov/west-coast/science-data/joint-us-canada-integrated-ecosystem-and-pacific-hake-acoustic-trawl-survey) by the [Fisheries Engineering and Acoustic Technologies Team](https://www.fisheries.noaa.gov/west-coast/sustainable-fisheries/fisheries-engineering-and-acoustic-technologies-team) at the NOAA Northwest Fisheries Science Center.

Go to https://uw-echospace.github.io/EchoPro/ to view Jupyter notebooks that demonstrate EchoPro functionality and typical workflows.

## Installation

Python EchoPro is not yet available for installation as a package on [PyPI](https://pypi.org/) or [conda-forge](https://conda-forge.org/). Until then, it must be installed either as a **"user"** from the GitHub repository or from the code (mainly for continued development purposes, as a **"developer"**) after "cloning" the GitHub repository using `git`. See the instructions below. Either way, we'll use [conda](https://docs.conda.io) to install EchoPro dependencies using the conda environment file [condaenvironment.yaml](https://github.com/uw-echospace/EchoPro/blob/master/condaenvironment.yaml) from the repository. Installation of EchoPro dependencies has been tested extensively with `conda`. 

There are [different ways of installing `conda`](https://oceanhackweek.org/resources/prep/conda.html#installing-conda), but we recommend the use of [Miniconda](https://docs.conda.io/en/latest/miniconda.html). `conda` can be used without administrative privileges.

### Installation as a user

This simpler installation method is recommended if you don't intend to work on developing the EchoPro code base itself. It'll install both EchoPro and its dependencies in one step.

1. Download the `condaenvironment.yaml` file. In https://github.com/uw-echospace/EchoPro/blob/master/condaenvironment.yaml, click on "Raw" (on the right) then save the file.
2. At the terminal (shell), change directory to where you've placed the `condaenvironment.yaml` file.
3. Install EchoPro and its dependencies, creating a new conda environment called "echopro": 
    ```bash
    conda env create -f condaenvironment.yaml
    ```
4. To use this new conda environment, simply activate it: 
    ```bash
    conda activate echopro
    ```
5. Install EchoPro from the GitHub repository:
    - To install the latest release using the [Python wheel](https://realpython.com/python-wheels/) file:
        ```bash
        pip install https://uw-echospace.github.io/EchoPro/EchoPro-latest-py3-none-any.whl
        ```
    - To intall the latest development version:
        ```bash
        pip install git+https://github.com/uw-echospace/EchoPro.git
        ```

In order to run EchoPro you will also need to **download the EchoPro configuration files** (see below). To run the example Jupyter notebooks for the sample 2019 inputs, you will also need to download the notebooks and the input data files.

### Installation as a developer

Follow these steps if you intend to make code contributions to EchoPro:

1. Clone the repository (alternatively, fork the repository first, then clone your fork):
    ```bash
    git clone https://github.com/uw-echospace/EchoPro.git
    ```
2. `cd` to the new `EchoPro` directory:
    ```bash
    cd EchoPro
    ```
3. Install the dependencies and create a new conda environment called "echopro": 
    ```bash
    conda env create -f condaenvironment.yaml
    ```
4. Activate the environment: 
    ```bash
    conda activate echopro
    ```
5. Install EchoPro in development mode:
    ```bash
    pip install -e .
    ```

The EchoPro configuration files and example Jupyter notebooks are available in the files you have cloned, but you will need to **download the 2019 sample input data files** (see below).

## Download 2019 sample input files

Download the folder [2019_consolidated_files](https://drive.google.com/drive/folders/13o1z5ebn3G05kAmfAVYJ3QqNEgxL8xxw?usp=sharing),
which contains all input files necessary to run the example notebooks. Note that this link has restricted access and the folder can only be downloaded by approved parties.

Note that these files incorporate modifications from the original input files provided the NWFSC FEAT team (Chu). These changes involve primarily column names, some file names, and the source of mappings between haul numbers and transects numbers.

## Running EchoPro

First, set `data_root_dir` in [survey_year_2019_config.yml](https://github.com/uw-echospace/EchoPro/blob/master/config_files/survey_year_2019_config.yml)
to the path to the directory `2019_consolidated_files`, which was downloaded in the previous step.

- For a **User** installation:
    - Download `survey_year_2019_config.yml` and `initialization_config.yml` from https://github.com/uw-echospace/EchoPro/blob/master/config_files/. After navigating to each file, click on "Raw" (on the right) then save the file. Create a directory called `config_files` to store these files.
    - Download the notebooks from https://uw-echospace.github.io/EchoPro/ and place them in a new directory named `example_notebooks`, at the same level as the `config_files` directory.
- For a **Developer** installation:
    - `survey_year_2019_config.yml` is found in the `config_files` directory.
    - The notebooks are found in the `example_notebooks` directory.

Now start Jupyter **notebook**:
```bash
cd example_notebooks
jupyter-notebook
```

Select the notebook you'd like to run. Once the notebook is open, set the "kernel" (conda environment) to "echopro" by going to the menu item `Kernel > Change kernel` and selecting "Python [conda env:echopro]".

Note that an interactive widget used for exploring and selecting the semi-variogram currently does not work correctly with **JupyterLab**. This wil be addressed in a future release.

## Tests

Python EchoPro contains several tests to continuously verify functionality as changes are implemented. These tests are found in [EchoPro/tests](https://github.com/uw-echospace/EchoPro/tree/master/EchoPro/tests). They are not currently set up to run automatically in the GitHub Actions Continuous Integration system. Instead, they must be run locally. They require a set of Excel output files from Matlab EchoPro for comparison of restults, available in [this Google Drive folder](https://drive.google.com/drive/folders/1_sUDUJY_e6M9cB3n5uibWeoukPqFP8ZO?usp=drive_link). This link has restricted access and the folder can only be downloaded by approved parties. After downloading the files, please set the file path to the containing folder in `tests/conftest.py`, [here](https://github.com/uw-echospace/EchoPro/blob/master/EchoPro/tests/conftest.py#L45).
