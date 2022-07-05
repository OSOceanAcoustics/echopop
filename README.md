# EchoPro

This repository is for the development of the Python incarnation of EchoPro. This program uses combined acoustic 
data analysis results with biological information from trawls (such as length, age, etc.) to produce biomass estimates 
of Pacific hake.

## Python EchoPro Workflow

We have compiled the Jupyter notebook [echopro_workflow.ipynb](https://github.com/uw-echospace/EchoPro/blob/master/echopro_workflow.ipynb) 
that highlights the current workflow of the Python version of EchoPro. It can be ran after cloning this repository and
setting up the Jupyter notebook. Below we outline the process of setting up the notebook using [Anaconda](https://www.anaconda.com/) 
and obtaining all necessary data. All steps contained in codeblocks should be performed using a terminal. 

1. Use [conda_install.yaml](https://github.com/uw-echospace/EchoPro/blob/master/EchoPro/conda_install.yaml)
to create the conda environment called `echopro_env`:
```
conda env create -f conda_install.yaml 
```
2. Activate the conda environment:
```
conda activate echopro_env
```
3. Construct a kernel for the Jupyter notebook called `echopro_env` using our conda environment:
```
python -m ipykernel install --user --name=echopro_env
```
4. Navigate to the directory that contains `echopro_workflow.ipynb` and open it:
```
jupyter-notebook echopro_workflow.ipynb 
```
5. Once the Jupyter notebook has been opened, you can set the kernel by going to the "Kernel" tab in the notebook, 
hovering over "Change kernel", and selecting "echopro_env".
6. Download the folder [2019_consolidated_files](https://drive.google.com/drive/folders/13o1z5ebn3G05kAmfAVYJ3QqNEgxL8xxw?usp=sharing), 
which contains all files necessary to run the notebook. Note that this link has restricted access and the folder can 
only be downloaded by approved parties. 
7. Modify `data_root_dir` in [survey_year_2019_config.yml](https://github.com/uw-echospace/EchoPro/blob/master/config_files/survey_year_2019_config.yml)
so that it points to the directory `2019_consolidated_files`, which was downloaded in the previous step.

## Project updates

* [April 7, 2022](https://github.com/uw-echospace/EchoPro/blob/master/project_docs/2022_04_07_update.md)


