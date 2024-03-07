# Installation

Echopop is not yet available to be installed on PyPI or conda. We plan to enable these distributions at release v0.4.2 after the current code refactoring is completed.

Until then, you can install Echopop from the repository following the steps below:

```shell
# create a conda environment using the supplied requirements files
conda create -c conda-forge -n echopop --yes python=3.9 --file requirements.txt

# switch to the newly built environment
conda activate echopop

# ipykernel is recommended, in order to use with JupyterLab and IPython
# to aid with development. We recommend you install JupyterLab separately
conda install -c conda-forge ipykernel

# install echopop in editable mode (-e)
pip install .
```

```{attention} 
We recommend using the ``libmamba`` solver instead of the classic solver.
   See instructions `here <https://conda.github.io/conda-libmamba-solver/getting-started/>`_
   for installation and usage.
```


If you want to create the full development environment so that you can make changes to the Echopop code, use below when creating the conda environment:
```shell
# note the last one docs/requirements.txt is only required for building docs
conda create -c conda-forge -n echopop --yes python=3.9 --file requirements.txt --file requirements-dev.txt --file docs/requirements.txt
```

and when installing echopop at the end, use:
```shell
# install echopop in editable mode
pip install -e .
```

The `-e` flag here means that you are installing Echopop in a development mode, which allows you to not only use but also develop the code.
