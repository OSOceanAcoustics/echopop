# Installation

Echopop is not yet available to be installed on PyPI or conda. We plan to enable these distributions at release v0.4.2 after the current refacotring code is completed.

Until then, you can install Echopop from the repository following the steps below.

## 1. Create a virtual environment

To keep your Echopop environment isolated, it's recommended to create a virtual environment using Conda or Python's built-in `venv` module. Here's an example using Conda:

```bash
conda create --name echopop-env
conda activate echopop-env
```

```{attention} 
We recommend using the ``libmamba`` solver instead of the classic solver.
   See instructions `here <https://conda.github.io/conda-libmamba-solver/getting-started/>`_
   for installation and usage.
```


Or, using Python's venv:

```bash
python -m venv echopop-env
source echopop-env/bin/activate  # On Windows, use `echopop-env\Scripts\activate`
```

## 2. Clone the repo
Now that you have a virtual environment set up, you can clone the Echopop project repository to your local machine using the following command:

```bash
git clone https://github.com/OSOceanAcoustics/echopop.git
```

## 3. Install the package

Navigate to the project directory you've just cloned and install the Echopop package:

```bash
cd <project_directory>
pip install -e .
```

The `-e` flag here means that you are installing Echopop in a development mode, which allows you to not only use but also develop the code.
