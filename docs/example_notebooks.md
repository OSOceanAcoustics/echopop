(user_guide)=
# User guide

The notebooks in this section demonstrate how to use Echopop functions to string together a workflow to infer biomass and other population variables using the combination of acoustic and biological data.  

It is often helpful to assign a `DATA_ROOT` directory path as a shortcut for writing the full filepaths of target datasets, as shown in the example. However, this depends on how the different input files are organized.  

```python
from pathlib import Path

DATA_ROOT = Path("C:/Data/")
```

`````{admonition} Function arguments
:class: tip
Use the `help()` function or hover tooltips in certain IDEs (e.g. `VSCode`) to investigate the arguments and scope of each function.
`````

The workflows have been broken up into the following sections:

- [](nasc-ingestion)
- [](biodata-ingestion)
- [](stratify-data)
- [](number-proportions)
- [](weight-proportions)
- [](basic-inversion)
- [](geostats)
- [](kriged-biomass-to-nasc)
- [](jolly-hampton)
- [](basic-plotting)
- [](report-generation)
- [](year-specific-workflows)
