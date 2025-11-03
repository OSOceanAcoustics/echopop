# Example notebooks

A variety of analyses will be available via `echopop`, including the example workflows provided here. It may be helpful to assign a `DATA_ROOT` directory path as a shortcut for writing the full filepaths of target datasets. However, this is dependent on how the different acoustic, biological, and other datasets are organized. 

```python
from pathlib import Path

DATA_ROOT = Path("C:/Data/")
```

`````{admonition} Function arguments
:class: tip
Use the `help()` function or hover tooltips in certain IDEs (e.g. `VSCode`) to investigate how certain functions are parameterized and used.
`````

The workflows have been broken up into the following sections:

- [](nasc-ingestion)