# Echopop

Echopop combines acoustic data collected by echosounders with biological "ground truth" information from trawls to produce biological estimates, such as biomass and abundance. Here, "pop" stands for the animal "population."

Current the processing is configured to work with Acoustic-Trawl survey data for Pacific hake, but we will soon add components to include Pacific krill into the package. The majority of the computational implementation is applicable for other fish and zooplankton species, and we plan to expand the package for general support in the near future.


## Contributors

Echopop development is currently co-led by Wu-Jung Lee (@leewujung) and Brandyn Lucca (@brandynluca). Brandon Reyes (@b-reyes) and Emilio Mayorga (@emiliom) contributed significantly to a previous version of Echopop.



<!-- ```{admonition} Glitches with some interactive graphical elements
While the notebooks in this site are rendered, there are some glitches in the display we're still working out. In particular, an [ipywidgets](https://ipywidgets.readthedocs.io/en/stable/) interactive graphical element in the semivariogram widget doesn't display correctly. The notebooks do run correctly when executed with Jupyter Notebook ("classic", not JupyterLab).
``` -->

<!-- Go to the individual example notebooks below or in the table of content on the left.

```{tableofcontents}
``` -->


## Acknowledgement

We thank Dezhang Chu (@DezhangChu) of the NOAA Northwest Fisheries Science Center (NWFSC)
for providing the Matlab EchoPro program he developed
that many elements of Echopop is based on,
as well as his detailed consultation for implementations specific to Pacific hake.

We thank Rebecca Thomas (@rebeccathomas-NOAA),
Beth Phillips (@ElizabethMPhillips),
Alicia Billings (@aliciabillings-noaa),
and Julia Clemons (@JuliaClemons-NOAA)
of the NWFSC Fisheries Engineering and Acoustics Team (FEAT)
for continuing discussions that make Echopop better.

This project is supported by NOAA Fisheries.

```{image} images/noaa_fisheries_logo.png
:alt: NOAA_fisheries_logo
:width: 230px
```


## License

Echopop is licensed under the open source [Apache 2.0 license](https://opensource.org/licenses/Apache-2.0).

---------------

Copyright (c) 2022-2024, Echopop Developers.
