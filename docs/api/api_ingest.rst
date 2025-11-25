.. |pd.DataFrame| replace:: pandas.DataFrame
.. |np.ndarray[int]| replace:: numpy.ndarray[int]
.. |np.ndarray[float]| replace:: numpy.ndarray[float]
.. |np.ndarray[np.number]| replace:: numpy.ndarray[numpy.number]


Data ingestion
==============

Acoustics
---------
.. autofunction:: echopop.ingest.sv.ingest_echoview_sv
.. autofunction:: echopop.ingest.nasc.merge_echoview_nasc
.. autofunction:: echopop.ingest.nasc.consolidate_echvoiew_nasc
.. autofunction:: echopop.ingest.nasc.generate_transect_region_haul_key
.. autofunction:: echopop.ingest.nasc.process_region_names
.. autofunction:: echopop.ingest.nasc.read_nasc_file
.. autofunction:: echopop.ingest.nasc.read_transect_region_haul_key

Biological
----------
.. autofunction:: echopop.ingest.load_biological_data

Strata and Geostrata
--------------------
.. autofunction:: echopop.ingest.load_geostrata
.. autofunction:: echopop.ingest.join_geostrata_by_latitude
.. autofunction:: echopop.ingest.load_strata
.. autofunction:: echopop.ingest.join_strata_by_haul

Spatial
-------
.. autofunction:: echopop.ingest.load_isobath_data
.. autofunction:: echopop.ingest.load_mesh_data
.. autofunction:: echopop.ingest.load_kriging_variogram_params