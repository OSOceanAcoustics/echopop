.. |pd.DataFrame| replace:: pandas.DataFrame
.. |gpd.GeoDataFrame| replace:: geopandas.GeoDataFrame
.. |np.ndarray[int]| replace:: numpy.ndarray[int]
.. |np.ndarray[float]| replace:: numpy.ndarray[float]
.. |np.ndarray[np.number]| replace:: numpy.ndarray[numpy.number]

Geostatistics
=============

Mesh Cropping and Coordinate Transformations
--------------------------------------------
.. autofunction:: echopop.geostatistics.hull_crop
.. autofunction:: echopop.geostatistics.transect_coordinate_centroid
.. autofunction:: echopop.geostatistics.transform_coordinates
.. autofunction:: echopop.geostatistics.transect_extent
.. autofunction:: echopop.geostatistics.utm_string_generator
.. autofunction:: echopop.geostatistics.wgs84_to_utm
.. autofunction:: echopop.geostatistics.projection.reproject_dataset

(Semi)Variogram and Covariance Models
-------------------------------------
.. autoclass:: echopop.geostatistics.Variogram
   :members: 
.. autofunction:: echopop.geostatistics.compute_variogram
.. autofunction:: echopop.geostatistics.fit_variogram   

Kriging and Spatial Interpolation
---------------------------------
.. autoclass:: echopop.geostatistics.Kriging
   :members:
.. autofunction:: echopop.geostatistics.uniform_search_strategy   