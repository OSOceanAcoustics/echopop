.. |pd.DataFrame| replace:: pandas.DataFrame
.. |np.ndarray[int]| replace:: numpy.ndarray[int]
.. |np.ndarray[float]| replace:: numpy.ndarray[float]
.. |np.ndarray[np.number]| replace:: numpy.ndarray[numpy.number]

Survey
======

Biological
----------
.. autofunction:: echopop.survey.biology.fit_length_weight_regression
.. autofunction:: echopop.survey.biology.quantize_length_data

Distributions and proportions
-----------------------------
.. autofunction:: echopop.survey.proportions.compute_binned_counts
.. autofunction:: echopop.survey.proportions.number_proportions
.. autofunction:: echopop.survey.proportions.binned_weights
.. autofunction:: echopop.survey.proportions.stratum_averaged_weight
.. autofunction:: echopop.survey.proportions.weight_proportions
.. autofunction:: echopop.survey.proportions.fitted_weight_proportions
.. autofunction:: echopop.survey.proportions.fitted_weight_proportions_combined
.. autofunction:: echopop.survey.proportions.get_nasc_proportions_slice
.. autofunction:: echopop.survey.proportions.get_number_proportions_slice
.. autofunction:: echopop.survey.proportions.get_weight_proportions_slice

Net selectivity correction
--------------------------
.. autofunction:: echopop.survey.selectivity.assign_selectivity_expansion

Statistics
----------
.. autofunction:: echopop.survey.statistics.confidence_interval
.. autoclass:: echopop.survey.stratified.JollyHampton
   :members: 

Transect processing
-------------------
.. autofunction:: echopop.survey.transect.compute_interval_distance