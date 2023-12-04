API reference
==============

**Content**

* `Survey`_
* `Data loading`_
    * :py:class:`KrigingMesh <EchoPro.data_loader.KrigingMesh>`
    * :py:class:`LoadBioData <EchoPro.data_loader.LoadBioData>`
    * :py:class:`LoadStrataData <EchoPro.data_loader.LoadStrataData>`
    * :py:class:`load_nasc_df <EchoPro.data_loader.load_nasc_df>`
* `Computational routines`_
    * :py:class:`Kriging <EchoPro.computation.Kriging>`
    * :py:class:`SemiVariogram <EchoPro.computation.SemiVariogram>`
    * :py:class:`ComputeTransectVariables <EchoPro.computation.ComputeTransectVariables>`
    * :py:class:`ComputeKrigingVariables <EchoPro.computation.ComputeKrigingVariables>`
    * :py:class:`Bootstrapping <EchoPro.computation.Bootstrapping>`
* `Reports`_
    * :py:class:`Reports <EchoPro.reports.Reports>`


Survey
------

.. automodule:: EchoPro
   :members: Survey


Data loading
------------

.. automodule:: EchoPro.data_loader
    :members: KrigingMesh, LoadBioData, LoadStrataData, load_nasc_df


Computational routines
----------------------

Kriging
^^^^^^^

.. automodule:: EchoPro.computation
    :members: Kriging, SemiVariogram

Computation of results
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: EchoPro.computation
    :members: ComputeTransectVariables, ComputeKrigingVariables

Bootstrapping
^^^^^^^^^^^^^

.. automodule:: EchoPro.computation
    :members: Bootstrapping


Reports
-------

.. automodule:: EchoPro.reports
    :members: Reports
