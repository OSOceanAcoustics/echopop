API reference
==============

**Content**

* `Survey`_
* `Data loading`_
    * :py:class:`KrigingMesh <echopop.data_loader.KrigingMesh>`
    * :py:class:`LoadBioData <echopop.data_loader.LoadBioData>`
    * :py:class:`LoadStrataData <echopop.data_loader.LoadStrataData>`
    * :py:class:`load_nasc_df <echopop.data_loader.load_nasc_df>`
* `Computational routines`_
    * :py:class:`Kriging <echopop.computation.Kriging>`
    * :py:class:`SemiVariogram <echopop.computation.SemiVariogram>`
    * :py:class:`ComputeTransectVariables <echopop.computation.ComputeTransectVariables>`
    * :py:class:`ComputeKrigingVariables <echopop.computation.ComputeKrigingVariables>`
    * :py:class:`Bootstrapping <echopop.computation.Bootstrapping>`
* `Reports`_
    * :py:class:`Reports <echopop.reports.Reports>`


Survey
------

.. automodule:: echopop
   :members: Survey


Data loading
------------

.. automodule:: echopop.data_loader
    :members: KrigingMesh, LoadBioData, LoadStrataData, load_nasc_df


Computational routines
----------------------

Kriging
^^^^^^^

.. automodule:: echopop.computation
    :members: Kriging, SemiVariogram

Computation of results
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: echopop.computation
    :members: ComputeTransectVariables, ComputeKrigingVariables

Bootstrapping
^^^^^^^^^^^^^

.. automodule:: echopop.computation
    :members: Bootstrapping


Reports
-------

.. automodule:: echopop.reports
    :members: Reports
