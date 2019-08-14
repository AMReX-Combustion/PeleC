
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _LES:

LES and Hybrid LES/DNS Support
------------------------------

.. note:: PeleC has support for LES subgrid models for *non-reacting* flows; documentation is a work in progress.


Models for large eddy simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning:: The LES source terms do not currently support EB cut cells.


PeleC currently supports two LES models, the constant and dynamic
Smagorinsky models. An extensive discussion of the compressible
version of these models can be found in Martín, M. Pino, U. Piomelli,
and G. V. Candler. "Subgrid-Scale Models for Compressible Large-Eddy
Simulations." Theoretical and Computational Fluid Dynamics 13, no. 5
(2000): 361–76. The constant Smagorinsky model was verified using the
method of manufactured solutions.


Using
#####

To use the LES models, the user must first activate the LES source
terms by doing: ``pelec.do_les = 1``. The default LES model is
constant Smagorinsky. The user can pick the LES model by setting
``pelec.les_model = NUM``:

* ``les_model = 0``: constant Smagorinsky model
* ``les_model = 1``: dynamic Smagorinsky model

For the constant Smagorinsky model, the user may define the model
coefficients: ``pelec.Cs``, ``pelec.CI``, and ``pelec.PrT``. These
must be set by the user in the input file as they default to ``Cs =
0``, ``CI = 0``, and ``PrT = 1.0``. To use a typical constant
Smagorinsky model, the user would add the following in the input file:

::

   pelec.do_les = 1
   pelec.les_model = 0
   # TURBULENCE PARAMETERS
   pelec.Cs = 0.16
   pelec.CI = 0.09
   pelec.PrT = 0.7


For the dynamic Smagorinsky model, the user can define a test filter
type (``les_test_filter_type``) and filter-grid ratio
(``les_test_filter_fgr``). A complete description of the filter types
and filter-grid ratios can be found in the section on filtering
below. The default test filter type is a box filter with a filter-grid
ratio of 2. An example for activating the dynamic Smagorinsky model is:

::

   pelec.do_les = 1
   pelec.les_model = 1
   pelec.les_test_filter_type = 3
   pelec.les_test_filter_fgr = 2


Developing
##########

The LES models are implemented as source terms in `PeleC_les.cpp` and
the kernels are implemented in ``Src_DIMd/lesterm_DIMd.f90``.


Explicit filtering of the hydrodynamic source terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning:: Explicit filtering is a work in progress. It should not be used with wall boundary conditions. Filtering at coarse-fine boundaries is experimental.

In some instances the user may want to explicitly filter the
Navier-Stokes equations to perform explicitly filtered large eddy
simulations. The objective of explicitly filtering the PDEs is to
remove wavenumber content in the solution that may be contaminated by
numerical error. To achieve this objective, the hydrodynamic source
terms, i.e. the source terms that may introduce wavenumbers that
cannot be represented accurately by the grid, are explicitly filtered
after each computation of the hydrodynamic source terms.

Using
#####

To explicitly filter the hydrodynamic source terms, the following
should option should be turned on in the input file:
``pelec.use_explicit_filter = 1``. The user specifies the filter-grid
ratio using ``pelec.les_filter_fgr = NUM``, where ``NUM`` is the
filter-grid ratio desired, e.g. ``pelec.les_filter_fgr = 2``. The user
also specifies a filter type through ``pelec.les_filter_type = NUM``:

* ``les_filter_type = 0``: no filtering
* ``les_filter_type = 1``: standard box filter
* ``les_filter_type = 2``: standard Gaussian filter

We have also implemented a set of filters defined in Sagaut & Grohens (1999) Int. J. Num. Meth. Fluids:

* ``les_filter_type = 3``: 3 point box filter approximation (Eq. 26)
* ``les_filter_type = 4``: 5 point box filter approximation (Eq. 27)
* ``les_filter_type = 5``: 3 point box filter optimized approximation (Table 1)
* ``les_filter_type = 6``: 5 point box filter optimized approximation (Table 1)
* ``les_filter_type = 7``: 3 point Gaussian filter approximation
* ``les_filter_type = 8``: 5 point Gaussian filter approximation (Eq. 29)
* ``les_filter_type = 9``: 3 point Gaussian filter optimized approximation (Table 1)
* ``les_filter_type = 10``: 5 point Gaussian filter optimized approximation (Table 1)

An example input file section for a Gaussian filter with a filter-grid
ration of 2 would be:

::

   pelec.use_explicit_filter=1
   pelec.les_filter_type=2
   pelec.les_filter_fgr=2


Developing
##########

The weights for these filters are set in ``Filter.cpp``. To add a
filter type, one needs to add an enum to the ``filter_types`` and
define a corresponding ``set_NAME_weights`` function to be called at
initialization.

The application of a filter can be done on a Fab or MultiFab. The
implementation is done in ``filter_DIMd.f90``. The loop nesting
ordering was chosen to be performant on existing HPC architectures and
discussed in PeleC milestone reports. An example call to the filtering operation is

::

   les_filter = Filter(les_filter_type, les_filter_fgr);
   ...
   les_filter.apply_filter(bxtmp, flux[i], filtered_flux[i], Density, NUM_STATE);

The user must ensure that the correct number of grow cells is present in the Fab or MultiFab.
