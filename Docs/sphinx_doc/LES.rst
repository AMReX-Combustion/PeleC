
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _LES:

LES and Hybrid LES/DNS Support
------------------------------

.. note:: PeleC has support for LES subgrid models for *non-reacting* flows; documentation is a work in progress.

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
