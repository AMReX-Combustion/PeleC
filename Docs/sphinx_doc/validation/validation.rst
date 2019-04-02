
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Validation:


.. highlight:: rst

Validation of PeleC
-------------------

Decay of homogeneous isotropic turbulence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _HIT:

Simulations were performed for turbulent Mach number = 0.1, Taylor
scale based Reynolds number = 100, Prandtl number = 0.71, and k0 = 4
at several different resolutions (uniform discretization, no AMR)
using a prescribed energy spectrum.

The definitions for the different quantities and reference data (in
black) can be found in `Johnsen et al. (2009) JCP
<http://dx.doi.org/10.1016/j.jcp.2009.10.028>`_ and `Movahed and
Johnsen (2015) JFM <http://dx.doi.org/10.1017/jfm.2015.200>`_. VisIt
was used to post-process some quantities using `visit_pp_aux_vars.py`
(located in the Exec/HIT folder) with the command

.. code-block:: bash

   $ visit -nowin -cli -s Exect/RegTests/visit_pp_aux_vars.py

Generating the initial conditions
#################################

The initial conditions for this validation problem are derived from
the following energy spectrum:

.. math::
   E(k) \sim k^4 \exp(-2 (k/k_0)^2), \frac{3 u'}{2} = \int_0^\infty E(k) \mathrm{d}k

and can be generated with the Python script in the HIT directory:

.. code-block:: bash

   $ python3 gen_hit_ic.py

Building and running
####################

The decay of homogeneous isotropic turbulence case can be found in ``Exec/RegTests/HIT``:

.. code-block:: bash

   $ make -j 16 DIM=3 USE_MPI=TRUE
   $ mpiexec -n 36 $EXECUTABLE inputs_3d

The user can run a convergence study by generating initial conditions
for higher resolutions and varying ``amr.ncell``.


Results
#######

As the resolution increases, there is good agreement between the Pele
data and reference data.

.. figure:: ./validation/hit/reynolds.png
   :align: center
   :figwidth: 40%

   Reynolds number as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: Johnsen et al. (2009) JCP.

.. figure:: ./validation/hit/mach.png
   :align: center
   :figwidth: 40%

   Mach number as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: Johnsen et al. (2009) JCP.

.. figure:: ./validation/hit/enstrophy.png
   :align: center
   :figwidth: 40%

   Enstrophy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: Johnsen et al. (2009) JCP.

.. figure:: ./validation/hit/KE.png
   :align: center
   :figwidth: 40%

   Kinetic energy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`.


Taylor-Green vortex breakdown
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This setup is one of the test problems outlined by `the High-Order CFD
workshop <https://www.grc.nasa.gov/hiocfd>`_. A complete description
of the problem can be found `at NASA HOCFDW website
<https://www.grc.nasa.gov/hiocfd/wp-content/uploads/sites/22/case_c3.3.pdf>`_
and the reference data is found `here
<https://www.grc.nasa.gov/wp-content/uploads/sites/22/C3.3_datafiles.zip>`_. More
details of the problem and methods used to obtain the reference data
can be found in `Bull and Jameson (2014) 7th AIAA Theoretical Fluid
Mechanics Conference (doi: 10.2514/6.2014-3210)` and `DeBonis (2013)
51st AIAA Aerospace Sciences Meeting (doi:10.2514/6.2013-382)`.

Building and running
####################

The Taylor-Green vortex breakdown case can be found in ``Exec/RegTests/TG``:

.. code-block:: bash

   $ make -j 16 DIM=3 USE_MPI=TRUE
   $ mpiexec -n 36 $EXECUTABLE inputs_3d amr.ncell=64 64 64

The user can run a convergence study by varying ``amr.ncell``.


Results
#######

As the resolution increases, there is good agreement between the Pele
data and reference data.

.. figure:: ./validation/tg/dissipation.png
   :align: center
   :figwidth: 40%

   Dissipation as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: HOCFDW.

.. figure:: ./validation/tg/enstrophy.png
   :align: center
   :figwidth: 40%

   Enstrophy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black squares: HOCFDW.

.. figure:: ./validation/tg/KE.png
   :align: center
   :figwidth: 40%

   Kinetic energy as a function of time. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black: HOCFDW.

.. figure:: ./validation/tg/spectrum.png
   :align: center
   :figwidth: 40%

   Spectrum at :math:`t=9 t_c`. Solid red: :math:`32^3`, dashed green :math:`64^3`, dash-dotted blue: :math:`128^3`, dotted orange: :math:`256^3`, black: HOCFDW.
