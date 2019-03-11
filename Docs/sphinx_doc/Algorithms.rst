
.. _GettingStarted:


.. Warning:: This documentation is a placeholder, and contents should currently be considered a work in progress, out of context, or just plain wrong until this note is removed!

Algorithms
==========

PeleC Timestepping
------------------

PeleC supports two options for timestepping: a second-order explicit method-of-lines approach (MOL), and an iterative scheme base on a spectral deferred correction approach (SDC). Both time-steppers share considerable code in order to actually perform the update.


Standard Time Advance
~~~~~~~~~~~~~~~~~~~~~
The MOL time stepper is a standard second order predictor-corrector approach with (optional) fixed point iteration to tightly couple the reaction and transport. The advection :math:`(A)` and diffusion :math:`(D)` terms are computed explicitly using a time-explicit finite-volume formulation; reaction trerms are integrated with VODE (DVODE or CVODE via SUNDIALS), with a forcing term that incorporates the (pointwise) influence of advection and diffusion (:math:`F_{AD}`). For cold start the reaction term (:math:`I_R`) is evaluated with (:math:`F_{AD} = 0`)

.. math::
   S^n &= AD(u^n) \hspace{2em} {\small \text{(stencils require grow-cell data at }t^{n}\text{)}}

   u^* &= u^n + \Delta t(S^n +I_R)

   S^{n+1} &= AD(u^*) \hspace{2em} {\small \text{(stencils require grow-cell data at }t^{n+1}\text{)}}

   u^{**} &= \frac{1}{2}(u^n+u^*) + \frac{1}{2}\left(S^{n+1}+I_R\right){\Delta t}

   F_{AD} &= \frac{1}{\Delta t} (u^{**} -u^n) - I_R

   I_R &= I_R(u^n, F_{AD})

   u^{n+1} &= u^n + \Delta t(F_{AD} +I_R)\text{.}


Note that the advection and diffusion terms (evaluated above at old and new time) require grow cells to be filled at the appropriate solution time.  The filling operation is orchestrated by the AMReX software framework via the `FillPatch` operation.  Grow cells from neighboring mesh patches (and through periodic/re-entrant boundaries) are copied on intersection.  User-specified functions provide data at the physical boundaries.  Cells along the coarse-fine boundary are interpolated in space and time from available coarse data (note that this requires that the fine data be "properly nested" in the coarser levels).  Also, because :math:`(A)` and :math:`(D)` are both time-explicit, they are computed together using grow cells filled by the same `FillPatch` operation.

With time-implicit reactions, the final update is iterated:

.. math::
   S^{n+1,k} &= AD(u^{n+1,k})

   F_{AD}^{k} &= \frac{1}{2}(S^n+S^{n+1,k})

   I_R^{k} &= I_R(u^n, F_{AD}^{k})

   u^{n+1,k+1} &= u^n + \Delta t(F_{AD}^{k} +I_R^{k})\text{.}


Hyperbolics
-----------

Two hyperbolic treatments are available.

PPM
~~~

The unsplit piecewise parabolic method is used for regular geometries and is the same algorithm used in several other AMReX codes including CASTRO and MAUI. 


Method of Lines with Characteristic Extrapolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative formulation well suited to Embedded Boundary geometry treatment and also avaialble for regular grids is available and based on a method of lines approach. 


Diffusion
---------

Diffusion is modeled with a mixture average formulation.

Reaction
--------

A chemical reaction network is evaluated to determine the reaction source term.  The reaction network is selected at build time by setting the `CHEMISTRY_MODEL` flag in the makefile, where the value refers to one of the models available in `PelePhysics`. New models can be generated using `Fuego`, currently not part of `PelePhysics` but slated for inclusion in the near future.


Equation of State
-----------------

Several equation of state models are available based on ideal gas, gamma law gas or non-ideal equation of state. 
