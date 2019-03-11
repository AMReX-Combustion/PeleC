
.. _GettingStarted:


.. Warning:: This documentation is a placeholder, and contents should currently be considered a work in progress, out of context, or just plain wrong until this note is removed!

Algorithms
==========

PeleC Timestepping
------------------

PeleC uses two method of timestepping: a second order explicit method, and a spectral deferred correction (SDC) approach. These approaches share several code modules to perform the update; both used an iteration to couple the various physics together.


Standard Time Advance
~~~~~~~~~~~~~~~~~~~~~
The standard time advance is a second order predictor-corrector approach with (optional) fixed point iteration to tightly couple the reaction and transport. The Advection and diffusion (:math:`AD`) terms are computed explicitly using a finite-volume formulation; reaction terms are integrated with VODE (SUNDIALS) with a forcing term that includes advection and diffusion (:math:`F_{AD}`). For cold start the reaction term (:math:`I_R`) is evaluated from the instantaneous state without a forcing term.

.. math::
   S^n = AD(\overbrace{u^n}^\text{FillPatch at $t^n$})

   u^* = u^n + dt(S^n +I_R)

   S^{n+1}= AD(\overbrace{u^*}^\text{FillPatch at $t^{n+1}$})

   u^{**} = \frac{1}{2}(u^n+u^*) + \frac{1}{2}\left(S^{n+1}+I_R\right){dt}

   F_{AD} = \frac{1}{dt} (u^{**} -u^n) - I_R

   \text{update } I_R(u^n, F_{AD}) \text{ and }  u^{n+1} = u^n + dt(F_{AD} +I_R)\text{.}


Without reaction this would be the end of the timestep; with reaction, we iterate :math:`mol\_iters` times:

.. math::
   S^{n+1}= AD(\overbrace{u^{n+1}}^\text{FillPatch at $t^{n+1}$})

   F_{AD} = \frac{1}{2}(S^n+S^{n+1})

   \text{update } I_R(u^n, F_{AD}) \text{ and }  u^{n+1} = u^n + dt(F_{AD} +I_R)\text{.}


Hyperbolics
-----------

Two hyperbolic treatments are available.

PPM
~~~

The unsplit piecewise parabolic method is used for regular geometries and is the same algorithm used in several other AMReX codes including CASTRO and MAUI. 


Method of Lines with Characteristic Extrapolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative formulation well suited to Embedded Boundary geometry treatment and also available for regular grids is available and based on a method of lines approach. 


Diffusion
---------

Diffusion is modeled with a mixture average formulation.

Reaction
--------

A chemical reaction network is evaluated to determine the reaction source term.  The reaction network is selected at build time by setting the `CHEMISTRY_MODEL` flag in the makefile, where the value refers to one of the models available in `PelePhysics`. New models can be generated using `Fuego`, currently not part of `PelePhysics` but slated for inclusion in the near future.


Equation of State
-----------------

Several equation of state models are available based on ideal gas, gamma law gas or non-ideal equation of state. 
