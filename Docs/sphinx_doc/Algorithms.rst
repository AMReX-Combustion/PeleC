.. highlight:: rst

 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

.. _GettingStarted:


Algorithms
==========

PeleC Timestepping
------------------

PeleC uses two method of timestepping: a second order explicit method, and a spectral defferred correction (SDC) approach. The SDC approach uses iteration to couple the various source terms (advection, diffusion, reaction, spray) together. These approaches share several code modules to perform the update. 


Spectral Deferred Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 .. note::
  	This section under development

.. math::

	\Omega_i \Rightarrow \text{source term for all components from } i^{\text{th}} \text{physics}

	S \Rightarrow \text{State}

	S(t) \Rightarrow \text{Fillpached state}

.. math::
	\Omega^\text{old} \leftarrow f(S(t))
	
.. math::
	\Omega^{\text{new}} \leftarrow 	\Omega^{\text{old}}

.. math::
	\Omega^H \leftarrow h(S(t))

.. math::
	S^\text{new} \leftarrow S^\text{old} + 0.5\Delta t \sum_i \left[\Omega_i^\text{old} + \Omega_i^\text{new} \right] + (\Delta t)\Omega^H + (\Delta t )I_R

.. math::
	\Omega^\text{new} \leftarrow f'(S(t))

.. math::
	S^\text{new} \leftarrow r( \Omega^\text{old}, \Omega^\text{new}, S^\text{old}, S^\text{new})


Reaction update:


.. math::
	(\rho e_k)^\dagger = \frac{1}{2\rho^\dagger}\rho u_i^\dagger \rho u_i^\dagger

.. math::
	\hat{\rho} = \sum_i (\rho Y_i)^\dagger

.. math::
	e = \frac{\rho^\dagger - (\rho e_k)^\dagger}{\rho^\dagger}

.. math::
	\dot{(\rho e_k)} = \frac{\left[ (\rho e)^\ddagger - (\rho e_k)^\ddagger \right] - (\hat{\rho}e)}{\Delta t}

.. math::
	\dot{(\rho Y)_{\text{ext},k}} = a_{\text{ext},k} \qquad k \in \text{species}

.. math::
	(\rho Y)_k^\ddagger, e^\ddagger, T^\ddagger \leftarrow R\left[T^\dagger, (\rho Y)^\dagger, \dot{(\rho e_k)}, \dot{(\rho Y)_{\text{ext},k}} \right]

.. math::
	\rho^\ddagger = \sum_k (\rho Y_k)^\ddagger

.. math::
	(\rho u_i)^\ddagger = (\rho u_i)^\dagger + \Delta t a_i \qquad i \in UMX,UMZ

.. math::
	(\rho e)^\ddagger = \rho^\ddagger e^\ddagger

.. math::
	(\rho E)^\ddagger = (\rho e)^\ddagger + \frac{1}{2}\rho^\ddagger (\rho u_i)^\ddagger  (\rho u_i)^\ddagger

If updating, copy to S...

.. math::
	I_{R,k} = \frac{\rho Y_k)^\ddagger - (\rho Y_l)^\dagger}{\Delta t} - a_k \qquad k \in \text{species}

.. math::
	I_{R,k} = \frac{(\rho e)^\ddagger - (\rho e)^\dagger}{\Delta t} - a_k \qquad k \in \text{UEDEN}


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
