
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Equations:



Equations
=========

Conservative system
-------------------

PeleC advances the following set of fully compressible equations for the conserved state vector: :math:`\mathbf{U} = (\rho, \rho \mathbf{u}, \rho E, \rho Y_k, \rho A_k, \rho B_k):`

.. math::
 
  \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u}) + S_{{\rm ext},\rho},

  \frac{\partial (\rho \mathbf{u})}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} \mathbf{u} -p{\mathcal I} + \Pi) - \nabla p +\rho \mathbf{g} + \mathbf{S}_{{\rm ext},\rho\mathbf{u}},

  \frac{\partial (\rho E)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} E + p \mathbf{u}) + \rho \mathbf{u} \cdot \mathbf{g} + \nabla\cdot \boldsymbol{\mathcal{Q}}+ S_{{\rm ext},\rho E},

  \frac{\partial (\rho Y_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} Y_k)
  - \nabla \cdot \boldsymbol{\mathcal{F}}_{k} + \rho \dot\omega_k + S_{{\rm ext},\rho Y_k},

  \frac{\partial (\rho A_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} A_k) + S_{{\rm ext},\rho A_k},

  \frac{\partial (\rho B_k)}{\partial t} &=& - \nabla \cdot (\rho \mathbf{u} B_k) + S_{{\rm ext},\rho B_k}.
  


Here :math:`\rho, \mathbf{u}, T`, and :math:`p` are the density, velocity,
temperature and pressure, respectively. :math:`E
= e + \mathbf{u} \cdot \mathbf{u} / 2` is the total energy with :math:`e` representing the
internal energy, which is defined as in the CHEMKIN standard to include both sensible
and chemical energy (species heats of formation) and is conserved across chemical reactions.
:math:`Y_k` is the mass fraction of the :math:`k^{\rm th}` species,
with associated production rate, :math:`\dot\omega_k`.  Here :math:`\mathbf{g}` is the gravitational vector, and
:math:`S_{{\rm ext},\rho}, \mathbf{S}_{{\rm ext},\rho\mathbf{u}}`, etc., are user-specified
source terms.  :math:`A_k` is an advected quantity, i.e., a tracer.  Also
:math:`\boldsymbol{\mathcal{F}}_{m}, \Pi`, and :math:`\boldsymbol{\mathcal{Q}}` are
the diffusive transport fluxes for species, momentum and heat.  Note that the internal
energy for species :math:`k` includes its heat of formation (and can therefore take on negative and
positive values).  The auxiliary fields, :math:`B_k`, have a user-defined
evolution equation, but by default are treated as advected quantities.

In the code we carry around :math:`T` and :math:`\rho e` in the
state vector even though they are redundant with the state since they may be derived from the other conserved
quantities.  The ordering of the elements within :math:`\mathbf{U}` is defined
by integer variables in the routine ``set_method_params`` in ``Src_nd/PeleC_nd.F90``.

Some notes:

* Regardless of the dimensionality of the problem, we always carry
  all 3 components of the velocity. You should always initialize all velocity components to zero, and
  always construct the kinetic energy with all three velocity components.

* There are ``NADV`` advected quantities, which range from :math:`{\tt
  UFA: UFA+nadv-1}`.  The advected quantities have no effect at all on
  the rest of the solution but can be useful as tracer quantities.

* There are ``NSPECIES`` species defined in the chemistry model, which range from :math:`{\tt UFS: UFS+nspecies-1}`.

* There are ``NAUX`` auxiliary variables, from :math:`{\tt UFX:UFX+naux-1}`. The auxiliary variables are passed into the equation
  of state routines along with the species.



Primitive Forms
---------------

PeleC uses the primitive form of the fluid equations, defined in terms of
the state :math:`\mathbf{Q} = (\rho, \mathbf{u}, p, \rho e, Y_k, A_k, B_k)`, to construct the
interface states that are input to the Riemann problem. All of the primitive variables are derived from the conservative state
vector. This task is performed in the routine ``ctoprim`` located in ``Src_nd/advection_util_nd.F90``.

The primitive variable equations for density, velocity, and pressure are:

.. math::
  
  \frac{\partial\rho}{\partial t} &=& -\mathbf{u}\cdot\nabla\rho - \rho\nabla\cdot\mathbf{u} + S_{{\rm ext},\rho}

  \frac{\partial\mathbf{u}}{\partial t} &=& -\mathbf{u}\cdot\nabla\mathbf{u} - \frac{1}{\rho}\nabla p + \mathbf{g} + 
  \frac{1}{\rho} (\mathbf{S}_{{\rm ext},\rho\mathbf{u}} - \mathbf{u} \; S_{{\rm ext},\rho})

  \frac{\partial p}{\partial t} &=& -\mathbf{u}\cdot\nabla p - \rho c^2\nabla\cdot\mathbf{u} +
  \left(\frac{\partial p}{\partial \rho}\right)_{e,Y}S_{{\rm ext},\rho}\nonumber

  &&+\  \frac{1}{\rho}\sum_k\left(\frac{\partial p}{\partial Y_k}\right)_{\rho,e,Y_j,j\neq k}\left(\rho\dot\omega_k + S_{{\rm ext},\rho Y_k} - Y_kS_{{\rm ext},\rho}\right)\nonumber

  && +\  \frac{1}{\rho}\left(\frac{\partial p}{\partial e}\right)_{\rho,Y}\left[-eS_{{\rm ext},\rho} - \sum_k\rho q_k\dot\omega_k + \nabla\cdot k_{\rm th}\nabla T \right.\nonumber

  && \quad\qquad\qquad\qquad+\ S_{{\rm ext},\rho E} - \mathbf{u}\cdot\left(\mathbf{S}_{{\rm ext},\rho\mathbf{u}} - \frac{\mathbf{u}}{2}S_{{\rm ext},\rho}\right)\Biggr] 
  

The advected quantities appear as:

.. math::
  
  \frac{\partial Y_k}{\partial t} &=& -\mathbf{u}\cdot\nabla Y_k + \dot\omega_k + \frac{1}{\rho}
                                     ( S_{{\rm ext},\rho Y_k}  - Y_k S_{{\rm ext},\rho} ),

  \frac{\partial A_k}{\partial t} &=& -\mathbf{u}\cdot\nabla A_k + \frac{1}{\rho}
                                     ( S_{{\rm ext},\rho A_k} - A_k S_{{\rm ext},\rho} ),

  \frac{\partial B_k}{\partial t} &=& -\mathbf{u}\cdot\nabla B_k + \frac{1}{\rho} 
                                     ( S_{{\rm ext},\rho B_k}  - B_k S_{{\rm ext},\rho} ).
  


When accessing the primitive variable state vector, the integer variable
keys for the different quantities are listed in the routine ``set_method_params`` in ``Src_nd/PeleC_nd.F90``.


Source Terms
------------

We now compute explicit source terms for each variable in :math:`\mathbf{Q}` and
:math:`\mathbf{U}`.  The primitive variable source terms will be used to construct
time-centered fluxes.  The conserved variable source will be used to
advance the solution. This task is performed in the routine ``srctoprim`` located in ``Src_nd/advection_util_nd.F90``. We neglect reaction source terms since they are
accounted for in the characteristic integration in the PPM algorithm.  The source terms are:

.. math::
  
    \mathbf{S}_{\mathbf{Q}}^n =
    \left(\begin{array}{c}
    S_\rho \\
    S_{\mathbf{u}} \\
    S_p \\
    S_{\rho e} \\
    S_{Y_k} \\
    S_{A_k} \\
    S_{B_k}
    \end{array}\right)^n  =  \left(\begin{array}{c}  S_{{\rm ext},\rho} \\
    \mathbf{g} + \frac{1}{\rho}\mathbf{S}_{{\rm ext},\rho\mathbf{u}} \\
    \frac{1}{\rho}\frac{\partial p}{\partial e}S_{{\rm ext},\rho E} + \frac{\partial p}{\partial\rho}S_{{\rm ext}\rho} \\
    \nabla\cdot k_{\rm th} \nabla T + S_{{\rm ext},\rho E} \\
    \frac{1}{\rho}S_{{\rm ext},\rho Y_k} \\
    \frac{1}{\rho}S_{{\rm ext},\rho A_k} \\
    \frac{1}{\rho}S_{{\rm ext},\rho B_k}
    \end{array}\right)^n,
    

.. math::
    
    \mathbf{S}_{\mathbf{U}}^n =
    \left(\begin{array}{c}
    \mathbf{S}_{\rho\mathbf{u}}\\
    S_{\rho E} \\
    S_{\rho Y_k} \\
    S_{\rho A_k} \\
    S_{\rho B_k}
    \end{array}\right)^n
    =
    \left(\begin{array}{c}
    \rho \mathbf{g} + \mathbf{S}_{{\rm ext},\rho\mathbf{u}} \\
    \rho \mathbf{u} \cdot \mathbf{g} + \nabla\cdot k_{\rm th} \nabla T + S_{{\rm ext},\rho E} \\
    S_{{\rm ext},\rho Y_k} \\
    S_{{\rm ext},\rho A_k} \\
    S_{{\rm ext},\rho B_k}
    \end{array}\right)^n.
