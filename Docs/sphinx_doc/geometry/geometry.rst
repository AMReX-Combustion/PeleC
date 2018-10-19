.. highlight:: rst



Geometry treatment in PeleC
===========================

Geometry is treated in PeleC using an embedded boundary (EB) formulation, based on datastructures and algorithmic components provided by AMReX.   In the EB formalism, geometry is represented by volumes (:math:`v_l`) and apertures (:math:`A_l^k`). See :ref:`eb_cell_fig` for an illustration where the grey area represents the region excluded from the solution domain and the arrows represent fluxes.

.. _eb_cell_fig:

.. figure:: volume.pdf
   :alt: EB Cell
   :width: 500

   Embedded boundary cell

The geometry components in AMReX are used in PeleC to implement a time-explicit integrator based on the method-of-lines.  For the advection and diffusion components of the PeleC time integrator, the time rate of change of the conserved fields, S, in cell :math:`l` can be written as 

.. math::
  \frac{dS_l}{dt} = \nabla \cdot F

where :math:`F` is the intensive flux of :math:`S` through the faces that bound the cell.

Here, we are concerned with how to compute the right-hand-side of this expression in the presence of an EB. This includes how to initialize and query the necessary AMReX-provided data structures containing the geometry information, and how to compute the PeleC-specific advection and diffusion operators.  The AMReX-EB calculation requires

1. Creation of a functional specification of the irregular geometry to embed in the uniform grid.
2. Construction of map of the (continuous) implicit representation of geometry onto the discrete mesh on all AMR levels.  This will be a large, complex, distributed data structure.
3. Communication of the subsets of this large data set to the local cores tasked with building the PeleC operators.
4. Actual construction of the diffusion and advection components of the PeleC time advance.

AMReX provides for the first 3 steps.  For step 4, PeleC constructs a number of helper classes to pre-compute and organize key components of numerical operators. These objects and functions are implemented in the following files:

* Source/PeleC.cpp
* Source/PeleC.H
* Source/PeleC_init_eb.cpp
* Source/EbStencilTypes.H
* Source/EBStencilTypes_mod.F90
* Source/Src_3d/PeleC_init_eb_3d.f90
* Source/Src_3d/Hyp_pele_MOL_3d.F90
* Source/Src_3d/slope_mol_3d_EB.f90
* Source/PeleC_diffusion.cpp

In the remainder of this chapter, the mathematical algorithms for computing the fluxes are described, and are followed by an overview of the data structures used to facilitate their calculation in the presence of an EB.

Computation of state time derivative
------------------------------------

Ultimately, the time derivative that is integrated (with RK2 presently) is the "hybrid divergence" (discussed below), which computed using diffusive and advective fluxes. These fluxes are computed in the routine PeleC_diffusion.cpp. 

Diffusive Flux Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: Placeholder for section to be written.

.. f:function:: nbrsTest_nd_module/pc_compute_tangential_vel_derivs_eb


The diffusion operator is implemented with the aid of the following functions(currently in PeleC_init_eb_3d.f90 and counterparts in _2d.f90)

.. f:function:: nbrsTest_nd_module/pc_fill_bndry_grad_stencil

.. f:function:: nbrsTest_nd_module/pc_apply_eb_boundry_flux_stencil

.. f:function:: nbrsTest_nd_module/pc_fill_bndry_grad_stencil


Advective Flux Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~

The advective (hyperbolic) fluxes computation is driven by the routine pc_hyp_mol_flux found in the file Hyp_pele_MOL_3d.F90, with call signature:

.. f:function:: hyp_advection_module/pc_hyp_mol_flux

    :p q: Input state
    :p qaux: Augmented state
    :p Ax: Apertures for X edges
    :p flux1: Flux in X direction on X edges
    :p Ay: Apertures for Y edges
    :p flux2: Flux in Y direction on Y edges
    :p Az: Apertures for Z edges
    :p flux3: Flux in Z direction on Z edges
    :p flatn: Flattening parameter (not used; passed to slope routines)
    :p V: Cell volumes
    :p D: Divergence (hyperbolic fluxes added to input divergence on output)
    :p flag: Cell type flag
    :p ebflux: Flux across EB face
    :p h: Grid spacing

Within this routine, for each direction, characteristic extrapolation is used to compute left and right states at the cell faces:

.. math::
  {u^l_\perp} = u^- + \frac{1}{2\rho^-}\left( \alpha^-_2 - \alpha^-_1\right)

  {p^l} = p^- + \frac{c}{2}\left( \alpha^-_2 +\alpha^-_1\right)

  u^l_{\parallel, 1} = v^- + \frac{1}{2} \alpha^-_3

  u^l_{\parallel, 2} = w^- + \frac{1}{2} \alpha^-_4

  \rho^l Y_k^l = Y_k^-\rho^- + \frac{1}{2c}\left[\alpha^-_{4+k} + Y_k^-\left(\alpha^-_1 + \alpha^-_2\right)\right]

  \rho^l = \sum{\rho^lY_k^l}

  Y_k^l = \frac{\rho^l Y_k^l}{\rho^l}

The right states are computed as:

.. math::
  {u^r_\perp} = u^+ - \frac{1}{2\rho^+}\left( \alpha^+_2 - \alpha^+_1\right)

  {p^r} = p^+ - \frac{c}{2}\left( \alpha^+_2 +\alpha^+_1\right)

  u^r_{\parallel, 1} = v^- - \frac{1}{2} \alpha^-_3

  u^r_{\parallel, 2} = w^- - \frac{1}{2} \alpha^-_4

  \rho^r Y_k^r = Y_k^+\rho^+ - \frac{1}{2c}\left[\alpha^+_{4+k} + Y_k^+\left(\alpha^+_1 + \alpha^+_2\right)\right]

  \rho^r = \sum{\rho^rY_k^r}

  Y_k^r = \frac{\rho^r Y_k^r}{\rho^r}

The computations in the y- and z- direction are analogous; the flux on an EB face to apply a no-slip boundary condition at a wall is somewhat different. In that case, the left and right states are taken as the state at the cell center, except for the velocity is reflected across the EB face. That is:

.. math:: 
  u^l_\perp = - u \cdot \vec{n}

  u^l_{\parallel, 1} = u^l_{\parallel_2} = 0.0

  p^l = p

  Y_k^l = Y_k

  \rho^l = \rho

and, as noted the right state is identical except for:

.. math::
  u^r_\perp = - u^l_\perp

Once the left and right states are computed, a Riemann solver (in this case one preserving the physical constraints on the intermediate state) is used to compute fluxes that are assembled into a conservative and non-conservative update for the regular and cut cells.

The characteristic extrapolation requires (slope limited) fluxes; these are found in the file slope_mol_3d_EB.f90. The call signature for the slope computation is:


.. f:function:: slope_module/slopex

    :p q: Input state
    :p flatn: Flattening coefficient (not used)
    :p qaux: Augmented state (used for sound speed)
    :p flag: Cell type flag

      
Which computes the slope routines compute (limited) slopes as:

.. math::
  \Delta_1^- = 0.5\frac{1}{c}\left(p-p^-\right) - 0.5 \rho \left( u - u^-\right)  

  \Delta_2^- = 0.5\frac{1}{c}\left(p-p^-\right) + 0.5 \rho \left( u - u^-\right)  

  \Delta_3^- = v - v^-

  \Delta_4^- = w - w^-

  \Delta^-_{k=5..nspec} = \rho Y_k - \rho^- Y_k^- - \frac{1}{c^2}Y_k \left(p-p^-\right)

If cell is irregular, or neighbor to left is irregular, :math:`\Delta^- = 0.0`.

.. math::
  \Delta_1^+ = 0.5\frac{1}{c}\left(p^+ - p\right) - 0.5\rho\left(u^+ - u\right)

  \Delta_2^+ = 0.5\frac{1}{c}\left(p^+ - p\right) + 0.5\rho\left(u^+ - u\right)

  \Delta_3^+ = v^+ - v

  \Delta_4^+ = w^+ - w

  \Delta_{5...nspc}^+ = \rho^+ Y_k^+ - \rho Y_k - \frac{Y_k}{c^2}\left(p^+ - p \right)

Again, if cell is irregular, or neighbor to right is irregular, :math:`\Delta^+ = 0.0`. Finally, the slopes are limited according to:

.. math::
  \Delta_i = \frac{1}{2}\left(\Delta_i^- + \Delta_i^+\right)


  \alpha_i^{\mathrm{lim}} = \mathrm{sign}\left\{\Delta_i \right\} \cdot \min\left\{ \Delta^{lim}_i, \left|\Delta_i \right|\right\}

where:

.. math::
  \Delta^{lim} = \left\{ \begin{aligned} {} 2 \min\left\{ |\Delta^-|,|\Delta^+|\right\} \quad& \mathrm{if} \Delta^- \cdot \Delta^+ \ge 0 \\ 0 & \quad \mathrm{otherwise}\end{aligned}\right.

The formulation of the y- and z-directions is analogous to the x-direction. 


Hybrid Divergence and Redistribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A straightforward implemention of the finite-volume advance of intensive conserved fields is numerically unstable (this is the well-known "small cell issue") due to presence of the cell volume in the denominator of the time derivative:

.. math::
  (DC)_l = \frac{1}{v_l} \sum_{k_l} \left( F_k \cdot n_k A_k \right),

where :math:`k_l` is the number of regular and cut faces surrounding cell :math:`l` and :math:`F_k` is the intensive flux at the centroid of face :math:`k`.  An alternative update takes the so-called "non-conservative" form, constructed using a weighted average of the conservative updates of neighboring cells:

.. math::
  (DNC)_l = \frac{1}{\sum_{n_l}N_n v_l} \sum_{n_l}N_n v_n (DC)_n,

where :math:`n_l` is the number of cells in the `neighborhood` of cut cell :math:`l`. :math:`N_n` takes the value of 0 or 1 depending if cell :math:`n` is connected to cell :math:`l`. While this update is numerically stable, it does not discretely conserve the field quantities.  In PeleC, we use a hybrid update strategy, a weighted average of the two that is numerically stable and "maximally conservative" locally, without violating CFL constraints based on the regular cells:

.. math::
  (HD)_l = v_l(DC)_l + (1-v_l)(DNC)_l.

In order to maintain global conservation, the difference between the hybrid divergence and conservative divergence is a correction that is distributed to neighboring cells:

.. math::
  \Delta_l^n = \frac{v_l(1-v_l)\left[(DC)_l - (DNC)_l\right]N_l^n W_l^n v_n^l}{\sum_{n_l}N_l^nW_l^nv_l^n}

In PeleC, this neighborhood is obtained by the AMReX function `get_neighbors`, which identifies all cells within a single step in each coordinate direction that is connected to cell :math:`l`. Two adjacent cells may be not connected if there is an embedded boundary section between them.

The redistribution is applied as:

.. math::
  (HD)_n^l = (HD)_n^l +  \frac{\Delta_l^n}{v_n^l},

and the hybrid divergence is integrated using RK2. 

The weights for redistribution :math:`W_l^n` can be set to any field in PeleC. We have found that setting the weights to the cell volumes is effective, while pure density weighting sometimes leads to stability issues when several very small cells share a neighborhood such as in a geometry corner.

This procedure is implemented in the `pc_fix_div_and_redistribute` routine:


.. f:function:: nbrsTest_nd_module/pc_fix_div_and_redistribute

    This performs four steps
        1. Recompute conservative divergence, DC, on cut cells...need DC in 2 grow cells for    final result
        2. Compute non-conservative and hybrid divergence, DNC and HD, and redistribution mass  dM in cut cells. We will need this in 1 grow cells (see below), so it depends on     having a conservative div in 2 grow cells
        3. Now that we finished computing HD and dM everywhere, it is safe to increment DC to   hold HD
        4. Redistribute dM - THIS REQUIRES THAT DC BE GOOD IN 1 GROW CELL

    This interpolates fluxes from face centers to the centroid of the uncovered part of the face 

    :p f0: Edge centered flux in x direction on x faces
    :p f1: Edge centered flux in y direction on y faces
    :p f2: Edge centered flux in z direction on z faces
    :p sv_ebg: Geometry information for cut cells
    :p ebflux: Flux through cut face
    :p DC: Divergence

.. f:function:: nbrsTest_nd_module/pc_apply_face_stencil
    This is used to apply a pre-filled stencil operation to face data.



Data Structures and utility functions
-------------------------------------

Several structures exist to store geometry dependent information. These are populated on creation of a new AMRLevel (described below) and stored in the PeleC object so that they are available for computation. These facilitate accessing the EB data from the fortran layer and have equivalent C++ struct and fortran types definitions so that they can be passed between the languages. The C++ struct definitions are in the file EBStencilTypes.H and the fortran type definitions are in the file EBStencilTypes_mod.F90 within the pelec_eb_stencil_types_module module. The datatypes are:

+----------------+----------------+
| C++ struct     | fortran type   |
+================+================+
| EBBoundaryGeom | eb_bndry_geom  |
+----------------+----------------+
| EBBndrySten    | eb_bndry_sten  |
+----------------+----------------+
| FaceSten       | face_sten      |
+----------------+----------------+

Routines to fill and apply these as necessary can be found in the dimension specific files in e.g. Source/Src_3d/PeleC_init_eb_3d.f90 within the `nbrsTest_nd_module` module.

These structs are defined below:


.. f:type:: pelec_eb_stencil_types_module/eb_bndry_geom

.. doxygenstruct:: EBBndryGeom
    :members:
    :undoc-members:

Similarly, two structs are used to cache boundary/face stencils


.. doxygenstruct:: EBBndrySten
    :members:
    :undoc-members:

.. f:type:: pelec_eb_stencil_types_module/eb_bndry_sten


.. doxygenstruct:: FaceSten
    :members:
    :undoc-members:
.. f:type:: pelec_eb_stencil_types_module/face_sten



Initialization
--------------

Creating an EB geometry also requires knowledge of the finest level that will be used so that geometries that 'telescope', i.e., coarser volume fractions are consistent with applying the coarsening operator to the finer volumes, can be created. To that end there is a global geometry creation step, facilitated by the `initialize_EBIS` function, as well as a step that happens when a new AMRLevel is created. The latter happens by a call to  `PeleC::initialize_eb_structs`  through `PeleC::init_eb` called from the PeleC constructor. Following construction of the geometry, the geometric information is copied into the structures described in the previous section and the various interpolation stencils are populated. 

Cartesian grid, embedded boundary (EB) methods are methods where the geometric description is formed by cutting a Cartesian mesh with surface of the geometry.  AMReX's methods to handle EB geometry information, and PeleC's treatment of the EB aware update could use many possible sources for geometric description. The necessary information is, on a per-cell basis:

* Apertures for faces intersected by cut cells,
* cut cell volumes that 'telescope', that is, volumes at a coarser level are consistent with averaging the volumes from finer levels,
* connectivity indicating which neighbor cells are connected to a given cell, and
* coordinates of cell and face centroids. 

Additionally, the algorithms ultimately require surface normals, but these can be trivially recomputed from the aperture. 

In the following subsections, we will first describe using geometry creation tools based on EB infrastructure in AMReX with origins in a fork off Chombo's \cite{Chombo} infrastructure for setting up EB calculations and then how that information is used to populate the datastructures used to compute the fluxes discussed above. 

.. include:: geometry_creation.rst


Populating PeleC specific geometric description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the geometry is created, the following are populated in the Pelc::init_eb routine:

* sv_eb_bndry_geom
* sv_eb_bndry_grad_stencil
* sv_eb_flux 
* sv_eb_bcval 
* flux_interp_stencil

At present, the geometry must be static, so the above structures are valid for the life of the PeleC AMRLevel object. 

The relevant functions are:


.. doxygenfunction:: PeleC::init_eb


.. doxygenfunction:: PeleC::initialize_eb_structs



AMReX Building Blocks
---------------------

.. AMReX Building Blocks

.. include:: amrex_geometry.rst



Basic work iterator for for EB geometry
---------------------------------------

First fillpatch


.. code-block:: c

    {
      FillPatchIterator fpi(*this, coeff_cc, nGrowTr, time, State_Type, 0, NUM_STATE);
      MultiFab& S = fpi.get_mf(); 
    
      cons_to_prim(S,Q,Qaux);
    
      if (level > 0) 
      {   
        const BoxArray& crse_grids = getLevel(level-1).boxArray();
        const DistributionMapping& dmc = getLevel(level-1).DistributionMap();
        MultiFab Sc(crse_grids,dmc,NUM_STATE,nGrowTr);
        FillPatch(getLevel(level-1),Sc,nGrowTr,time,State_Type,0,NUM_STATE);
    
        Qc.define(crse_grids,dmc,QVAR,nGrowTr);
        Qcaux.define(crse_grids,dmc,NQAUX>0?NQAUX:1,nGrowTr);
        cons_to_prim(Sc,Qc,Qcaux);
      }
    }




Then iterate over `MultiFab`s:


.. code-block:: c

    EBFArrayBox& feb = static_cast<EBFArrayBox&>(Q[mfi]);
    const auto& flag_fab = feb.getEBCellFlagFab();
    FabType typ = flag_fab.getType(cbox);
    if (typ == FabType::covered)
    {
    }
    else if (typ == FabType::singlevalued)
    {
      pc_compute_tangential_vel_derivs_eb(cbox.loVect(), cbox.hiVect(),
                                          BL_TO_FORTRAN_3D(Q[mfi]),
                                          BL_TO_FORTRAN_3D(tander_ec[d]),
                                          BL_TO_FORTRAN_ANYD(flag_fab),
                                          geom.CellSize(),&d);
    }
    else if (typ == FabType::multivalued)
    {
      amrex::Abort("multi-valued eb tangential derivatives to be implemented");
    }
    else
    {
      pc_compute_tangential_vel_derivs(cbox.loVect(), cbox.hiVect(),
                                       BL_TO_FORTRAN_3D(Q[mfi]),
                                       BL_TO_FORTRAN_3D(tander_ec[d]),
                                       geom.CellSize(),&d);
    }



