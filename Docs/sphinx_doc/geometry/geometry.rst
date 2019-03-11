.. highlight:: rst



Geometry treatment in PeleC
===========================


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







Initialization
--------------

Creating an EB geometry also requires knowledge of the finest level that will be used so that geometries that 'telescope', i.e., coarser volume fractions are consistent with applying the coarsening operator to the finer volumes, can be created. To that end there is a global geometry creation step, facilitated by the `initialize_EB2` function, as well as a step that happens when a new AMRLevel is created. The latter happens by a call to  `PeleC::initialize_eb2_structs`  through `PeleC::init_eb` called from the PeleC constructor. Following construction of the geometry, the geometric information is copied into the structures described in the previous section and the various interpolation stencils are populated. 

Cartesian grid, embedded boundary (EB) methods are methods where the geometric description is formed by cutting a Cartesian mesh with surface of the geometry.  AMReX's methods to handle EB geometry information, and PeleC's treatment of the EB aware update could use many possible sources for geometric description. The necessary information is, on a per-cell basis:

* Apertures for faces intersected by cut cells,
* cut cell volumes that 'telescope', that is, volumes at a coarser level are consistent with averaging the volumes from finer levels,
* connectivity indicating which neighbor cells are connected to a given cell, and
* coordinates of cell and face centroids. 

Additionally, the algorithms ultimately require surface normals, but these can be trivially recomputed from the aperture. 

In the following subsections, we will first describe using geometry creation tools based on EB infrastructure in AMReX with origins in a fork off Chombo's \cite{Chombo} infrastructure for setting up EB calculations and then how that information is used to populate the datastructures used to compute the fluxes discussed above. 

.. include:: geometry_creation.rst




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



