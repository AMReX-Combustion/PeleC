.. highlight:: rst

Geometry initialization
=======================

Creating an EB geometry also requires knowledge of the finest level that will be used so that geometries that 'telescope', 
i.e., coarser volume fractions are consistent with applying the coarsening operator to the finer volumes, can be created. 
To that end there is a global geometry creation step, facilitated by the `initialize_EB2` function, as well as a step that 
happens when a new AMRLevel is created. The latter happens by a call to  `PeleC::initialize_eb2_structs`  through `PeleC::init_eb` 
called from the PeleC constructor. Following construction of the geometry, the geometric information is 
copied into the structures described in the previous section and the various interpolation stencils are populated. 

Cartesian grid, embedded boundary (EB) methods are methods where the geometric description is formed by cutting a Cartesian 
mesh with surface of the geometry.  AMReX's methods to handle EB geometry information, and PeleC's treatment of the
EB aware update could use many possible sources for geometric description. The necessary information is, on a per-cell basis:

* Apertures for faces intersected by cut cells,
* cut cell volumes that 'telescope', that is, volumes at a coarser level are consistent with averaging the volumes from finer levels,
* connectivity indicating which neighbor cells are connected to a given cell, and
* coordinates of cell and face centroids. 

Additionally, the algorithms ultimately require surface normals, but these can be trivially recomputed from the aperture. 

In the following subsections, we will first describe using geometry creation tools based on EB infrastructure in AMReX with origins in a fork off Chombo's \cite{Chombo} infrastructure 
for setting up EB calculations and then how that information is used to populate the datastructures used to compute the fluxes discussed above. 

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



