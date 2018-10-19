
AMReX functions useful for Pele development
-------------------------------------------

Pele is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. Below is a quick reference to many of the AMReX tools used to build up Pele. The full AMReX documentation can be found `here <https://amrex-codes.github.io/AMReXUsersGuide.pdf>`_. 


.. 
.. 
.. Solution environment
.. --------------------
.. 
.. .. doxygenclass:: amrex::BoxArray
..    :project: PeleC
.. 
.. 
.. .. doxygenclass:: amrex::StateData
.. 
.. .. doxygenfunction:: amrex::StateData::copyOld
.. 
.. .. doxygenfunction:: amrex::StateData::setTimeLevel
.. 
.. .. doxygenfunction:: amrex::StateData::hasOldData
.. 
.. .. doxygenfunction:: amrex::StateData::swapTimeLevels
.. 
.. .. .. Need desc_list and derive_lst type...
.. 
.. Data structures
.. ---------------
.. 
.. .. doxygenclass:: amrex::Array
.. 
.. .. doxygenclass:: amrex::FArrayBox
.. 
.. .. doxygenclass:: BaseFab
.. .. 
.. .. .. doxygenfunction:: BaseFab< T >::dataPtr
.. .. 
.. .. .. doxygenfunction:: BaseFab::loVect
.. .. 
.. .. .. doxygenfunction:: BaseFab::hiVect
.. .. 
.. .. .. doxygenfunction:: BaseFab::nComp
.. .. 
.. .. .. doxygenfunction:: BaseFab::setVal(Real)
.. .. 
.. .. .. doxygenfunction:: amrex::FArrayBox::resize
.. .. 
.. .. .. doxygenfunction:: amrex::FArrayBox::plus
.. .. 
.. .. .. doxygenclass:: TagBox
.. .. 
.. .. .. doxygenclass:: Box
.. .. 
.. .. .. doxygenfunction:: amrex::Box::surroundingNodes(const Box&, int)
.. .. 
.. .. .. doxygenfunction:: amrex::grow(const Box&, int)
.. .. 
.. .. .. doxygenclass:: iMultiFab
.. .. 
.. .. .. doxygenclass:: amrex::MultiFab<FArrayBox>
.. .. 
.. .. .. doxygenfunction:: amrex::MultiFab::SumBoundary(int, int, const Periodicity&)
.. .. 
.. .. .. doxygenfunction:: amrex::MultiFab::setVal(value_type)
.. .. 
.. .. .. doxygenfunction:: Real amrex::MultiFab::sum
.. 
.. .. doxygenclass:: amrex::RealBox
.. 
.. 
.. .. Problem Description
.. .. -------------------
.. .. 
.. .. .. doxygenclass:: Geometry
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: Geometry::isPeriodic
.. ..    :project: PeleC
.. .. 
.. .. .. .. doxygenfunction:: Geometry::isCartesian
.. ..    :project: PeleC
.. .. 
.. .. PeleC Implementation 
.. .. --------------------
.. .. 
.. .. .. doxygenclass:: PeleC
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenclass:: AmrLevel
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenclass:: StateDescriptor
.. ..    :project: PeleC
.. .. 
.. .. .. .. doxygenfunction:: ctoprim
.. .. .. .. doxygenfunction:: srctoprim
.. .. .. .. doxygenfunction:: cons_to_prim
.. .. .. .. doxygenfunction:: check_for_nan
.. .. 
.. .. Tools for iterating over solution
.. .. ---------------------------------
.. .. 
.. .. .. doxygenclass:: MFIter
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenclass:: FillPatchIterator
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: FillPatch(AmrLevel&, MultiFab&, int, Real, int, int, int, int)
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: MultiFab::Copy
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: MultiFab::Add
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: MultiFab::define(const BoxArray&, int, int, FabAlloc, const IntVect&, ParallelDescriptor::Color)
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: MultiFab::clear
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: MultiFab::FillBoundary(bool)
.. ..    :project: PeleC
.. .. 
.. .. .. doxygenfunction:: amrex::MultiFab::Saxpy
.. ..    :project: PeleC
.. .. .. 
.. .. .. Multilevel tools
.. .. .. ----------------
.. .. .. .. doxygenfunction:: amrex::average_down(MultiFab&, MultiFab&, const Geometry&, const Geometry&, int, int, const int)
.. .. ..    :project: PeleC
.. .. .. 
.. .. .. .. doxygenclass:: FluxRegister
.. .. ..    :project: PeleC
.. .. .. 
.. .. .. .. doxygenfunction:: FluxRegister::FineAdd(const MultiFab&, int, int, int, int, Real)
.. .. ..    :project: PeleC
.. .. .. 
.. .. .. .. doxygenfunction:: FluxRegister::CrseInit(const MultiFab&, int, int, int, int, Real, FrOp)
.. .. ..    :project: PeleC
.. .. .. 
.. .. Boundary Condition Tools
.. .. ------------------------
.. .. .. doxygenclass:: BCRec
.. .. 
.. .. .. .. doxygenclass:: DiffusionBndry
.. .. 
.. .. .. doxygenfunction:: DiffusionBndry::setBndryDataGivenS
.. 
.. .. doxygenclass:: amrex::InterpBndryData
.. .. 
.. .. Linear Solvers
.. .. --------------
.. .. .. doxygenclass:: ABecLaplacian
.. .. 
.. .. .. doxygenfunction:: ABecLaplacian::applyBC





