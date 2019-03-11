
AMReX functions useful for Pele development
-------------------------------------------

Pele is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. Below is a quick reference list with links to many of the AMReX tools used to build up PeleC. The full AMReX documentation can be found `here <https://amrex-codes.github.io/AMReXUsersGuide.pdf>`_. 


Solution environment
~~~~~~~~~~~~~~~~~~~~

* `amrex::BoxArray <https://amrex-codes.github.io/amrex/docs_html/Basics.html#boxarray>`_
* `amrex::StateData <https://amrex-codes.github.io/amrex/docs_html/AmrLevel.html?highlight=statedata#statedata>`_

  * amrex::StateData::copyOld
  * amrex::StateData::setTimeLevel
  * amrex::StateData::hasOldData
  * amrex::StateData::swapTimeLevels

Data structures
~~~~~~~~~~~~~~~

* `amrex::Array <https://amrex-codes.github.io/amrex/docs_html/Basics.html#vector-and-array>`_
* `amrex::FArrayBox <https://amrex-codes.github.io/amrex/docs_html/Basics.html#basefab-farraybox-and-iarraybox>`_

  * BaseFab< T >::dataPtr
  * BaseFab< T >::loVect
  * BaseFab< T >::hiVect
  * BaseFab< T >::nComp
  * BaseFab< T >::setVal(Real)
  * amrex::FArrayBox::resize
  * amrex::FArrayBox::plus

* `TagBox <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=tagbox#tagbox-and-cluster>`_
* `Box <https://amrex-codes.github.io/amrex/docs_html/Basics.html#box-intvect-and-indextype>`_

  * amrex::Box::surroundingNodes(const Box&, int)
  * amrex::grow(const Box&, int)

* `iMultiFab, amrex::MultiFab<FArraybox> <https://amrex-codes.github.io/amrex/docs_html/Basics.html#fabarray-multifab-and-imultifab>`_
   * amrex::MultiFab::SumBoundary(int, int, const Periodicity&)
   * amrex::MultiFab::setVal(value_type)
   * amrex::MultiFab::Copy
   * amrex::MultiFab::Add
   * amrex::MultiFab::define(const BoxArray&, int, int, FabAlloc, const IntVect&, ParallelDescriptor::Color)
   * amrex::MultiFab::sum
   * amrex::MultiFab::clear
   * amrex::MultiFab::FillBoundary
   * amrex::MultiFab::Saxpy
   * amrex::RealBox

PeleC Implementation 
~~~~~~~~~~~~~~~~~~~~

* `amrex::AmrLevel <https://amrex-codes.github.io/amrex/docs_html/AmrLevel.html#amrlevel-class>`_
* StateDescriptor
* `MFIter <https://amrex-codes.github.io/amrex/docs_html/Basics.html#mfiter-and-tiling>`_
* `FillPatch(AmrLevel&, MultiFab&, int, Real, int, int, int, int) <https://amrex-codes.github.io/amrex/docs_html/AsyncIter.html?highlight=fillpatch>`_



Multilevel tools
~~~~~~~~~~~~~~~~
* `amrex::average_down(MultiFab&, MultiFab&, const Geometry&, const Geometry&, int, int, const int) <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=averagedown>`_
* `FluxRegister <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=fluxregister#using-fluxregisters>`_
   * FluxRegister::FineAdd(const MultiFab&, int, int, int, int, Real)
   * FluxRegister::CrseInit(const MultiFab&, int, int, int, int, Real, FrOp)
 




