// Problem-specific derives
#ifndef PROBLEMDERIVE_H
#define PROBLEMDERIVE_H

#include <AMReX_Derive.H>
#include <AMReX_FArrayBox.H>

#include "PeleC.H"
#include "IndexDefines.H"

struct EmptyProbDeriveStruct
{
  static void
  add(amrex::DeriveList& /*derive_lst*/, amrex::DescriptorList& /*desc_lst*/)
  {
    // Add derives as follows and define the derive function below:
    // derive_lst.add(
    //  "varname", amrex::IndexType::TheCellType(), 1, pc_varname,
    //  the_same_box);
    // derive_lst.addComponent("varname", desc_lst, State_Type, 0, NVAR);
  }

  static void pc_varname(
    const amrex::Box& /*bx*/,
    amrex::FArrayBox& /*derfab*/,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& /*datfab*/,
    const amrex::Geometry& /*geomdata*/,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
    // auto const dat = datfab.array();
    // auto arr = derfab.array();
    // amrex::ParallelFor(
    //  bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // do something with arr
    //      });
  }
};

template <typename ProbDeriveStruct, typename... Args>
void
add_problem_derives(Args&&... args)
{
  ProbDeriveStruct::add(std::forward<Args>(args)...);
}
#endif
