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
  }
};

template <typename ProbDeriveStruct, typename... Args>
void
add_problem_derives(Args&&... args)
{
  ProbDeriveStruct::add(std::forward<Args>(args)...);
}
#endif
