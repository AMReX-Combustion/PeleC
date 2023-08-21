#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <AMReX_ErrorList.H>
#include "prob_parm.H"

struct TaggingParm
{
  amrex::Real denerr = 1.0e10;
  amrex::Real dengrad = 1.0e10;
  amrex::Real denratio = 1.0e10;
  amrex::Real presserr = 1.0e10;
  amrex::Real pressgrad = 1.0e10;
  amrex::Real velerr = 1.0e10;
  amrex::Real velgrad = 1.0e10;
  amrex::Real vorterr = 1.0e10;
  amrex::Real temperr = 1.0e10;
  amrex::Real lotemperr = -1.0e10;
  amrex::Real tempgrad = 1.0e10;
  amrex::Real ftracerr = 1.0e10;
  amrex::Real ftracgrad = 1.0e10;
  amrex::Real vfracerr = 1.0e10;
  amrex::Real detag_eb_factor = 3.0;

  std::string eb_refine_type = "static";

  int max_denerr_lev = 10;
  int max_dengrad_lev = 10;
  int max_denratio_lev = 10;
  int max_presserr_lev = 10;
  int max_pressgrad_lev = 10;
  int max_velerr_lev = 10;
  int max_velgrad_lev = 10;
  int max_vorterr_lev = 10;
  int max_temperr_lev = 10;
  int max_lotemperr_lev = 10;
  int max_tempgrad_lev = 10;
  int max_ftracerr_lev = 10;
  int max_ftracgrad_lev = 10;
  int max_vfracerr_lev = 10;
  int max_eb_refine_lev = 10;
  int min_eb_refine_lev = -1;
  int adapt_eb_refined_lev = -1;

  amrex::Vector<amrex::AMRErrorTag> err_tags;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if (field(iv) >= fielderr) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_loerror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of low field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if (field(iv) <= fielderr) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_graderror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldgrad,
  char tagval) noexcept
{
  // Tag on regions of high field gradient
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  AMREX_D_TERM(const auto ivp0 = iv + amrex::IntVect::TheDimensionVector(0);
               const auto ivm0 = iv - amrex::IntVect::TheDimensionVector(0);
               const amrex::Real axp = std::abs(field(ivp0) - field(iv));
               const amrex::Real axm = std::abs(field(iv) - field(ivm0));
               const amrex::Real ax = amrex::max<amrex::Real>(axp, axm);
               , const auto ivp1 = iv + amrex::IntVect::TheDimensionVector(1);
               const auto ivm1 = iv - amrex::IntVect::TheDimensionVector(1);
               const amrex::Real ayp = std::abs(field(ivp1) - field(iv));
               const amrex::Real aym = std::abs(field(iv) - field(ivm1));
               const amrex::Real ay = amrex::max<amrex::Real>(ayp, aym);
               , const auto ivp2 = iv + amrex::IntVect::TheDimensionVector(2);
               const auto ivm2 = iv - amrex::IntVect::TheDimensionVector(2);
               const amrex::Real azp = std::abs(field(ivp2) - field(iv));
               const amrex::Real azm = std::abs(field(iv) - field(ivm2));
               const amrex::Real az = amrex::max<amrex::Real>(azp, azm);)
#if AMREX_SPACEDIM > 1
  if (amrex::max<amrex::Real>(AMREX_D_DECL(ax, ay, az)) >= fieldgrad)
#elif AMREX_SPACEDIM == 1
  if (ax >= fieldgrad)
#endif
  {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_ratioerror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldratio,
  char tagval) noexcept
{
  // Tag on regions of high ratio of adjacent cell values in field
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  AMREX_D_TERM(const auto ivp0 = iv + amrex::IntVect::TheDimensionVector(0);
               const auto ivm0 = iv - amrex::IntVect::TheDimensionVector(0);
               amrex::Real axp1 = std::abs(field(ivp0) / field(iv));
               axp1 = amrex::max<amrex::Real>(axp1, 1.0 / axp1);
               amrex::Real axm1 = std::abs(field(ivm0) / field(iv));
               axm1 = amrex::max<amrex::Real>(axm1, 1.0 / axm1);
               const amrex::Real ax = amrex::max<amrex::Real>(axm1, axp1);
               , const auto ivp1 = iv + amrex::IntVect::TheDimensionVector(1);
               const auto ivm1 = iv - amrex::IntVect::TheDimensionVector(1);
               amrex::Real ayp1 = std::abs(field(ivp1) / field(iv));
               ayp1 = amrex::max<amrex::Real>(ayp1, 1.0 / ayp1);
               amrex::Real aym1 = std::abs(field(ivm1) / field(iv));
               aym1 = amrex::max<amrex::Real>(aym1, 1.0 / aym1);
               const amrex::Real ay = amrex::max<amrex::Real>(aym1, ayp1);
               , const auto ivp2 = iv + amrex::IntVect::TheDimensionVector(2);
               const auto ivm2 = iv - amrex::IntVect::TheDimensionVector(2);
               amrex::Real azp1 = std::abs(field(ivp2) / field(iv));
               azp1 = amrex::max<amrex::Real>(azp1, 1.0 / azp1);
               amrex::Real azm1 = std::abs(field(ivm2) / field(iv));
               azm1 = amrex::max<amrex::Real>(azm1, 1.0 / azm1);
               const amrex::Real az = amrex::max<amrex::Real>(azm1, azp1);)
#if AMREX_SPACEDIM > 1
  if (amrex::max<amrex::Real>(AMREX_D_DECL(ax, ay, az)) >= fieldratio)
#elif AMREX_SPACEDIM == 1
  if (ax >= fieldratio)
#endif
  {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_abserror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if (std::abs(field(iv)) >= fielderr) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error_bounds(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real lbnd,
  const amrex::Real ubnd,
  char tagval) noexcept
{
  // Tag on regions inside bounds
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if ((lbnd < field(iv)) && (field(iv) < ubnd)) {
    tag(iv) = tagval;
  }
}

struct EmptyProbTagStruct
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int /*i*/,
    const int /*j*/,
    const int /*k*/,
    amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
    amrex::Array4<char> const& /*tag*/,
    amrex::Array4<amrex::Real const> const& /*field*/,
    char /*tagval*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
    const amrex::Real /*time*/,
    const int /*level*/,
    ProbParmDevice const& /*d_prob_parm_device*/) noexcept
  {
  }
};

template <typename ProbTagStruct, typename... Args>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
set_problem_tags(Args&&... args)
{
  ProbTagStruct::set_problem_tags(std::forward<Args>(args)...);
}

#endif
