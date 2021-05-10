// clang-format off
#ifdef PELEC_USE_EB

#include <iamr_redistribution.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 1)

void
Redistribution::MakeITracker ( Box const& bx,
                               Array4<Real const> const& apx,
                               Array4<Real const> const& apy,
                               Array4<Real const> const& vfrac,
                               Array4<int> const& itracker,
                               Geometry const& lev_geom,
                               std::string redist_type)
{
  amrex::Abort("MakeITracker not supported in 1D.");
}
#endif
#endif
// clang-format on
