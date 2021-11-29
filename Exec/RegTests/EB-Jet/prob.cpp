#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_init", PeleC::h_prob_parm_device->p_init);
  pp.query("T_init", PeleC::h_prob_parm_device->T_init);
  pp.query("T_inner", PeleC::h_prob_parm_device->T_inner);
  pp.query("T_outer", PeleC::h_prob_parm_device->T_outer);
  pp.query("rad_temp", PeleC::h_prob_parm_device->rad_temp);
  pp.query("rad_jet", PeleC::h_prob_parm_device->rad_jet);
  pp.query("u_jet", PeleC::h_prob_parm_device->u_jet);
  pp.query("T_jet", PeleC::h_prob_parm_device->T_jet);

  // Initial values
  PeleC::h_prob_parm_device->massfrac[N2_ID] = 1.0;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->T_init, PeleC::h_prob_parm_device->rho_init,
    PeleC::h_prob_parm_device->e_init);
}
}

void
PeleC::problem_pre_advance()
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
    int local_i = mfi.LocalIndex();
    int Nvals = (!eb_in_domain) ? 0 : sv_eb_bcval[local_i].numPts();
    auto* ebg = (Nvals > 0 ? sv_eb_bndry_geom[local_i].data() : nullptr);
    auto* bcvals = sv_eb_bcval[local_i].dataPtr();
    const ProbParmDevice* pp = d_prob_parm_device;

    amrex::ParallelFor(Nvals, [=] AMREX_GPU_DEVICE(int L) {
      const amrex::IntVect iv = ebg[L].iv;
      const auto& prob_lo = geom.ProbLoArray();
      const auto& prob_hi = geom.ProbHiArray();
      const auto& dx = geom.CellSizeArray();
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> x({AMREX_D_DECL(
        prob_lo[0] + (iv[0] + 0.5) * dx[0], prob_lo[1] + (iv[1] + 0.5) * dx[1],
        prob_lo[2] + (iv[2] + 0.5) * dx[2])});

      amrex::Real ebnorm[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};
      const amrex::Real ebnorm_mag = std::sqrt(AMREX_D_TERM(
        ebnorm[0] * ebnorm[0], +ebnorm[1] * ebnorm[1], +ebnorm[2] * ebnorm[2]));
      for (amrex::Real& dir : ebnorm) {
        dir /= ebnorm_mag;
      }
      const amrex::Real theta = std::atan2(ebnorm[1], ebnorm[0]);
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xp({AMREX_D_DECL(
        std::cos(theta) * x[0] - std::sin(theta) * x[1],
        std::sin(theta) * x[0] + std::cos(theta) * x[1], x[2])});
      const amrex::Real xpc = 0.5 * (prob_hi[0] - prob_lo[0]) * std::cos(theta);
      const amrex::Real zc = 0.0;
      const amrex::Real rad =
        std::sqrt((xp[0] - xpc) * (xp[0] - xpc) + (xp[2] - zc) * (xp[2] - zc));

      bcvals[QTEMP * Nvals + L] = eb_boundary_T;
      bcvals[QU * Nvals + L] = 0.0;
      bcvals[QV * Nvals + L] = 0.0;
      bcvals[QW * Nvals + L] = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        bcvals[(QFS + n) * Nvals + L] = pp->massfrac[n];
      }

      if (rad < pp->rad_temp) {
        bcvals[QTEMP * Nvals + L] = pp->T_inner;
      } else {
        bcvals[QTEMP * Nvals + L] = pp->T_outer;
      }

      if (rad < pp->rad_jet) {
        bcvals[QTEMP * Nvals + L] = pp->T_jet;
        bcvals[QU * Nvals + L] = pp->u_jet * std::cos(theta);
        bcvals[QV * Nvals + L] = pp->u_jet * std::sin(theta);
        bcvals[QW * Nvals + L] = 0.0;
        // for (int n = 0; n < NUM_SPECIES; n++) {
        //   bcvals[(QFS + n) * Nvals + L] = 0.0;
        // }
        // bcvals[(QFS + O2_ID) * Nvals + L] = 1.0;
      }
    });
  }
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
