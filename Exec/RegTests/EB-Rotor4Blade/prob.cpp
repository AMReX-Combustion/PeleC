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
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p0", PeleC::h_prob_parm_device->p0);
  pp.query("T0", PeleC::h_prob_parm_device->T0);
  pp.query("d_inner", PeleC::h_prob_parm_device->d_inner);
  pp.query("d_outer", PeleC::h_prob_parm_device->d_outer);
  pp.query("viscosity", PeleC::h_prob_parm_device->mu);
  pp.query("thermal_conductivity", PeleC::h_prob_parm_device->kappa);

  // check for rotational frame
  amrex::ParmParse pp_pelec("pelec");
  pp_pelec.query("do_rf", PeleC::h_prob_parm_device->do_rf);
  pp_pelec.query("rf_omega", PeleC::h_prob_parm_device->rf_omega);
  pp_pelec.query("rf_axis", PeleC::h_prob_parm_device->rf_axis);
  pp_pelec.query("rf_axis_x", PeleC::h_prob_parm_device->rf_axis_x);
  pp_pelec.query("rf_axis_y", PeleC::h_prob_parm_device->rf_axis_y);
  pp_pelec.query("rf_axis_z", PeleC::h_prob_parm_device->rf_axis_z);

  auto& trans_parm = PeleC::trans_parms.host_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::h_prob_parm_device->mu;
  trans_parm.const_conductivity = PeleC::h_prob_parm_device->kappa;
  PeleC::trans_parms.sync_to_device();
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

void
EBConcCylinders::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  auto plo = geom.ProbLoArray();
  auto phi = geom.ProbHiArray();
  amrex::ParmParse ppprob("prob");
  ppprob.query("d_inner", PeleC::h_prob_parm_device->d_inner);
  ppprob.query("d_outer", PeleC::h_prob_parm_device->d_outer);

  ProbParmDevice const* pp = PeleC::h_prob_parm_device;

  amrex::Real r_in = 0.5 * pp->d_inner;
  amrex::Real bladelen = 1.5 * r_in;
  amrex::Real bladethick = 0.4 * r_in;
  amrex::Real zlen = phi[2] - plo[2];

  // infinite cylinders
  amrex::EB2::CylinderIF inner(r_in, 2, {AMREX_D_DECL(0.0, 0, 0)}, false);

  amrex::RealArray lo, hi;

  lo[0] = 0.5 * r_in;
  hi[0] = lo[0] + bladelen;
  lo[1] = -0.5 * bladethick;
  hi[1] = 0.5 * bladethick;
  lo[2] = plo[2] - 2.0 * zlen;
  hi[2] = phi[2] + 2.0 * zlen;
  amrex::EB2::BoxIF bf1(lo, hi, false);

  hi[0] = -0.5 * r_in;
  lo[0] = hi[0] - bladelen;
  lo[1] = -0.5 * bladethick;
  hi[1] = 0.5 * bladethick;
  lo[2] = plo[2] - 2.0 * zlen;
  hi[2] = phi[2] + 2.0 * zlen;
  amrex::EB2::BoxIF bf2(lo, hi, false);

  lo[1] = 0.5 * r_in;
  hi[1] = lo[1] + bladelen;
  lo[0] = -0.5 * bladethick;
  hi[0] = 0.5 * bladethick;
  lo[2] = plo[2] - 2.0 * zlen;
  hi[2] = phi[2] + 2.0 * zlen;
  amrex::EB2::BoxIF bf3(lo, hi, false);

  hi[1] = -0.5 * r_in;
  lo[1] = hi[1] - bladelen;
  lo[0] = -0.5 * bladethick;
  hi[0] = 0.5 * bladethick;
  lo[2] = plo[2] - 2.0 * zlen;
  hi[2] = phi[2] + 2.0 * zlen;
  amrex::EB2::BoxIF bf4(lo, hi, false);

  amrex::EB2::CylinderIF outer(
    0.5 * pp->d_outer, 2, {AMREX_D_DECL(0.0, 0, 0)}, true);

  auto conc_cyl = amrex::EB2::makeUnion(inner, bf1, bf2, bf3, bf4, outer);
  auto gshop = amrex::EB2::makeShop(conc_cyl);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}
