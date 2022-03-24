#include "prob.H"

void
pc_prob_close()
{
}

void
parse_params(ProbParmDevice* prob_parm_device)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("d_cyl", prob_parm_device->d_cyl);
  pp.query("l_nozzle", prob_parm_device->l_nozzle);
  pp.query("d_nozzle", prob_parm_device->d_nozzle);
  pp.query("l_combustor", prob_parm_device->l_combustor);
  pp.query("l_exit", prob_parm_device->l_exit);
  pp.query("d_swirler", prob_parm_device->d_swirler);
  pp.query("u_swirl_ax", prob_parm_device->u_swirl_ax);
  pp.query("T_swirl", prob_parm_device->T_swirl);
  pp.query("p0", prob_parm_device->p0);
  pp.query("T0", prob_parm_device->T0);
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
  parse_params(PeleC::h_prob_parm_device);
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
EBConvergingNozzle::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  parse_params(PeleC::h_prob_parm_device);
  ProbParmDevice const* pp = PeleC::h_prob_parm_device;

  amrex::EB2::CylinderIF main(
    0.5 * pp->d_cyl, 2.0 * pp->l_combustor, 0,
    {AMREX_D_DECL(0.5 * pp->l_combustor, 0, 0)}, true);

  amrex::Real slope_nozzle =
    (0.5 * pp->d_cyl - 0.5 * pp->d_nozzle) / pp->l_nozzle;
  amrex::Real norm = -1.0 / slope_nozzle;
  amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
  amrex::EB2::PlaneIF nozzle_plane(
    {AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(1 / nmag, slope_nozzle / nmag, 0.0)},
    true);
  auto nozzle = amrex::EB2::translate(
    amrex::EB2::rotate(amrex::EB2::lathe(nozzle_plane), 90 * M_PI / 180, 1),
    {AMREX_D_DECL(pp->l_combustor + 0.5 * pp->d_nozzle / slope_nozzle, 0, 0)});

  amrex::EB2::CylinderIF exit(
    0.5 * pp->d_nozzle, 4 * pp->l_exit, 0,
    {AMREX_D_DECL(pp->l_combustor, 0, 0)}, true);
  auto nozzle_exit = amrex::EB2::makeIntersection(nozzle, exit);

  auto polys = amrex::EB2::makeUnion(main, nozzle_exit);
  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}
