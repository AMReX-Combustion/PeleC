#include "prob.H"

void
EBLinePistonCylinder::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::EB2::SplineIF Piston;
  std::vector<amrex::RealVect> lnpts;
  amrex::RealVect p;

  const amrex::Real scaleFact = 0.25;

  p = amrex::RealVect(
    AMREX_D_DECL(49.0 * 0.1 * scaleFact, 7.8583 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  p = amrex::RealVect(
    AMREX_D_DECL(36.193 * 0.1 * scaleFact, 7.8583 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  Piston.addLineElement(lnpts);
  lnpts.clear();

  p = amrex::RealVect(
    AMREX_D_DECL(36.193 * 0.1 * scaleFact, 7.8583 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  p = amrex::RealVect(
    AMREX_D_DECL(24.035 * 0.1 * scaleFact, -7.8586 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  Piston.addLineElement(lnpts);
  lnpts.clear();

  p = amrex::RealVect(
    AMREX_D_DECL(24.035 * 0.1 * scaleFact, -7.8586 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  p = amrex::RealVect(
    AMREX_D_DECL(20.0 * 0.1 * scaleFact, -7.8586 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  Piston.addLineElement(lnpts);
  lnpts.clear();

  p = amrex::RealVect(
    AMREX_D_DECL(20.0 * 0.1 * scaleFact, -7.8586 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  p = amrex::RealVect(
    AMREX_D_DECL(1.9934 * 0.1 * scaleFact, 3.464 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  Piston.addLineElement(lnpts);
  lnpts.clear();

  p = amrex::RealVect(
    AMREX_D_DECL(1.9934 * 0.1 * scaleFact, 3.464 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  p = amrex::RealVect(
    AMREX_D_DECL(0.09061 * 0.1 * scaleFact, 3.464 * 0.1 * scaleFact, 0.0));
  lnpts.push_back(p);
  Piston.addLineElement(lnpts);

  amrex::EB2::CylinderIF cylinder(
    48.0 * 0.1 * scaleFact, 70.0 * 0.1 * scaleFact, 2,
    {0.0, 0.0, -10.0 * 0.1 * scaleFact}, true);

  auto revolvePiston = amrex::EB2::lathe(Piston);
  auto PistonCylinder = amrex::EB2::makeUnion(revolvePiston, cylinder);
  auto gshop = amrex::EB2::makeShop(PistonCylinder);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  amrex::ParmParse pp("prob");

  // Chamber conditions
  pp.query("P_mean", PeleC::h_prob_parm_device->P_mean);
  pp.query("T_mean", PeleC::h_prob_parm_device->T_mean);
  pp.query("Y_CH4_chamber", PeleC::h_prob_parm_device->Y_CH4_chamber);
  pp.query("Y_O2_chamber", PeleC::h_prob_parm_device->Y_O2_chamber);

  // Injection parameters
  pp.query("nholes", PeleC::h_prob_parm_device->nholes);
  pp.query("cone_angle", PeleC::h_prob_parm_device->cone_angle);
  pp.query("centx", PeleC::h_prob_parm_device->centx);
  pp.query("centy", PeleC::h_prob_parm_device->centy);
  pp.query("r_circ", PeleC::h_prob_parm_device->r_circ);
  pp.query("r_hole", PeleC::h_prob_parm_device->r_hole);
  pp.query("T_jet", PeleC::h_prob_parm_device->T_jet);
  pp.query("vel_jet", PeleC::h_prob_parm_device->vel_jet);
  pp.query("injection_start", PeleC::h_prob_parm_device->inj_start);
  pp.query("injection_duration", PeleC::h_prob_parm_device->inj_dur);
  pp.query("tau", PeleC::h_prob_parm_device->tau);
  pp.query("Z", PeleC::h_prob_parm_device->Z);

  amrex::ParmParse ppic("ic");
  ppic.query("hitIC", PeleC::h_prob_parm_device->hitIC);
  ppic.query("inres", PeleC::h_prob_parm_device->inres);
  ppic.query("binfmt", PeleC::h_prob_parm_device->binfmt);
  ppic.query("iname", PeleC::prob_parm_host->iname);
  ppic.query("uin_norm", PeleC::h_prob_parm_device->uin_norm);
  ppic.query("lscale", PeleC::h_prob_parm_device->lscale);
  ppic.query("offset", PeleC::h_prob_parm_device->offset);
  ppic.query("urms0", PeleC::h_prob_parm_device->urms0);
  amrex::Vector<amrex::Real> win_lo(
    AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::lowest());
  amrex::Vector<amrex::Real> win_hi(
    AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::max());
  ppic.queryarr("win_lo", win_lo, 0, AMREX_SPACEDIM);
  ppic.queryarr("win_hi", win_hi, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    PeleC::h_prob_parm_device->win_lo[i] = win_lo[i];
    PeleC::h_prob_parm_device->win_hi[i] = win_hi[i];
  }
  ppic.query("win_slope", PeleC::h_prob_parm_device->win_slope);

  if (not PeleC::h_prob_parm_device->hitIC) {
    amrex::Print() << "Skipping HIT IC input file reading and assuming restart."
                   << std::endl;
  } else {
    const size_t nx = PeleC::h_prob_parm_device->inres;
    const size_t ny = PeleC::h_prob_parm_device->inres;
    const size_t nz = PeleC::h_prob_parm_device->inres;
    amrex::Vector<double> data(nx * ny * nz * 6);
    if (PeleC::h_prob_parm_device->binfmt) {
      read_binary(PeleC::prob_parm_host->iname, nx, ny, nz, 6, data);
    } else {
      read_csv(PeleC::prob_parm_host->iname, nx, ny, nz, data);
    }

    // Extract position and velocities
    PeleC::prob_parm_host->h_xinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->xinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_uinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->uinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_vinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->vinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_winput.resize(nx * ny * nz);
    PeleC::prob_parm_host->winput.resize(nx * ny * nz);
    for (int i = 0; i < PeleC::prob_parm_host->h_xinput.size(); i++) {
      PeleC::prob_parm_host->h_xinput[i] =
        (data[0 + i * 6] + PeleC::h_prob_parm_device->offset) /
        PeleC::h_prob_parm_device->lscale;
      PeleC::prob_parm_host->h_uinput[i] = data[3 + i * 6] *
                                           PeleC::h_prob_parm_device->urms0 /
                                           PeleC::h_prob_parm_device->uin_norm;
      PeleC::prob_parm_host->h_vinput[i] = data[4 + i * 6] *
                                           PeleC::h_prob_parm_device->urms0 /
                                           PeleC::h_prob_parm_device->uin_norm;
      PeleC::prob_parm_host->h_winput[i] = data[5 + i * 6] *
                                           PeleC::h_prob_parm_device->urms0 /
                                           PeleC::h_prob_parm_device->uin_norm;
    }

    // Get the xarray table and the differences.
    PeleC::prob_parm_host->h_xarray.resize(nx);
    PeleC::prob_parm_host->xarray.resize(nx);
    for (int i = 0; i < PeleC::prob_parm_host->h_xarray.size(); i++) {
      PeleC::prob_parm_host->h_xarray[i] = PeleC::prob_parm_host->h_xinput[i];
    }
    PeleC::prob_parm_host->h_xdiff.resize(nx);
    PeleC::prob_parm_host->xdiff.resize(nx);
    std::adjacent_difference(
      PeleC::prob_parm_host->h_xarray.begin(),
      PeleC::prob_parm_host->h_xarray.end(),
      PeleC::prob_parm_host->h_xdiff.begin());
    PeleC::prob_parm_host->h_xdiff[0] = PeleC::prob_parm_host->h_xdiff[1];

    // Dimensions of the input box.
    PeleC::h_prob_parm_device->Linput =
      PeleC::prob_parm_host->h_xarray[nx - 1] +
      0.5 * PeleC::prob_parm_host->h_xdiff[nx - 1];
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_xinput.begin(),
    PeleC::prob_parm_host->h_xinput.end(),
    PeleC::prob_parm_host->xinput.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_uinput.begin(),
    PeleC::prob_parm_host->h_uinput.end(),
    PeleC::prob_parm_host->uinput.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_vinput.begin(),
    PeleC::prob_parm_host->h_vinput.end(),
    PeleC::prob_parm_host->vinput.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_winput.begin(),
    PeleC::prob_parm_host->h_winput.end(),
    PeleC::prob_parm_host->winput.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_xarray.begin(),
    PeleC::prob_parm_host->h_xarray.end(),
    PeleC::prob_parm_host->xarray.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_xdiff.begin(),
    PeleC::prob_parm_host->h_xdiff.end(), PeleC::prob_parm_host->xdiff.begin());

  // Get pointers to the data
  PeleC::h_prob_parm_device->d_xinput = PeleC::prob_parm_host->xinput.data();
  PeleC::h_prob_parm_device->d_uinput = PeleC::prob_parm_host->uinput.data();
  PeleC::h_prob_parm_device->d_vinput = PeleC::prob_parm_host->vinput.data();
  PeleC::h_prob_parm_device->d_winput = PeleC::prob_parm_host->winput.data();
  PeleC::h_prob_parm_device->d_xarray = PeleC::prob_parm_host->xarray.data();
  PeleC::h_prob_parm_device->d_xdiff = PeleC::prob_parm_host->xdiff.data();
}
}

void
pc_prob_close()
{
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
