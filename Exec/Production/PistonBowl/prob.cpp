#include "prob.H"

void
ReadPipeInflow(const std::string iname)
{
  std::ifstream infile(iname, std::ios::in | std::ios::binary);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }

  double tmp_double = 0.0; /*needs to be double */
  PeleC::h_prob_parm_device->inflowNtime = read_binary_int(infile);
  PeleC::h_prob_parm_device->nr = read_binary_int(infile);
  PeleC::h_prob_parm_device->nt = read_binary_int(infile);
  PeleC::h_prob_parm_device->inflowFreq = read_binary_double(infile);
  PeleC::h_prob_parm_device->timeFromInflow = read_binary_double(infile);

  // Time grid
  PeleC::prob_parm_host->h_timeInput.resize(
    PeleC::h_prob_parm_device->inflowNtime);
  PeleC::prob_parm_host->timeInput.resize(
    PeleC::h_prob_parm_device->inflowNtime);
  for (int i = 0; i < PeleC::prob_parm_host->h_timeInput.size(); i++) {
    PeleC::prob_parm_host->h_timeInput[i] =
      i * PeleC::h_prob_parm_device->inflowFreq;
  }
  PeleC::h_prob_parm_device->timeInflowMax = *std::max_element(
    PeleC::prob_parm_host->h_timeInput.begin(),
    PeleC::prob_parm_host->h_timeInput.end());

  // Radial grid
  PeleC::prob_parm_host->rInput.resize(PeleC::h_prob_parm_device->nr + 1, 0.0);
  for (int i = 0; i < PeleC::prob_parm_host->rInput.size(); i++) {
    PeleC::prob_parm_host->rInput[i] = read_binary_double(infile);
  }
  PeleC::prob_parm_host->h_rM.resize(PeleC::h_prob_parm_device->nr, 0.0);
  PeleC::prob_parm_host->rM.resize(PeleC::h_prob_parm_device->nr, 0.0);
  for (int i = 0; i < PeleC::prob_parm_host->h_rM.size(); i++) {
    PeleC::prob_parm_host->h_rM[i] =
      PeleC::prob_parm_host->rInput[i + 1] * 0.5 +
      PeleC::prob_parm_host->rInput[i] * 0.5;
  }

  // Azimuthal grid
  PeleC::prob_parm_host->thetaInput.resize(
    PeleC::h_prob_parm_device->nt + 1, 0.0);
  for (int i = 0; i < PeleC::prob_parm_host->thetaInput.size(); i++) {
    PeleC::prob_parm_host->thetaInput[i] = read_binary_double(infile);
  }
  PeleC::prob_parm_host->h_thetaM.resize(PeleC::h_prob_parm_device->nt, 0.0);
  PeleC::prob_parm_host->thetaM.resize(PeleC::h_prob_parm_device->nt, 0.0);
  for (int i = 0; i < PeleC::prob_parm_host->h_thetaM.size(); i++) {
    PeleC::prob_parm_host->h_thetaM[i] =
      0.5 * PeleC::prob_parm_host->thetaInput[i + 1] +
      0.5 * PeleC::prob_parm_host->thetaInput[i];
  }
  PeleC::h_prob_parm_device->thetaMax = *std::max_element(
    PeleC::prob_parm_host->h_thetaM.begin(),
    PeleC::prob_parm_host->h_thetaM.end());

  // Read in velocities
  PeleC::prob_parm_host->h_Uz.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  PeleC::prob_parm_host->Uz.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  for (int k = 0; k < PeleC::h_prob_parm_device->inflowNtime; k++) {
    for (int j = 0; j < PeleC::h_prob_parm_device->nr; j++) {
      for (int i = 0; i < PeleC::h_prob_parm_device->nt; i++) {
        PeleC::prob_parm_host->h_Uz
          [(k * PeleC::h_prob_parm_device->nr + j) *
             PeleC::h_prob_parm_device->nt +
           i] = read_binary_double(infile);
      }
    }
  }

  PeleC::prob_parm_host->h_Ur.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  PeleC::prob_parm_host->Ur.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  for (int k = 0; k < PeleC::h_prob_parm_device->inflowNtime; k++) {
    for (int j = 0; j < PeleC::h_prob_parm_device->nr; j++) {
      for (int i = 0; i < PeleC::h_prob_parm_device->nt; i++) {
        PeleC::prob_parm_host->h_Ur
          [(k * PeleC::h_prob_parm_device->nr + j) *
             PeleC::h_prob_parm_device->nt +
           i] = read_binary_double(infile);
      }
    }
  }

  PeleC::prob_parm_host->h_Ut.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  PeleC::prob_parm_host->Ut.resize(
    PeleC::h_prob_parm_device->inflowNtime * PeleC::h_prob_parm_device->nr *
      PeleC::h_prob_parm_device->nt,
    0.0);
  for (int k = 0; k < PeleC::h_prob_parm_device->inflowNtime; k++) {
    for (int j = 0; j < PeleC::h_prob_parm_device->nr; j++) {
      for (int i = 0; i < PeleC::h_prob_parm_device->nt; i++) {
        PeleC::prob_parm_host->h_Ut
          [(k * PeleC::h_prob_parm_device->nr + j) *
             PeleC::h_prob_parm_device->nt +
           i] = read_binary_double(infile);
      }
    }
  }

  infile.close();
}

void
EBLinePistonCylinder(
  const amrex::Geometry& geom, const int required_level, const int max_level)
{
  amrex::EB2::SplineIF Piston;
  std::vector<amrex::RealVect> lnpts;
  amrex::RealVect p;
  int max_coarsening_level = max_level; // Because there are no mg solvers here

  amrex::Real scaleFact;
  scaleFact = 0.25;

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
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("Pres_domain", PeleC::h_prob_parm_device->Pres_domain);
  pp.query("Temp_domain", PeleC::h_prob_parm_device->Temp_domain);
  pp.query("Yfuel_domain", PeleC::h_prob_parm_device->Yfuel_domain);
  pp.query("Yox_domain", PeleC::h_prob_parm_device->Yox_domain);
  pp.query("YN2_domain", PeleC::h_prob_parm_device->YN2_domain);
  pp.query("dens_jet", PeleC::h_prob_parm_device->dens_jet);
  pp.query("vel_jet", PeleC::h_prob_parm_device->vel_jet);
  pp.query("Yox_jet", PeleC::h_prob_parm_device->Yox_jet);
  pp.query("Yfuel_jet", PeleC::h_prob_parm_device->Yfuel_jet);
  pp.query("YN2_jet", PeleC::h_prob_parm_device->YN2_jet);
  pp.query("centx", PeleC::h_prob_parm_device->centx);
  pp.query("centz", PeleC::h_prob_parm_device->centz);
  pp.query("r_circ", PeleC::h_prob_parm_device->r_circ);
  pp.query("r_hole", PeleC::h_prob_parm_device->r_hole);
  pp.query("nholes", PeleC::h_prob_parm_device->nholes);
  pp.query("cone_angle", PeleC::h_prob_parm_device->cone_angle);
  pp.query("inj_time", PeleC::h_prob_parm_device->inj_time);

  amrex::ParmParse ppic("ic");
  ppic.query("hitIC", PeleC::h_prob_parm_device->hitIC);
  ppic.query("inres", PeleC::h_prob_parm_device->inres);
  ppic.query("binfmt", PeleC::h_prob_parm_device->binfmt);
  ppic.query("iname", PeleC::prob_parm_host->iname);
  ppic.query("uin_norm", PeleC::h_prob_parm_device->uin_norm);
  ppic.query("lscale", PeleC::h_prob_parm_device->lscale);
  ppic.query("offset", PeleC::h_prob_parm_device->offset);
  ppic.query("urms0", PeleC::h_prob_parm_device->urms0);

  amrex::ParmParse ppin("inflow");
  ppin.query("inflowProfileFile", PeleC::prob_parm_host->inflowProfileFile);
  ppin.query("turb_inflow_type", PeleC::h_prob_parm_device->turb_inflow_type);

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

  // Read inflow
  if (PeleC::h_prob_parm_device->turb_inflow_type == 1) {
    ReadPipeInflow(PeleC::prob_parm_host->inflowProfileFile);
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_timeInput.begin(),
    PeleC::prob_parm_host->h_timeInput.end(),
    PeleC::prob_parm_host->timeInput.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_rM.begin(),
    PeleC::prob_parm_host->h_rM.end(), PeleC::prob_parm_host->rM.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_thetaM.begin(),
    PeleC::prob_parm_host->h_thetaM.end(),
    PeleC::prob_parm_host->thetaM.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_Uz.begin(),
    PeleC::prob_parm_host->h_Uz.end(), PeleC::prob_parm_host->Uz.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_Ur.begin(),
    PeleC::prob_parm_host->h_Ur.end(), PeleC::prob_parm_host->Ur.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_Ut.begin(),
    PeleC::prob_parm_host->h_Ut.end(), PeleC::prob_parm_host->Ut.begin());
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
  PeleC::h_prob_parm_device->d_timeInput =
    PeleC::prob_parm_host->timeInput.data();
  PeleC::h_prob_parm_device->d_rM = PeleC::prob_parm_host->rM.data();
  PeleC::h_prob_parm_device->d_thetaM = PeleC::prob_parm_host->thetaM.data();
  PeleC::h_prob_parm_device->d_Uz = PeleC::prob_parm_host->Uz.data();
  PeleC::h_prob_parm_device->d_Ur = PeleC::prob_parm_host->Ur.data();
  PeleC::h_prob_parm_device->d_Ut = PeleC::prob_parm_host->Ut.data();
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
