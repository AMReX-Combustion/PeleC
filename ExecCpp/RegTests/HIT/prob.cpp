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
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("iname", PeleC::prob_parm_host->iname);
    pp.query("binfmt", PeleC::prob_parm_device->binfmt);
    pp.query("restart", PeleC::prob_parm_device->restart);
    pp.query("lambda0", PeleC::prob_parm_device->lambda0);
    pp.query("reynolds_lambda0", PeleC::prob_parm_device->reynolds_lambda0);
    pp.query("mach_t0", PeleC::prob_parm_device->mach_t0);
    pp.query("prandtl", PeleC::prob_parm_device->prandtl);
    pp.query("inres", PeleC::prob_parm_device->inres);
    pp.query("uin_norm", PeleC::prob_parm_device->uin_norm);
  }

  {
    amrex::ParmParse pp("forcing");
    pp.query("u0", PeleC::prob_parm_device->forcing_u0);
    pp.query("v0", PeleC::prob_parm_device->forcing_v0);
    pp.query("w0", PeleC::prob_parm_device->forcing_w0);
    pp.query("forcing", PeleC::prob_parm_device->forcing_force);
  }

  // Define the length scale
  PeleC::prob_parm_device->L_x = probhi[0] - problo[0];
  PeleC::prob_parm_device->L_y = probhi[1] - problo[1];
  PeleC::prob_parm_device->L_z = probhi[2] - problo[2];

  // Wavelength associated to Taylor length scale
  PeleC::prob_parm_device->k0 = 2.0 / PeleC::prob_parm_device->lambda0;

  // Initial density, velocity, and material properties
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(
    PeleC::prob_parm_device->p0, massfrac, PeleC::prob_parm_device->T0,
    PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->eint0);
  EOS::RTY2Cs(
    PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->T0, massfrac, cs);
  EOS::TY2Cp(PeleC::prob_parm_device->T0, massfrac, cp);

  PeleC::prob_parm_device->urms0 =
    PeleC::prob_parm_device->mach_t0 * cs / sqrt(3.0);
  PeleC::prob_parm_device->tau =
    PeleC::prob_parm_device->lambda0 / PeleC::prob_parm_device->urms0;

  TransParm trans_parm;

  // Default
  trans_parm.const_viscosity = 0.0;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_conductivity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  // User-specified
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", trans_parm.const_viscosity);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_conductivity", trans_parm.const_conductivity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);
  }

  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = PeleC::prob_parm_device->rho0 *
                               PeleC::prob_parm_device->urms0 *
                               PeleC::prob_parm_device->lambda0 /
                               PeleC::prob_parm_device->reynolds_lambda0;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::prob_parm_device->prandtl;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "lambda0, k0, rho0, urms0, tau, p0, T0, gamma, mu, k, c_s0, Reynolds, "
       "Mach, Prandtl, u0, v0, w0, forcing"
    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << PeleC::prob_parm_device->lambda0 << "," << PeleC::prob_parm_device->k0
    << "," << PeleC::prob_parm_device->rho0 << ","
    << PeleC::prob_parm_device->urms0 << "," << PeleC::prob_parm_device->tau
    << "," << PeleC::prob_parm_device->p0 << "," << PeleC::prob_parm_device->T0
    << "," << EOS::gamma << "," << trans_parm.const_viscosity << ","
    << trans_parm.const_conductivity << "," << cs << ","
    << PeleC::prob_parm_device->reynolds_lambda0 << ","
    << PeleC::prob_parm_device->mach_t0 << ","
    << PeleC::prob_parm_device->prandtl << ","
    << PeleC::prob_parm_device->forcing_u0 << ","
    << PeleC::prob_parm_device->forcing_v0 << ","
    << PeleC::prob_parm_device->forcing_w0 << ","
    << PeleC::prob_parm_device->forcing_force << std::endl;
  ofs.close();

  // Load velocity fields from file. Assume data set ordered in Fortran
  // format and reshape the data accordingly. One thing to keep in mind
  // is that this contains the entire input data. We will interpolate
  // this data later to just match our box. Another assumption is that
  // the input data is a periodic cube. If the input cube is smaller
  // than our domain size, the cube will be repeated throughout the
  // domain (hence the mod operations in the interpolation).
  if (PeleC::prob_parm_device->restart) {
    amrex::Print() << "Skipping input file reading and assuming restart."
                   << std::endl;
  } else {
#ifdef AMREX_USE_FLOAT
    amrex::Abort("HIT cannot run in single precision at the moment.");
#else
    const size_t nx = PeleC::prob_parm_device->inres;
    const size_t ny = PeleC::prob_parm_device->inres;
    const size_t nz = PeleC::prob_parm_device->inres;
    amrex::Vector<amrex::Real> data(
      nx * ny * nz * 6); /* this needs to be double */
    if (PeleC::prob_parm_device->binfmt) {
      read_binary(PeleC::prob_parm_host->iname, nx, ny, nz, 6, data);
    } else {
      read_csv(PeleC::prob_parm_host->iname, nx, ny, nz, data);
    }

    // Extract position and velocities
    PeleC::prob_parm_host->h_xinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_uinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_vinput.resize(nx * ny * nz);
    PeleC::prob_parm_host->h_winput.resize(nx * ny * nz);
    for (long i = 0; i < PeleC::prob_parm_host->h_xinput.size(); i++) {
      PeleC::prob_parm_host->h_xinput[i] = data[0 + i * 6];
      PeleC::prob_parm_host->h_uinput[i] = data[3 + i * 6] *
                                           PeleC::prob_parm_device->urms0 /
                                           PeleC::prob_parm_device->uin_norm;
      PeleC::prob_parm_host->h_vinput[i] = data[4 + i * 6] *
                                           PeleC::prob_parm_device->urms0 /
                                           PeleC::prob_parm_device->uin_norm;
      PeleC::prob_parm_host->h_winput[i] = data[5 + i * 6] *
                                           PeleC::prob_parm_device->urms0 /
                                           PeleC::prob_parm_device->uin_norm;
    }

    // Get the xarray table and the differences.
    PeleC::prob_parm_host->h_xarray.resize(nx);
    for (long i = 0; i < PeleC::prob_parm_host->h_xarray.size(); i++) {
      PeleC::prob_parm_host->h_xarray[i] = PeleC::prob_parm_host->h_xinput[i];
    }
    PeleC::prob_parm_host->h_xdiff.resize(nx);
    std::adjacent_difference(
      PeleC::prob_parm_host->h_xarray.begin(),
      PeleC::prob_parm_host->h_xarray.end(),
      PeleC::prob_parm_host->h_xdiff.begin());
    PeleC::prob_parm_host->h_xdiff[0] = PeleC::prob_parm_host->h_xdiff[1];

    // Make sure the search array is increasing
    if (not std::is_sorted(
          PeleC::prob_parm_host->h_xarray.begin(),
          PeleC::prob_parm_host->h_xarray.end())) {
      amrex::Abort("Error: non ascending x-coordinate array.");
    }

    // Get pointer to the data
    PeleC::prob_parm_host->xinput.resize(
      PeleC::prob_parm_host->h_xinput.size());
    PeleC::prob_parm_host->uinput.resize(
      PeleC::prob_parm_host->h_uinput.size());
    PeleC::prob_parm_host->vinput.resize(
      PeleC::prob_parm_host->h_vinput.size());
    PeleC::prob_parm_host->winput.resize(
      PeleC::prob_parm_host->h_winput.size());
    PeleC::prob_parm_host->xarray.resize(
      PeleC::prob_parm_host->h_xarray.size());
    PeleC::prob_parm_host->xdiff.resize(PeleC::prob_parm_host->h_xdiff.size());
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
      PeleC::prob_parm_host->h_xdiff.end(),
      PeleC::prob_parm_host->xdiff.begin());

    PeleC::prob_parm_device->d_xinput = PeleC::prob_parm_host->xinput.data();
    PeleC::prob_parm_device->d_uinput = PeleC::prob_parm_host->uinput.data();
    PeleC::prob_parm_device->d_vinput = PeleC::prob_parm_host->vinput.data();
    PeleC::prob_parm_device->d_winput = PeleC::prob_parm_host->winput.data();
    PeleC::prob_parm_device->d_xarray = PeleC::prob_parm_host->xarray.data();
    PeleC::prob_parm_device->d_xdiff = PeleC::prob_parm_host->xdiff.data();

    // Dimensions of the input box.
    PeleC::prob_parm_device->Linput =
      PeleC::prob_parm_host->h_xarray[nx - 1] +
      0.5 * PeleC::prob_parm_host->h_xdiff[nx - 1];
#endif
  }
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
