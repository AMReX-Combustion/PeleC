#include "prob.H"

namespace ProbParm {
std::string iname = "";
AMREX_GPU_DEVICE_MANAGED bool binfmt = false;
AMREX_GPU_DEVICE_MANAGED bool restart = false;
AMREX_GPU_DEVICE_MANAGED amrex::Real lambda0 = 0.5;
AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds_lambda0 = 100.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real mach_t0 = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 0.71;
AMREX_GPU_DEVICE_MANAGED int inres = 0;
AMREX_GPU_DEVICE_MANAGED amrex::Real uin_norm = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Linput = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real k0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real urms0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real tau = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real eint0 = 0.0;
amrex::Gpu::ManagedVector<amrex::Real>* v_xinput = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_uinput = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_vinput = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_winput = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_xarray = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_xdiff = nullptr;

AMREX_GPU_DEVICE_MANAGED amrex::Real* xinput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* uinput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* vinput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* winput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* xarray = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* xdiff = nullptr;

} // namespace ProbParm

void
pc_prob_close()
{
  delete ProbParm::v_xinput;
  delete ProbParm::v_uinput;
  delete ProbParm::v_vinput;
  delete ProbParm::v_winput;
  delete ProbParm::v_xarray;
  delete ProbParm::v_xdiff;

  ProbParm::v_xinput = nullptr;
  ProbParm::v_uinput = nullptr;
  ProbParm::v_vinput = nullptr;
  ProbParm::v_winput = nullptr;
  ProbParm::v_xarray = nullptr;
  ProbParm::v_xdiff = nullptr;
  ProbParm::xinput = nullptr;
  ProbParm::uinput = nullptr;
  ProbParm::vinput = nullptr;
  ProbParm::winput = nullptr;
  ProbParm::xarray = nullptr;
  ProbParm::xdiff = nullptr;
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
  amrex::ParmParse pp("prob");
  pp.query("iname", ProbParm::iname);
  pp.query("binfmt", ProbParm::binfmt);
  pp.query("restart", ProbParm::restart);
  pp.query("lambda0", ProbParm::lambda0);
  pp.query("reynolds_lambda0", ProbParm::reynolds_lambda0);
  pp.query("mach_t0", ProbParm::mach_t0);
  pp.query("prandtl", ProbParm::prandtl);
  pp.query("inres", ProbParm::inres);
  pp.query("uin_norm", ProbParm::uin_norm);

  amrex::ParmParse ppf("forcing");
  ppf.query("u0", forcing_params::u0);
  ppf.query("v0", forcing_params::v0);
  ppf.query("w0", forcing_params::w0);
  ppf.query("forcing", forcing_params::forcing);

  // Define the length scale
  ProbParm::L_x = probhi[0] - problo[0];
  ProbParm::L_y = probhi[1] - problo[1];
  ProbParm::L_z = probhi[2] - problo[2];

  // Wavelength associated to Taylor length scale
  ProbParm::k0 = 2.0 / ProbParm::lambda0;

  // Initial density, velocity, and material properties
  amrex::Real cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(
    ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, ProbParm::eint0);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);

  ProbParm::urms0 = ProbParm::mach_t0 * cs / sqrt(3.0);
  ProbParm::tau = ProbParm::lambda0 / ProbParm::urms0;

  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity = ProbParm::rho0 * ProbParm::urms0 *
                                      ProbParm::lambda0 /
                                      ProbParm::reynolds_lambda0;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / ProbParm::prandtl;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs)
    << "lambda0, k0, rho0, urms0, tau, p0, T0, gamma, mu, k, c_s0, Reynolds, "
       "Mach, Prandtl, u0, v0, w0, forcing"
    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::lambda0 << "," << ProbParm::k0 << "," << ProbParm::rho0 << ","
    << ProbParm::urms0 << "," << ProbParm::tau << "," << ProbParm::p0 << ","
    << ProbParm::T0 << "," << EOS::gamma << ","
    << transport_params::const_viscosity << ","
    << transport_params::const_conductivity << "," << cs << ","
    << ProbParm::reynolds_lambda0 << "," << ProbParm::mach_t0 << ","
    << ProbParm::prandtl << "," << forcing_params::u0 << ","
    << forcing_params::v0 << "," << forcing_params::w0 << ","
    << forcing_params::forcing << std::endl;
  ofs.close();

  // Load velocity fields from file. Assume data set ordered in Fortran
  // format and reshape the data accordingly. One thing to keep in mind
  // is that this contains the entire input data. We will interpolate
  // this data later to just match our box. Another assumption is that
  // the input data is a periodic cube. If the input cube is smaller
  // than our domain size, the cube will be repeated throughout the
  // domain (hence the mod operations in the interpolation).
  if (ProbParm::restart) {
    amrex::Print() << "Skipping input file reading and assuming restart."
                   << std::endl;
  } else {
    const size_t nx = ProbParm::inres;
    const size_t ny = ProbParm::inres;
    const size_t nz = ProbParm::inres;
    amrex::Vector<double> data(nx * ny * nz * 6); /* this needs to be double */
    if (ProbParm::binfmt) {
      read_binary(ProbParm::iname, nx, ny, nz, 6, data);
    } else {
      read_csv(ProbParm::iname, nx, ny, nz, data);
    }

    // Extract position and velocities
    ProbParm::v_uinput = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_vinput = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_winput = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_xarray = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_xinput = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_xdiff = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::v_xinput->resize(nx * ny * nz);
    ProbParm::v_uinput->resize(nx * ny * nz);
    ProbParm::v_vinput->resize(nx * ny * nz);
    ProbParm::v_winput->resize(nx * ny * nz);
    for (unsigned long i = 0; i < ProbParm::v_xinput->size(); i++) {
      (*ProbParm::v_xinput)[i] = data[0 + i * 6];
      (*ProbParm::v_uinput)[i] =
        data[3 + i * 6] * ProbParm::urms0 / ProbParm::uin_norm;
      (*ProbParm::v_vinput)[i] =
        data[4 + i * 6] * ProbParm::urms0 / ProbParm::uin_norm;
      (*ProbParm::v_winput)[i] =
        data[5 + i * 6] * ProbParm::urms0 / ProbParm::uin_norm;
    }

    // Get the xarray table and the differences.
    ProbParm::v_xarray->resize(nx);
    for (unsigned long i = 0; i < ProbParm::v_xarray->size(); i++) {
      (*ProbParm::v_xarray)[i] = (*ProbParm::v_xinput)[i];
    }
    ProbParm::v_xdiff->resize(nx);
    std::adjacent_difference(
      ProbParm::v_xarray->begin(), ProbParm::v_xarray->end(),
      ProbParm::v_xdiff->begin());
    (*ProbParm::v_xdiff)[0] = (*ProbParm::v_xdiff)[1];

    // Make sure the search array is increasing
    if (not std::is_sorted(
          ProbParm::v_xarray->begin(), ProbParm::v_xarray->end()))
      amrex::Abort("Error: non ascending x-coordinate array.");

    // Get pointer to the data
    ProbParm::xinput = ProbParm::v_xinput->dataPtr();
    ProbParm::uinput = ProbParm::v_uinput->dataPtr();
    ProbParm::vinput = ProbParm::v_vinput->dataPtr();
    ProbParm::winput = ProbParm::v_winput->dataPtr();
    ProbParm::xarray = ProbParm::v_xarray->dataPtr();
    ProbParm::xdiff = ProbParm::v_xdiff->dataPtr();

    // Dimensions of the input box.
    ProbParm::Linput =
      (*ProbParm::v_xarray)[nx - 1] + 0.5 * (*ProbParm::v_xdiff)[nx - 1];
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
