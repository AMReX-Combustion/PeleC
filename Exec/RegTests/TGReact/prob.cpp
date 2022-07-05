#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 1.0e-3;
AMREX_GPU_DEVICE_MANAGED amrex::Real u0 = 100.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real mu = 0.00625;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_x = 2.0 * PI;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_y = 2.0 * PI;
AMREX_GPU_DEVICE_MANAGED amrex::Real omega_z = 2.0 * PI;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
std::string case_type = "cold";
AMREX_GPU_DEVICE_MANAGED int case_type_int = 0;
std::string iname = "";
AMREX_GPU_DEVICE_MANAGED int nx = 0;
AMREX_GPU_DEVICE_MANAGED int nvars = 15;
amrex::Gpu::ManagedVector<amrex::Real>* v_input = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_xarray = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* v_dxinput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* input = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* dxinput = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* xarray = nullptr;
AMREX_GPU_DEVICE_MANAGED int input_x = 0;
AMREX_GPU_DEVICE_MANAGED int input_H2_init = 1;
AMREX_GPU_DEVICE_MANAGED int input_O2_init = 2;
AMREX_GPU_DEVICE_MANAGED int input_N2_init = 3;
AMREX_GPU_DEVICE_MANAGED int input_Tad = 4;
AMREX_GPU_DEVICE_MANAGED int input_p = 5;
AMREX_GPU_DEVICE_MANAGED int input_H2 = 6;
AMREX_GPU_DEVICE_MANAGED int input_O2 = 7;
AMREX_GPU_DEVICE_MANAGED int input_N2 = 8;
AMREX_GPU_DEVICE_MANAGED int input_H2O = 9;
AMREX_GPU_DEVICE_MANAGED int input_H = 10;
AMREX_GPU_DEVICE_MANAGED int input_O = 11;
AMREX_GPU_DEVICE_MANAGED int input_OH = 12;
AMREX_GPU_DEVICE_MANAGED int input_HO2 = 13;
AMREX_GPU_DEVICE_MANAGED int input_H2O2 = 14;
} // namespace ProbParm

void
pc_prob_close()
{
  delete ProbParm::v_input;
  delete ProbParm::v_dxinput;
  delete ProbParm::v_xarray;

  ProbParm::v_input = nullptr;
  ProbParm::v_dxinput = nullptr;
  ProbParm::v_xarray = nullptr;
  ProbParm::input = nullptr;
  ProbParm::dxinput = nullptr;
  ProbParm::xarray = nullptr;
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
  pp.query("rho0", ProbParm::rho0);
  pp.query("u0", ProbParm::u0);
  pp.query("mu", ProbParm::mu);
  pp.query("iname", ProbParm::iname);
  pp.query("case_type", ProbParm::case_type);
  pp.query("nx", ProbParm::nx);
  pp.query("nvars", ProbParm::nvars);

  // Convert case to int (can't have strings on GPU)
  if (ProbParm::case_type == "cold") {
    ProbParm::case_type_int = 0;
  } else if (ProbParm::case_type == "nonreacting") {
    ProbParm::case_type_int = 1;
  } else if (ProbParm::case_type == "reacting") {
    ProbParm::case_type_int = 2;
  }

  // Define the length scale
  ProbParm::L_x = probhi[0] - problo[0];
  ProbParm::L_y = probhi[1] - problo[1];
  ProbParm::L_z = probhi[2] - problo[2];

  // Load x profiles from file
#ifdef PELEC_USE_REACTIONS
  amrex::Vector<double> data(
    ProbParm::nx * ProbParm::nvars); /* this needs to be double */
  read_csv(ProbParm::iname, ProbParm::nx, 1, 1, data);

  ProbParm::v_input = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::v_dxinput = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::v_xarray = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::v_input->resize(ProbParm::nx * ProbParm::nvars);
  ProbParm::v_dxinput->resize(ProbParm::nx);
  ProbParm::v_xarray->resize(ProbParm::nx);
  for (int i = 0; i < ProbParm::v_input->size(); i++) {
    (*ProbParm::v_input)[i] = data[i];
  }
  const amrex::Real m2cm = 100.0;
  for (int i = 0; i < ProbParm::nx; i++) {
    (*ProbParm::v_input)[ProbParm::input_x + i * ProbParm::nvars] *= m2cm;
    (*ProbParm::v_xarray)[i] =
      (*ProbParm::v_input)[ProbParm::input_x + i * ProbParm::nvars];
  }
  std::adjacent_difference(
    ProbParm::v_xarray->begin(), ProbParm::v_xarray->end(),
    ProbParm::v_dxinput->begin());
  (*ProbParm::v_dxinput)[0] = (*ProbParm::v_dxinput)[1];

  // Get pointer to the data
  ProbParm::input = ProbParm::v_input->dataPtr();
  ProbParm::dxinput = ProbParm::v_dxinput->dataPtr();
  ProbParm::xarray = ProbParm::v_xarray->dataPtr();

  // Write IC to file
  const amrex::Real tau = ProbParm::L_x / ProbParm::u0;

  // Output IC
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, tau, u0, p0, T0, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::L_x << "," << tau << "," << ProbParm::u0 << "," << ProbParm::p0
    << "," << ProbParm::T0 << "," << ProbParm::omega_x << ","
    << ProbParm::omega_y << "," << ProbParm::omega_z << std::endl;
  ofs.close();

#else

  // Initial density, velocity, and material properties
  amrex::Real cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;
  EOS::RTY2P(ProbParm::rho0, ProbParm::T0, massfrac, ProbParm::p0);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);

  transport_params::const_bulk_viscosity = 0.0;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity = ProbParm::mu;
  const amrex::Real reynolds =
    ProbParm::rho0 * ProbParm::u0 * ProbParm::L_x / ProbParm::mu;
  const amrex::Real mach = ProbParm::u0 / cs;
  const amrex::Real prandtl = 0.71;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / prandtl;
  const amrex::Real tau = ProbParm::L_x / ProbParm::u0;

  // Write IC to file
  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "L, tau, rho0, u0, p0, T0, gamma, mu, k, c_s0, "
                       "Reynolds, Mach, Prandtl, omega_x, omega_y, omega_z"
                    << std::endl;
  amrex::Print(ofs).SetPrecision(17)
    << ProbParm::L_x << "," << tau << "," << ProbParm::rho0 << ","
    << ProbParm::u0 << "," << ProbParm::p0 << "," << ProbParm::T0 << ","
    << EOS::gamma << "," << ProbParm::mu << ","
    << transport_params::const_conductivity << "," << cs << "," << reynolds
    << "," << mach << "," << prandtl << "," << ProbParm::omega_x << ","
    << ProbParm::omega_y << "," << ProbParm::omega_z << std::endl;
  ofs.close();
#endif
}
}

#ifdef DO_PROBLEM_POST_TIMESTEP
void
PeleC::problem_post_timestep()
{

  if ((verbose <= 0))
    return;

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real max_temp = 0.0;
  int datwidth = 14;
  int datprecision = 6;

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Problem post timestep" << std::endl;
    }

    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      max_temp = maxDerive("Temp", time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealMax(
      &max_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " MAX TEMP  = " << max_temp << '\n';

      if (parent->NumDataLogs() > 1) {

        std::ostream& data_log2 = parent->DataLog(1);

        // Write the quantities at this time
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << max_temp;
        data_log2 << std::endl;
      }
    }
  }
}
#endif

#ifdef DO_PROBLEM_POST_INIT
void
PeleC::problem_post_init()
{

  if ((verbose <= 0))
    return;

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real max_temp = 0.0;
  int datwidth = 14;
  int datprecision = 6;

  if (level == 0) {
    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      max_temp = maxDerive("Temp", time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealMax(
      &max_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {

      if (parent->NumDataLogs() > 1) {

        std::ostream& data_log2 = parent->DataLog(1);
        if (time == 0.0) {
          data_log2 << std::setw(datwidth) << "          time";
          data_log2 << std::setw(datwidth) << "       maxtemp";
          data_log2 << std::endl;
          data_log2 << std::setw(datwidth) << time;
          data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                    << max_temp;
          data_log2 << std::endl;
        }
      }
    }
  }
}
#endif
