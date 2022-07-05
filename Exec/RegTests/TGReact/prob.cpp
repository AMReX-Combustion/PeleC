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
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho0", PeleC::h_prob_parm_device->rho0);
  pp.query("u0", PeleC::h_prob_parm_device->u0);
  pp.query("mu", PeleC::h_prob_parm_device->mu);
  pp.query("iname", PeleC::prob_parm_host->iname);
  pp.query("case_type", PeleC::prob_parm_host->case_type);
  pp.query("nx", PeleC::h_prob_parm_device->nx);
  pp.query("nvars", PeleC::h_prob_parm_device->nvars);

  // Convert case to int (can't have strings on GPU)
  if (PeleC::prob_parm_host->case_type == "cold") {
    PeleC::h_prob_parm_device->case_type_int = 0;
  } else if (PeleC::prob_parm_host->case_type == "nonreacting") {
    PeleC::h_prob_parm_device->case_type_int = 1;
  } else if (PeleC::prob_parm_host->case_type == "reacting") {
    PeleC::h_prob_parm_device->case_type_int = 2;
  }

  // Define the length scale
  PeleC::h_prob_parm_device->L_x = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L_y = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L_z = probhi[2] - problo[2];

  // Load x profiles from file
  if (
    (PeleC::h_prob_parm_device->case_type_int == 1) ||
    (PeleC::h_prob_parm_device->case_type_int == 2)) {
    amrex::Vector<double> data(
      PeleC::h_prob_parm_device->nx *
      PeleC::h_prob_parm_device->nvars); /* this needs to be double */
    read_csv(
      PeleC::prob_parm_host->iname, PeleC::h_prob_parm_device->nx, 1, 1, data);

    PeleC::prob_parm_host->h_input.resize(
      PeleC::h_prob_parm_device->nx * PeleC::h_prob_parm_device->nvars);
    PeleC::prob_parm_host->h_dxinput.resize(PeleC::h_prob_parm_device->nx);
    PeleC::prob_parm_host->h_xarray.resize(PeleC::h_prob_parm_device->nx);
    for (int i = 0; i < PeleC::prob_parm_host->h_input.size(); i++) {
      PeleC::prob_parm_host->h_input[i] = data[i];
    }
    const amrex::Real m2cm = 100.0;
    for (int i = 0; i < PeleC::h_prob_parm_device->nx; i++) {
      PeleC::prob_parm_host->h_input
        [PeleC::h_prob_parm_device->input_x +
         i * PeleC::h_prob_parm_device->nvars] *= m2cm;
      PeleC::prob_parm_host->h_xarray[i] =
        PeleC::prob_parm_host->h_input
          [PeleC::h_prob_parm_device->input_x +
           i * PeleC::h_prob_parm_device->nvars];
    }
    std::adjacent_difference(
      PeleC::prob_parm_host->h_xarray.begin(),
      PeleC::prob_parm_host->h_xarray.end(),
      PeleC::prob_parm_host->h_dxinput.begin());
    PeleC::prob_parm_host->h_dxinput[0] = PeleC::prob_parm_host->h_dxinput[1];

    // Get pointer to the data
    PeleC::prob_parm_host->input.resize(PeleC::prob_parm_host->h_input.size());
    PeleC::prob_parm_host->xarray.resize(
      PeleC::prob_parm_host->h_xarray.size());
    PeleC::prob_parm_host->dxinput.resize(
      PeleC::prob_parm_host->h_dxinput.size());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_input.begin(),
      PeleC::prob_parm_host->h_input.end(),
      PeleC::prob_parm_host->input.begin());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_xarray.begin(),
      PeleC::prob_parm_host->h_xarray.end(),
      PeleC::prob_parm_host->xarray.begin());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_dxinput.begin(),
      PeleC::prob_parm_host->h_dxinput.end(),
      PeleC::prob_parm_host->dxinput.begin());
    PeleC::h_prob_parm_device->d_input = PeleC::prob_parm_host->input.data();
    PeleC::h_prob_parm_device->d_xarray = PeleC::prob_parm_host->xarray.data();
    PeleC::h_prob_parm_device->d_dxinput =
      PeleC::prob_parm_host->dxinput.data();

    // Write IC to file
    const amrex::Real tau =
      PeleC::h_prob_parm_device->L_x / PeleC::h_prob_parm_device->u0;

    // Output IC
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "L, tau, u0, p0, T0, omega_x, omega_y, omega_z"
                      << std::endl;
    amrex::Print(ofs).SetPrecision(17)
      << PeleC::h_prob_parm_device->L_x << "," << tau << ","
      << PeleC::h_prob_parm_device->u0 << "," << PeleC::h_prob_parm_device->p0
      << "," << PeleC::h_prob_parm_device->T0 << ","
      << PeleC::h_prob_parm_device->omega_x << ","
      << PeleC::h_prob_parm_device->omega_y << ","
      << PeleC::h_prob_parm_device->omega_z << std::endl;
    ofs.close();
  } else {

#ifdef USE_CONSTANT_TRANSPORT
    // Initial density, velocity, and material properties
    amrex::Real cs, cp;
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[0] = 1.0;
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2P(
      PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
      PeleC::h_prob_parm_device->p0);
    eos.RTY2Cs(
      PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
      cs);
    eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);

    auto& trans_parm = PeleC::trans_parms.host_trans_parm();
    trans_parm.const_bulk_viscosity = 0.0;
    trans_parm.const_diffusivity = 0.0;
    trans_parm.const_viscosity = PeleC::h_prob_parm_device->mu;
    const amrex::Real reynolds =
      PeleC::h_prob_parm_device->rho0 * PeleC::h_prob_parm_device->u0 *
      PeleC::h_prob_parm_device->L_x / PeleC::h_prob_parm_device->mu;
    const amrex::Real mach = PeleC::h_prob_parm_device->u0 / cs;
    const amrex::Real prandtl = 0.71;
    trans_parm.const_conductivity = trans_parm.const_viscosity * cp / prandtl;
    PeleC::trans_parms.sync_to_device();
    const amrex::Real tau =
      PeleC::h_prob_parm_device->L_x / PeleC::h_prob_parm_device->u0;

    // Write IC to file
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "L, tau, rho0, u0, p0, T0, gamma, mu, k, c_s0, "
                         "Reynolds, Mach, Prandtl, omega_x, omega_y, omega_z"
                      << std::endl;
    amrex::Print(ofs).SetPrecision(17)
      << PeleC::h_prob_parm_device->L_x << "," << tau << ","
      << PeleC::h_prob_parm_device->rho0 << "," << PeleC::h_prob_parm_device->u0
      << "," << PeleC::h_prob_parm_device->p0 << ","
      << PeleC::h_prob_parm_device->T0 << "," << eos.gamma << ","
      << PeleC::h_prob_parm_device->mu << "," << trans_parm.const_conductivity
      << "," << cs << "," << reynolds << "," << mach << "," << prandtl << ","
      << PeleC::h_prob_parm_device->omega_x << ","
      << PeleC::h_prob_parm_device->omega_y << ","
      << PeleC::h_prob_parm_device->omega_z << std::endl;
    ofs.close();
#else
    amrex::Abort("This needs constant transport");
#endif
  }
}
}

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

void
PeleC::problem_post_restart()
{
}
