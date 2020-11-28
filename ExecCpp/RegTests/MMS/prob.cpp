#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real reynolds = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real mach = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real prandtl = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_x_fact = 0.1;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_y_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_z_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u_0_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u_x_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u_y_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u_z_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v_0_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v_x_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v_y_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v_z_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real w_0_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real w_x_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real w_y_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real w_z_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p_x_fact = 0.2;
AMREX_GPU_DEVICE_MANAGED amrex::Real p_y_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p_z_fact = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhox = 2.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoy = 4.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_rhoz = 6.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_ux = 2.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_uy = 4.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_uz = 6.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_vx = 4.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_vy = 2.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_vz = 6.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_wx = 6.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_wy = 4.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_wz = 2.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_px = 6.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_py = 2.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real a_pz = 4.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_x = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_y = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real L_z = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real u0 = 0.0;
} // namespace ProbParm

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
  amrex::ParmParse pp("prob");
  pp.query("reynolds", ProbParm::reynolds);
  pp.query("mach", ProbParm::mach);
  pp.query("prandtl", ProbParm::prandtl);
  pp.query("rho_x_fact", ProbParm::rho_x_fact);
  pp.query("rho_y_fact", ProbParm::rho_y_fact);
  pp.query("rho_z_fact", ProbParm::rho_z_fact);
  pp.query("u_0_fact", ProbParm::u_0_fact);
  pp.query("u_x_fact", ProbParm::u_x_fact);
  pp.query("u_y_fact", ProbParm::u_y_fact);
  pp.query("u_z_fact", ProbParm::u_z_fact);
  pp.query("v_0_fact", ProbParm::v_0_fact);
  pp.query("v_x_fact", ProbParm::v_x_fact);
  pp.query("v_y_fact", ProbParm::v_y_fact);
  pp.query("v_z_fact", ProbParm::v_z_fact);
  pp.query("w_0_fact", ProbParm::w_0_fact);
  pp.query("w_x_fact", ProbParm::w_x_fact);
  pp.query("w_y_fact", ProbParm::w_y_fact);
  pp.query("w_z_fact", ProbParm::w_z_fact);
  pp.query("p_x_fact", ProbParm::p_x_fact);
  pp.query("p_y_fact", ProbParm::p_y_fact);
  pp.query("p_z_fact", ProbParm::p_z_fact);
  pp.query("a_rhox", ProbParm::a_rhox);
  pp.query("a_rhoy", ProbParm::a_rhoy);
  pp.query("a_rhoz", ProbParm::a_rhoz);
  pp.query("a_ux", ProbParm::a_ux);
  pp.query("a_uy", ProbParm::a_uy);
  pp.query("a_uz", ProbParm::a_uz);
  pp.query("a_vx", ProbParm::a_vx);
  pp.query("a_vy", ProbParm::a_vy);
  pp.query("a_vz", ProbParm::a_vz);
  pp.query("a_wx", ProbParm::a_wx);
  pp.query("a_wy", ProbParm::a_wy);
  pp.query("a_wz", ProbParm::a_wz);
  pp.query("a_px", ProbParm::a_px);
  pp.query("a_py", ProbParm::a_py);
  pp.query("a_pz", ProbParm::a_pz);

  // Define the length scale
  AMREX_D_TERM(ProbParm::L_x = probhi[0] - problo[0];
               , ProbParm::L_y = probhi[1] - problo[1];
               , ProbParm::L_z = probhi[2] - problo[2];)

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, eint);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);

  ProbParm::u0 = ProbParm::mach * cs;
  transport_params::const_bulk_viscosity =
    ProbParm::rho0 * ProbParm::u0 * ProbParm::L_x / ProbParm::reynolds;
  transport_params::const_diffusivity = 0.0;
  transport_params::const_viscosity =
    ProbParm::rho0 * ProbParm::u0 * ProbParm::L_x / ProbParm::reynolds;
  transport_params::const_conductivity =
    transport_params::const_viscosity * cp / ProbParm::prandtl;

  // clang-format off
  // MASA parameters for the following functions
  // rho = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  // u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * u_z * cos(a_uz * PI * z / L);
  // v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * v_z * cos(a_vz * PI * z / L);
  // w = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * w_z * cos(a_wz * PI * z / L);
  // p = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * p_z * cos(a_pz * PI * z / L);
  // clang-format on
  masa_set_param("L", ProbParm::L_x);
  masa_set_param("R", EOS::RU / EOS::AIRMW);
  masa_set_param("k", transport_params::const_conductivity);
  masa_set_param("Gamma", EOS::gamma);
  masa_set_param("mu", transport_params::const_viscosity);
  masa_set_param("mu_bulk", transport_params::const_bulk_viscosity);
  masa_set_param("rho_0", ProbParm::rho0);
  masa_set_param("rho_x", ProbParm::rho_x_fact * ProbParm::rho0);
  masa_set_param("rho_y", ProbParm::rho_y_fact);
  masa_set_param("rho_z", ProbParm::rho_z_fact);
  masa_set_param("u_0", ProbParm::u_0_fact * ProbParm::u0);
  masa_set_param("u_x", ProbParm::u_x_fact * ProbParm::u0);
  masa_set_param("u_y", ProbParm::u_y_fact);
  masa_set_param("u_z", ProbParm::u_z_fact);
  masa_set_param("v_0", ProbParm::v_0_fact * ProbParm::u0);
  masa_set_param("v_x", ProbParm::v_x_fact);
  masa_set_param("v_y", ProbParm::v_y_fact * ProbParm::u0);
  masa_set_param("v_z", ProbParm::v_z_fact);
  masa_set_param("w_0", ProbParm::w_0_fact * ProbParm::u0);
  masa_set_param("w_x", ProbParm::w_x_fact);
  masa_set_param("w_y", ProbParm::w_y_fact);
  masa_set_param("w_z", ProbParm::w_z_fact * ProbParm::u0);
  masa_set_param("p_0", ProbParm::p0);
  masa_set_param("p_x", ProbParm::p_x_fact * ProbParm::p0);
  masa_set_param("p_y", ProbParm::p_y_fact);
  masa_set_param("p_z", ProbParm::p_z_fact);
  masa_set_param("a_rhox", ProbParm::a_rhox);
  masa_set_param("a_rhoy", ProbParm::a_rhoy);
  masa_set_param("a_rhoz", ProbParm::a_rhoz);
  masa_set_param("a_ux", ProbParm::a_ux);
  masa_set_param("a_uy", ProbParm::a_uy);
  masa_set_param("a_uz", ProbParm::a_uz);
  masa_set_param("a_vx", ProbParm::a_vx);
  masa_set_param("a_vy", ProbParm::a_vy);
  masa_set_param("a_vz", ProbParm::a_vz);
  masa_set_param("a_wx", ProbParm::a_wx);
  masa_set_param("a_wy", ProbParm::a_wy);
  masa_set_param("a_wz", ProbParm::a_wz);
  masa_set_param("a_px", ProbParm::a_px);
  masa_set_param("a_py", ProbParm::a_py);
  masa_set_param("a_pz", ProbParm::a_pz);

  // Display and check
  if (amrex::ParallelDescriptor::IOProcessor()) {
    masa_display_param();
  }
  masa_sanity_check();
}
}

void
PeleC::problem_post_timestep()
{
  if ((verbose <= 0) || (!do_mms))
    return;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real rho_mms_err = 0.0;
  amrex::Real u_mms_err = 0.0;
  amrex::Real v_mms_err = 0.0;
  amrex::Real w_mms_err = 0.0;
  amrex::Real p_mms_err = 0.0;
  amrex::Real rho_residual = 0.0;
  amrex::Real rhou_residual = 0.0;
  amrex::Real rhov_residual = 0.0;
  amrex::Real rhow_residual = 0.0;
  amrex::Real rhoE_residual = 0.0;

#ifdef PELEC_USE_MASA
  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... MMS problem post timestep" << std::endl;
    }

    // Calculate the errors and residuals
    for (int lev = 0; lev <= finest_level; lev++) {
      PeleC& pc_lev = getLevel(lev);

      bool local_flag = true;
      rho_mms_err += pc_lev.volWgtSquaredSum("rhommserror", time, local_flag);
      AMREX_D_TERM(
        u_mms_err += pc_lev.volWgtSquaredSum("ummserror", time, local_flag);
        , v_mms_err += pc_lev.volWgtSquaredSum("vmmserror", time, local_flag);
        , w_mms_err += pc_lev.volWgtSquaredSum("wmmserror", time, local_flag);)
      p_mms_err += pc_lev.volWgtSquaredSum("pmmserror", time, local_flag);

      rho_residual += pc_lev.volWgtSquaredSumDiff(Density, time, local_flag);
      AMREX_D_TERM(
        rhou_residual += pc_lev.volWgtSquaredSumDiff(Xmom, time, local_flag);
        , rhov_residual += pc_lev.volWgtSquaredSumDiff(Ymom, time, local_flag);
        , rhow_residual += pc_lev.volWgtSquaredSumDiff(Zmom, time, local_flag);)
      rhoE_residual += pc_lev.volWgtSquaredSumDiff(Eden, time, local_flag);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    AMREX_D_TERM(
      amrex::ParallelDescriptor::ReduceRealSum(
        &u_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
      , amrex::ParallelDescriptor::ReduceRealSum(
          &v_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
      , amrex::ParallelDescriptor::ReduceRealSum(
          &w_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());)
    amrex::ParallelDescriptor::ReduceRealSum(
      &p_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    AMREX_D_TERM(
      amrex::ParallelDescriptor::ReduceRealSum(
        &rhou_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
      , amrex::ParallelDescriptor::ReduceRealSum(
          &rhov_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
      , amrex::ParallelDescriptor::ReduceRealSum(
          &rhow_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());)
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhoE_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Get the norm and normalize it
    amrex::Real V = volume.sum(0, false);
    rho_mms_err = std::sqrt(rho_mms_err / V);
    AMREX_D_TERM(u_mms_err = std::sqrt(u_mms_err / V);
                 , v_mms_err = std::sqrt(v_mms_err / V);
                 , w_mms_err = std::sqrt(w_mms_err / V);)
    p_mms_err = std::sqrt(p_mms_err / V);
    rho_residual = std::sqrt(rho_residual / V);
    rhou_residual = std::sqrt(rhou_residual / V);
    rhov_residual = std::sqrt(rhov_residual / V);
    rhow_residual = std::sqrt(rhow_residual / V);
    rhoE_residual = std::sqrt(rhoE_residual / V);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " RHO MMS ERROR  = " << rho_mms_err
                     << '\n';
      AMREX_D_TERM(amrex::Print() << "TIME= " << time
                                  << " U MMS ERROR    = " << u_mms_err << '\n';
                   , amrex::Print() << "TIME= " << time << " V MMS ERROR    = "
                                    << v_mms_err << '\n';
                   , amrex::Print() << "TIME= " << time << " W MMS ERROR    = "
                                    << w_mms_err << '\n';)
      amrex::Print() << "TIME= " << time << " P MMS ERROR    = " << p_mms_err
                     << '\n';
      amrex::Print() << "TIME= " << time << " RHO RESIDUAL   = " << rho_residual
                     << '\n';
      AMREX_D_TERM(amrex::Print() << "TIME= " << time << " RHO*U RESIDUAL = "
                                  << rhou_residual << '\n';
                   , amrex::Print() << "TIME= " << time << " RHO*V RESIDUAL = "
                                    << rhov_residual << '\n';
                   , amrex::Print() << "TIME= " << time << " RHO*W RESIDUAL = "
                                    << rhow_residual << '\n';)
      amrex::Print() << "TIME= " << time
                     << " RHO*E RESIDUAL = " << rhoE_residual << '\n';

      if (parent->NumDataLogs() > 1) {

        std::ostream& data_log2 = parent->DataLog(1);

        // Write the quantities at this time
        const int datwidth = 14;
        const int datprecision = 6;
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_mms_err;
        AMREX_D_TERM(
          data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                    << u_mms_err;
          , data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                      << v_mms_err;
          , data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                      << w_mms_err;)
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << p_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_residual;
        AMREX_D_TERM(
          data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                    << rhou_residual;
          , data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                      << rhov_residual;
          , data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                      << rhow_residual;)
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhoE_residual;
        data_log2 << std::endl;
      }
    }
  }

#else
  Error("MASA is not turned on. Turn on with USE_MASA=TRUE.");
#endif
}

void
PeleC::problem_post_init()
{

  if ((verbose <= 0) || (!do_mms))
    return;

  amrex::Real time = state[State_Type].curTime();

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {

      if (parent->NumDataLogs() > 1) {

        std::ostream& data_log2 = parent->DataLog(1);
        if (time == 0.0) {
          const int datwidth = 14;
          data_log2 << std::setw(datwidth) << "          time";
          data_log2 << std::setw(datwidth) << "   rho_mms_err";
          AMREX_D_TERM(data_log2 << std::setw(datwidth) << "     u_mms_err";
                       , data_log2 << std::setw(datwidth) << "     v_mms_err";
                       , data_log2 << std::setw(datwidth) << "     w_mms_err";)
          data_log2 << std::setw(datwidth) << "     p_mms_err";
          data_log2 << std::setw(datwidth) << "  rho_residual";
          AMREX_D_TERM(data_log2 << std::setw(datwidth) << " rhou_residual";
                       , data_log2 << std::setw(datwidth) << " rhov_residual";
                       , data_log2 << std::setw(datwidth) << " rhow_residual";)
          data_log2 << std::setw(datwidth) << " rhoE_residual";
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
