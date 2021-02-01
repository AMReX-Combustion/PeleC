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
    pp.query("reynolds", PeleC::prob_parm_device->reynolds);
    pp.query("mach", PeleC::prob_parm_device->mach);
    pp.query("prandtl", PeleC::prob_parm_device->prandtl);
    pp.query("rho_x_fact", PeleC::prob_parm_device->rho_x_fact);
    pp.query("rho_y_fact", PeleC::prob_parm_device->rho_y_fact);
    pp.query("rho_z_fact", PeleC::prob_parm_device->rho_z_fact);
    pp.query("u_0_fact", PeleC::prob_parm_device->u_0_fact);
    pp.query("u_x_fact", PeleC::prob_parm_device->u_x_fact);
    pp.query("u_y_fact", PeleC::prob_parm_device->u_y_fact);
    pp.query("u_z_fact", PeleC::prob_parm_device->u_z_fact);
    pp.query("v_0_fact", PeleC::prob_parm_device->v_0_fact);
    pp.query("v_x_fact", PeleC::prob_parm_device->v_x_fact);
    pp.query("v_y_fact", PeleC::prob_parm_device->v_y_fact);
    pp.query("v_z_fact", PeleC::prob_parm_device->v_z_fact);
    pp.query("w_0_fact", PeleC::prob_parm_device->w_0_fact);
    pp.query("w_x_fact", PeleC::prob_parm_device->w_x_fact);
    pp.query("w_y_fact", PeleC::prob_parm_device->w_y_fact);
    pp.query("w_z_fact", PeleC::prob_parm_device->w_z_fact);
    pp.query("p_x_fact", PeleC::prob_parm_device->p_x_fact);
    pp.query("p_y_fact", PeleC::prob_parm_device->p_y_fact);
    pp.query("p_z_fact", PeleC::prob_parm_device->p_z_fact);
    pp.query("a_rhox", PeleC::prob_parm_device->a_rhox);
    pp.query("a_rhoy", PeleC::prob_parm_device->a_rhoy);
    pp.query("a_rhoz", PeleC::prob_parm_device->a_rhoz);
    pp.query("a_ux", PeleC::prob_parm_device->a_ux);
    pp.query("a_uy", PeleC::prob_parm_device->a_uy);
    pp.query("a_uz", PeleC::prob_parm_device->a_uz);
    pp.query("a_vx", PeleC::prob_parm_device->a_vx);
    pp.query("a_vy", PeleC::prob_parm_device->a_vy);
    pp.query("a_vz", PeleC::prob_parm_device->a_vz);
    pp.query("a_wx", PeleC::prob_parm_device->a_wx);
    pp.query("a_wy", PeleC::prob_parm_device->a_wy);
    pp.query("a_wz", PeleC::prob_parm_device->a_wz);
    pp.query("a_px", PeleC::prob_parm_device->a_px);
    pp.query("a_py", PeleC::prob_parm_device->a_py);
    pp.query("a_pz", PeleC::prob_parm_device->a_pz);
  }

  // Define the length scale
  AMREX_D_TERM(PeleC::prob_parm_device->L_x = probhi[0] - problo[0];
               , PeleC::prob_parm_device->L_y = probhi[1] - problo[1];
               , PeleC::prob_parm_device->L_z = probhi[2] - problo[2];)

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  EOS::PYT2RE(
    PeleC::prob_parm_device->p0, massfrac, PeleC::prob_parm_device->T0,
    PeleC::prob_parm_device->rho0, eint);
  EOS::RTY2Cs(
    PeleC::prob_parm_device->rho0, PeleC::prob_parm_device->T0, massfrac, cs);
  EOS::TY2Cp(PeleC::prob_parm_device->T0, massfrac, cp);

  PeleC::prob_parm_device->u0 = PeleC::prob_parm_device->mach * cs;

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

  trans_parm.const_bulk_viscosity =
    PeleC::prob_parm_device->rho0 * PeleC::prob_parm_device->u0 *
    PeleC::prob_parm_device->L_x / PeleC::prob_parm_device->reynolds;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    PeleC::prob_parm_device->rho0 * PeleC::prob_parm_device->u0 *
    PeleC::prob_parm_device->L_x / PeleC::prob_parm_device->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / PeleC::prob_parm_device->prandtl;

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
  std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif

  // clang-format off
  // MASA parameters for the following functions
  // rho = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
  // u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * u_z * cos(a_uz * PI * z / L);
  // v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * v_z * cos(a_vz * PI * z / L);
  // w = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * w_z * cos(a_wz * PI * z / L);
  // p = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * p_z * cos(a_pz * PI * z / L);
  // clang-format on
  masa_set_param("L", PeleC::prob_parm_device->L_x);
  masa_set_param("R", EOS::RU / EOS::AIRMW);
  masa_set_param("k", trans_parm.const_conductivity);
  masa_set_param("Gamma", EOS::gamma);
  masa_set_param("mu", trans_parm.const_viscosity);
  masa_set_param("mu_bulk", trans_parm.const_bulk_viscosity);
  masa_set_param("rho_0", PeleC::prob_parm_device->rho0);
  masa_set_param(
    "rho_x",
    PeleC::prob_parm_device->rho_x_fact * PeleC::prob_parm_device->rho0);
  masa_set_param("rho_y", PeleC::prob_parm_device->rho_y_fact);
  masa_set_param("rho_z", PeleC::prob_parm_device->rho_z_fact);
  masa_set_param(
    "u_0", PeleC::prob_parm_device->u_0_fact * PeleC::prob_parm_device->u0);
  masa_set_param(
    "u_x", PeleC::prob_parm_device->u_x_fact * PeleC::prob_parm_device->u0);
  masa_set_param("u_y", PeleC::prob_parm_device->u_y_fact);
  masa_set_param("u_z", PeleC::prob_parm_device->u_z_fact);
  masa_set_param(
    "v_0", PeleC::prob_parm_device->v_0_fact * PeleC::prob_parm_device->u0);
  masa_set_param("v_x", PeleC::prob_parm_device->v_x_fact);
  masa_set_param(
    "v_y", PeleC::prob_parm_device->v_y_fact * PeleC::prob_parm_device->u0);
  masa_set_param("v_z", PeleC::prob_parm_device->v_z_fact);
  masa_set_param(
    "w_0", PeleC::prob_parm_device->w_0_fact * PeleC::prob_parm_device->u0);
  masa_set_param("w_x", PeleC::prob_parm_device->w_x_fact);
  masa_set_param("w_y", PeleC::prob_parm_device->w_y_fact);
  masa_set_param(
    "w_z", PeleC::prob_parm_device->w_z_fact * PeleC::prob_parm_device->u0);
  masa_set_param("p_0", PeleC::prob_parm_device->p0);
  masa_set_param(
    "p_x", PeleC::prob_parm_device->p_x_fact * PeleC::prob_parm_device->p0);
  masa_set_param("p_y", PeleC::prob_parm_device->p_y_fact);
  masa_set_param("p_z", PeleC::prob_parm_device->p_z_fact);
  masa_set_param("a_rhox", PeleC::prob_parm_device->a_rhox);
  masa_set_param("a_rhoy", PeleC::prob_parm_device->a_rhoy);
  masa_set_param("a_rhoz", PeleC::prob_parm_device->a_rhoz);
  masa_set_param("a_ux", PeleC::prob_parm_device->a_ux);
  masa_set_param("a_uy", PeleC::prob_parm_device->a_uy);
  masa_set_param("a_uz", PeleC::prob_parm_device->a_uz);
  masa_set_param("a_vx", PeleC::prob_parm_device->a_vx);
  masa_set_param("a_vy", PeleC::prob_parm_device->a_vy);
  masa_set_param("a_vz", PeleC::prob_parm_device->a_vz);
  masa_set_param("a_wx", PeleC::prob_parm_device->a_wx);
  masa_set_param("a_wy", PeleC::prob_parm_device->a_wy);
  masa_set_param("a_wz", PeleC::prob_parm_device->a_wz);
  masa_set_param("a_px", PeleC::prob_parm_device->a_px);
  masa_set_param("a_py", PeleC::prob_parm_device->a_py);
  masa_set_param("a_pz", PeleC::prob_parm_device->a_pz);

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
  if ((verbose <= 0) || (!do_mms)) {
    return;
  }

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

  if ((verbose <= 0) || (!do_mms)) {
    return;
  }

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
