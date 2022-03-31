#include "prob.H"

void
pc_prob_close()
{
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
  // For now, this method is restricted to 3D. Could be portable to 2D in the future with a little bit of work.

  // Set seed for noise. This makes it deterministic and supports restarts.
  // Change to AMReX random functions. Initialize turbulence inflow.
  const amrex::ULong deterministicSeed = 100;
  amrex::ResetRandomSeed(deterministicSeed);

  // We need to determine cell size via ParmParse here since geom data is not passed.
  amrex::Vector<int> n_cells;
  amrex::ParmParse ppamr("amr");
  ppamr.getarr("n_cell", n_cells);

  // consider using the finest level of AMR if we are using it.
  // For now we will just use the base level.
  int max_level = 0;
  ppamr.query("max_level", max_level);
  
  // Parse params
  amrex::ParmParse pp("prob");
  // pp.query("turb_length_scale", PeleC::h_prob_parm_device->turb_length_scale);
  pp.query("turb_intensity", PeleC::h_prob_parm_device->turb_intensity);
  pp.query("u_in", PeleC::h_prob_parm_device->u_in);
  
  // Compute turbulence length scale based on inlet domain size. Assume x-dimension is much larger than y or z.
  // Currently computing as 1/4th of the domain size, can change if necessary
  PeleC::h_prob_parm_device->turb_length_scale = 0.25 * amrex::min<amrex::Real>(probhi[1]-problo[1], probhi[2]-problo[2]);

  // Convert dimensions from PeleC default cgs units to mks units for method
  amrex::Real u_in_ms = PeleC::h_prob_parm_device->u_in / 100;
  amrex::Real L_in_ms = PeleC::h_prob_parm_device->turb_length_scale / 100;

  // largest length scale computed based on turb length scale 
  const amrex::Real k_min = 1.0 / PeleC::h_prob_parm_device->turb_length_scale;

  // smallest length scale. 
  const amrex::Real dx_base = (probhi[1] - problo[1]) / n_cells[1];
  const amrex::Real k_max = 1.0 / (dx_base);

  // Compute dk based on a discretization of M different discrete wavenumbers
  const amrex::Real dk = (k_max - k_min) / PeleC::h_prob_parm_device->turb_num_modes;

  // Initialize a few relevant variables.
  amrex::Vector<amrex::Real> zeta_cross_k(3);
  amrex::Vector<amrex::Real> xi_cross_k(3);
  amrex::Real k_mag, theta, phi, f, E, cross_mag, root_arg, std_omega, randa, omega_upper_bound;

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Checkpoint #1 (before loop) " << '\n';
  }
  
  // Loop to initialize and save arrays for zeta and xi to be used in next loop
  for (int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {
      // populate xi and zeta
      for (int i = 0; i < 3; i++) {
        PeleC::h_prob_parm_device->xi[i][n] = amrex::RandomNormal(0.0, 1.0);
        PeleC::h_prob_parm_device->zeta[i][n] = amrex::RandomNormal(0.0, 1.0);
      }
  }

  // Loop to initialize and save arrays for p, q, k, omega
  for (int m = 0; m < PeleC::h_prob_parm_device->turb_num_modes; m++) {

    // wave number magnitude
    k_mag = k_min + (m * dk);

    // associated frequency. Should this be turbulence velocity???? SHRW to address.
    f = k_mag * u_in_ms;

    // Generate E for this wavenumber
    E = 1.5 * (4.0 * std::pow(PeleC::h_prob_parm_device->turb_intensity * u_in_ms, 2.0)
               * (L_in_ms / u_in_ms)
         / std::pow(1.0 + 70.8 * std::pow(k_mag * L_in_ms, 2.0), 5.0/6.0));


    for (int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {

      // Generate random angles for isotropic sphere
      theta = amrex::Random() * 2.0 * constants::PI();
      phi   = amrex::Random() * constants::PI();

      // Convert to cartesian coords
      PeleC::h_prob_parm_device->k[0][n][m] = k_mag * cos(theta) * sin(phi);
      PeleC::h_prob_parm_device->k[1][n][m] = k_mag * sin(theta) * sin(phi);
      PeleC::h_prob_parm_device->k[2][n][m] = k_mag * cos(phi);

      // Generate random values for getting p and q
      randa = amrex::Random();
      omega_upper_bound = 2.0 * constants::PI() * f;
      PeleC::h_prob_parm_device->omega[n][m] = amrex::RandomNormal(0.0, omega_upper_bound);

      // Calculate p, q

      xi_cross_k[0] = (PeleC::h_prob_parm_device->xi[1][n] * PeleC::h_prob_parm_device->k[2][n][m]
                       - PeleC::h_prob_parm_device->k[1][n][m] * PeleC::h_prob_parm_device->xi[2][n]);
      xi_cross_k[1] = (PeleC::h_prob_parm_device->xi[2][n] * PeleC::h_prob_parm_device->k[0][n][m]
                       - PeleC::h_prob_parm_device->k[2][n][m] * PeleC::h_prob_parm_device->xi[0][n]);
      xi_cross_k[2] = (PeleC::h_prob_parm_device->xi[0][n] * PeleC::h_prob_parm_device->k[1][n][m]
                       - PeleC::h_prob_parm_device->k[0][n][m] * PeleC::h_prob_parm_device->xi[1][n]);

      cross_mag = std::sqrt(std::pow(xi_cross_k[0],2.0) +
                            std::pow(xi_cross_k[1],2.0) + std::pow(xi_cross_k[2],2.0));

      root_arg = randa * 4.0 * E / PeleC::h_prob_parm_device->sampling_number;

      for (int i = 0; i < 3; i++) {
        PeleC::h_prob_parm_device->p[i][n][m] = (xi_cross_k[i] / cross_mag) * std::sqrt(root_arg);
      }
      

      zeta_cross_k[0] = (PeleC::h_prob_parm_device->zeta[1][n] * PeleC::h_prob_parm_device->k[2][n][m]
                       - PeleC::h_prob_parm_device->k[1][n][m] * PeleC::h_prob_parm_device->zeta[2][n]);
      zeta_cross_k[1] = (PeleC::h_prob_parm_device->zeta[2][n] * PeleC::h_prob_parm_device->k[0][n][m]
                       - PeleC::h_prob_parm_device->k[2][n][m] * PeleC::h_prob_parm_device->zeta[0][n]);
      zeta_cross_k[2] = (PeleC::h_prob_parm_device->zeta[0][n] * PeleC::h_prob_parm_device->k[1][n][m]
                       - PeleC::h_prob_parm_device->k[0][n][m] * PeleC::h_prob_parm_device->zeta[1][n]);

      cross_mag = std::sqrt(std::pow(zeta_cross_k[0],2.0) +
                            std::pow(zeta_cross_k[1],2.0) + std::pow(zeta_cross_k[2],2.0));

      root_arg = (1.0 - randa) * 4.0 * E / PeleC::h_prob_parm_device->sampling_number;
      // amrex::Print() << "Root Arg = " << root_arg << '\n' << "Random# = " << randa << '\n' << "E = " << E << '\n';
      for (int i = 0; i < 3; i++) {
        PeleC::h_prob_parm_device->q[i][n][m] = (zeta_cross_k[i] / cross_mag) * std::sqrt(root_arg);
      }

      // Normalize k to k_tilde using the lowest wave number (equation 21 Huang 2010)
      for (int i = 0; i < 3; i++) {
	      PeleC::h_prob_parm_device->k[i][n][m] /= k_min;

      }

    }

  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Checkpoint #4 (after loop) " << '\n';
  }
  
  pp.query("p_init", PeleC::h_prob_parm_device->p_init);
  pp.query("T_init", PeleC::h_prob_parm_device->T_init);
  pp.query("phi_init", PeleC::h_prob_parm_device->phi_in);

  pp.query("p_in", PeleC::h_prob_parm_device->p_in);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);

  for (int n = 0; n < NUM_SPECIES; n++)
    PeleC::h_prob_parm_device->molefrac[n] = 0.0;

  // for CH4-air
  const amrex::Real a = 5.0;
  
  PeleC::h_prob_parm_device->molefrac[O2_ID] = 1.0 / (1.0 + (PeleC::h_prob_parm_device->phi_in / a)  + (0.79 / 0.21));
  PeleC::h_prob_parm_device->molefrac[C3H8_ID] = PeleC::h_prob_parm_device->phi_in * PeleC::h_prob_parm_device->molefrac[O2_ID] / a;
  PeleC::h_prob_parm_device->molefrac[N2_ID] = 1.0 - PeleC::h_prob_parm_device->molefrac[C3H8_ID] - PeleC::h_prob_parm_device->molefrac[O2_ID];

  // for initializing the domain with "engineering air"
  // ProbParm::molefrac_init[O2_ID] = 0.21;
  // ProbParm::molefrac_init[N2_ID] = 0.79;

  // Convert mole fracs to mass fracs
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(PeleC::h_prob_parm_device->molefrac.begin(), PeleC::h_prob_parm_device->massfrac.begin());

  // Initialize density and energy from mass fractions, T and P.
  eos.PYT2RE(
	     PeleC::h_prob_parm_device->p_in, PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T_in,
	     PeleC::h_prob_parm_device->rho_in, PeleC::h_prob_parm_device->e_in);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Checkpoint #5 (end of probinit) " << '\n';
  }
}
}


void
PeleC::problem_post_timestep()
{
   if (amrex::ParallelDescriptor::IOProcessor()) {
       amrex::Print() << "TIME= " << time << " p[0]  = " << PeleC::h_prob_parm_device->p[0][5][10] << '\n';
       amrex::Print() << "TIME= " << time << " q[0]  = " << PeleC::h_prob_parm_device->q[0][5][10] << '\n';
       amrex::Print() << "TIME= " << time << " k[0]  = " << PeleC::h_prob_parm_device->k[0][5][10] << '\n';
       amrex::Print() << "TIME= " << time << " k[1]  = " << PeleC::h_prob_parm_device->k[1][5][10] << '\n';
       amrex::Print() << "TIME= " << time << " omega  = " << PeleC::h_prob_parm_device->omega[5][10] << '\n';
   }
 }


void
PeleC::problem_post_init()
{
}


void
PeleC::problem_post_restart()
{
}
