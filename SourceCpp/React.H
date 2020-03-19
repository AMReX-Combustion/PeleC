#ifndef _REACT_H_
#define _REACT_H_

#include <iomanip>

#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "EOS.H"

#ifndef AMREX_USE_GPU
const amrex::Real alpha_rk64[6] = {
  0.218150805229859,  //            3296351145737.0/15110423921029.0,
  0.256702469801519,  //            1879360555526.0/ 7321162733569.0,
  0.527402592007520,  //            10797097731880.0/20472212111779.0,
  0.0484864267224467, //            754636544611.0/15563872110659.0,
  1.24517071533530,   //            3260218886217.0/ 2618290685819.0,
  0.412366034843237,  //            5069185909380.0/12292927838509.0
};

const amrex::Real beta_rk64[6] = {
  -0.113554138044166,  //-1204558336989.0/10607789004752.0,
  -0.215118587818400,  //-3028468927040.0/14078136890693.0,
  -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
  -1.07992686223881,   //-17275898420483.0/15997285579755.0,
  -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
  0.0};

const amrex::Real err_rk64[6] = {
  -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
  0.158481845574980,   // 473021958881.0/ 2984707536468.0,
  -0.0905918835751907, //-947229622805.0/10456009803779.0,
  -0.219084567203338,  //-2921473878215.0/13334914072261.0,
  0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
  0.0426421977505659   // 167623581683.0/ 3930932046784.0
};
#endif

// Timestep adapter
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
adapt_timestep(
  const amrex::Real rk_err[],
  const amrex::Real dt_rk4max,
  amrex::Real& dt,
  const amrex::Real dt_rk4min,
  const amrex::Real tol)
{
  amrex::Real change_factor;
  amrex::Real err = amrex::max(rk_err[0], rk_err[1]);
  for (int n = 2; n < NVAR; n++)
    err = amrex::max(rk_err[n], err);

  if (err < tol) {
    err = amrex::max(err, tol * 1.e-4);
    change_factor = std::sqrt(
      std::sqrt(tol / err)); // sqrt sqrt is faster than pow(var, 0.25)
    dt = amrex::min(dt_rk4max, dt * change_factor);
  } else {
    change_factor = std::pow(tol / err, 0.2);
    dt = amrex::max(dt_rk4min, dt * change_factor);
  }
}

// Do the reactions, here uout and IR change
// Rk integrator
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_expl_reactions(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& uold,
  amrex::Array4<amrex::Real> const& unew,
  amrex::Array4<const amrex::Real> const& asrc,
  amrex::Array4<amrex::Real> const& cost,
  amrex::Array4<amrex::Real> const& IR,
  const amrex::Real dt_react,
  const int nsteps_min,
  const int nsteps_max,
  const int nsteps_guess,
  const amrex::Real errtol,
  const int do_update)
{
#ifdef AMREX_USE_GPU
  // having a global __constant__ variable is slower than having this in local
  // memory.
  const amrex::Real alpha_rk64[6] = {
    0.218150805229859,  //            3296351145737.0/15110423921029.0,
    0.256702469801519,  //            1879360555526.0/ 7321162733569.0,
    0.527402592007520,  //            10797097731880.0/20472212111779.0,
    0.0484864267224467, //            754636544611.0/15563872110659.0,
    1.24517071533530,   //            3260218886217.0/ 2618290685819.0,
    0.412366034843237,  //            5069185909380.0/12292927838509.0
  };

  const amrex::Real beta_rk64[6] = {
    -0.113554138044166,  //-1204558336989.0/10607789004752.0,
    -0.215118587818400,  //-3028468927040.0/14078136890693.0,
    -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
    -1.07992686223881,   //-17275898420483.0/15997285579755.0,
    -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
    0.0};

  const amrex::Real err_rk64[6] = {
    -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
    0.158481845574980,   // 473021958881.0/ 2984707536468.0,
    -0.0905918835751907, //-947229622805.0/10456009803779.0,
    -0.219084567203338,  //-2921473878215.0/13334914072261.0,
    0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
    0.0426421977505659   // 167623581683.0/ 3930932046784.0
  };
#endif

  // compute rhoe_ext/rhoy_ext
  amrex::Real rhou = uold(i, j, k, UMX), rhov = uold(i, j, k, UMY),
              rhow = uold(i, j, k, UMZ);
  amrex::Real rho_old = uold(i, j, k, URHO);
  amrex::Real rhoInv = 1.0 / rho_old;

  amrex::Real rho = 0.;

  for (int n = UFS; n < UFS + NUM_SPECIES; n++)
    rho += uold(i, j, k, n);

  amrex::Real nrg =
    (uold(i, j, k, UEDEN) -
     (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) *
    rhoInv;
  rhou = unew(i, j, k, UMX);
  rhov = unew(i, j, k, UMY);
  rhow = unew(i, j, k, UMZ);
  rhoInv = 1.0 / unew(i, j, k, URHO);

  amrex::Real rhoedot_ext =
    ((unew(i, j, k, UEDEN) -
      (0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv)) -
     rho * nrg) /
    dt_react;

  amrex::Real rhoydot_ext[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++)
    rhoydot_ext[n] = asrc(i, j, k, UFS + n);

  // RK dts
  amrex::Real dt_rk = dt_react / nsteps_guess;
  const amrex::Real dt_min = dt_react / nsteps_max;
  const amrex::Real dt_max = dt_react / nsteps_min;
  amrex::Real updt_time = 0.0;
  amrex::Real urk[NVAR];
  amrex::Real urk_old[NVAR] = {};
  for (int n = 0; n < NVAR; ++n)
    urk[n] = uold(i, j, k, n);
  amrex::Real urk_err[NVAR] = {};

  // Do the RK!
  int steps = 0;
  while (updt_time < dt_react) {
    for (int n = 0; n < NVAR; n++) {
      urk_old[n] = urk[n];
      urk_err[n] = 0.0;
    }
    amrex::Real wdot[NUM_SPECIES] = {};
    amrex::Real massfrac[NUM_SPECIES] = {};
    amrex::Real ei[NUM_SPECIES] = {};
    updt_time += dt_rk;
    steps += 1;
    for (int stage = 0; stage < 6; stage++) {
      rhoInv = 1.0 / urk[URHO];
      for (int n = 0; n < NUM_SPECIES; ++n)
        massfrac[n] = urk[UFS + n] * rhoInv;
      EOS::mpr2wdot(
        urk[URHO], urk[UTEMP], massfrac, wdot); // this is wdot*molecularweights
      for (int n = 0; n < NUM_SPECIES; ++n)
        wdot[n] += rhoydot_ext[n];
      EOS::get_ei(urk[UTEMP], ei);
      amrex::Real tempsrc = rhoedot_ext;
      for (int n = 0; n < NUM_SPECIES; ++n)
        tempsrc -= wdot[n] * ei[n];
      amrex::Real cv;
      EOS::get_cv(massfrac, urk[UTEMP], cv);
      tempsrc /= (urk[URHO] * cv);

      /*================== Update urk_err =================== */
      // Species
      for (int n = 0; n < NUM_SPECIES; ++n)
        urk_err[UFS + n] += err_rk64[stage] * dt_rk * wdot[n];
      // Temperature
      urk_err[UTEMP] += err_rk64[stage] * dt_rk * tempsrc;

      /*================== Update Stage solution =================== */
      // Species
      for (int n = 0; n < NUM_SPECIES; ++n) {
        urk[UFS + n] = urk_old[UFS + n] + alpha_rk64[stage] * dt_rk * wdot[n];
      }
      // Temperature
      urk[UTEMP] = urk_old[UTEMP] + alpha_rk64[stage] * dt_rk * tempsrc;

      /*================== Update urk_old =========================== */
      // Species
      for (int n = 0; n < NUM_SPECIES; ++n)
        urk_old[UFS + n] = urk[UFS + n] + beta_rk64[stage] * dt_rk * wdot[n];
      // Temperature
      urk_old[UTEMP] = urk[UTEMP] + beta_rk64[stage] * dt_rk * tempsrc;

      /*================= Update urk[rho] ========================= */
      urk[URHO] = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n)
        urk[URHO] += urk[UFS + n];

      /*================ Adapt Time step! ======================== */
    } // end rk stages
    adapt_timestep(urk_err, dt_max, dt_rk, dt_min, errtol);
  } // end timestep loop
  cost(i, j, k) = steps;

  // Add drhoY/dt to reactions MultiFab and update unew if needed
  amrex::Real umnew = uold(i, j, k, UMX) + dt_react * asrc(i, j, k, UMX);
  amrex::Real vmnew = uold(i, j, k, UMY) + dt_react * asrc(i, j, k, UMY);
  amrex::Real wmnew = uold(i, j, k, UMZ) + dt_react * asrc(i, j, k, UMZ);

  if (do_update) {
    unew(i, j, k, URHO) = urk[URHO];
    unew(i, j, k, UMX) = umnew;
    unew(i, j, k, UMY) = vmnew;
    unew(i, j, k, UMZ) = wmnew;
    unew(i, j, k, UTEMP) = urk[UTEMP];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      unew(i, j, k, UFS + n) = urk[UFS + n];
    }
  }

  for (int n = 0; n < NUM_SPECIES; ++n) {
    IR(i, j, k, n) =
      (urk[UFS + n] - uold(i, j, k, UFS + n)) / dt_react - rhoydot_ext[n];
  }
  IR(i, j, k, NUM_SPECIES) =
    ((nrg * rho_old) + dt_react * rhoedot_ext +
     0.5 * (umnew * umnew + vmnew * vmnew + wmnew * wmnew) / urk[URHO] -
     uold(i, j, k, UEDEN)) /
      dt_react -
    asrc(i, j, k, UEDEN);
  /*    if(i==22 && j == 23 && k == 208){
          std::cout<<std::setprecision(17);
          amrex::Print()<< IR(i,j,k,NUM_SPECIES) << '\t' << (nrg*rho_old) +
     dt_react*rhoedot_ext
                                                  + 0.5*(umnew*umnew +
     vmnew*vmnew + wmnew*wmnew)/urk[URHO] << '\t'
                                                 << -uold(i,j,k,UEDEN) << '\t'
     << -asrc(i,j,k,UEDEN) << std::endl; std::cin.get();
      } */
}

#endif
