#include "Transport.H"

void
pc_transport_init()
{
  transport_params::init();
}

void
pc_transport_close()
{
  transport_params::finalize();
}

AMREX_GPU_DEVICE
void
pc_get_transport_coeffs(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& D)
{

  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);

  bool wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag;

  wtr_get_xi = true;
  wtr_get_mu = true;
  wtr_get_lam = true;
  wtr_get_Ddiag = true;

  amrex::Real T;
  amrex::Real rho;
  amrex::Real massloc[NUM_SPECIES];

  amrex::Real muloc, xiloc, lamloc;
  amrex::Real Ddiag[NUM_SPECIES];

  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {

        T = q(i, j, k, QTEMP);
        rho = q(i, j, k, QRHO);
        for (int n = 0; n < NUM_SPECIES; ++n) {
          massloc[n] = q(i, j, k, QFS + n);
        }

        pc_actual_transport(
          wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, T, rho, massloc,
          Ddiag, muloc, xiloc, lamloc);

        //   mu, xi and lambda are stored after D in the diffusion multifab
        for (int n = 0; n < NUM_SPECIES; ++n) {
          D(i, j, k, n) = Ddiag[n];
        }

        D(i, j, k, dComp_mu) = muloc;
        D(i, j, k, dComp_xi) = xiloc;
        D(i, j, k, dComp_lambda) = lamloc;
      }
    }
  }
}
