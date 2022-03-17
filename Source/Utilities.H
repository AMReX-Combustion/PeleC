#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <AMReX_IArrayBox.H>
#include <AMReX_FArrayBox.H>
#include "Constants.H"
#include "IndexDefines.H"
#include "PelePhysics.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_cmpTemp(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& S)
{
  amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
  amrex::Real T = S(i, j, k, UTEMP);
  amrex::Real e = S(i, j, k, UEINT) * rhoInv;
  amrex::Real massfrac[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = S(i, j, k, UFS + n) * rhoInv;
  }
  amrex::Real rho = S(i, j, k, URHO);
  auto eos = pele::physics::PhysicsType::eos();
  eos.REY2T(rho, e, massfrac, T);
  S(i, j, k, UTEMP) = T;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_rst_int_e(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& S,
  const bool allow_small_energy,
  const bool allow_negative_energy,
  const bool dual_energy_update_E_from_e,
  const amrex::Real dual_energy_eta2,
  const int /*verbose*/)
{
  if (std::is_same<
        pele::physics::PhysicsType::eos_type, pele::physics::eos::SRK>::value) {
    if ((!allow_small_energy) || (!allow_negative_energy)) {
      amrex::Error("Small/negative energy must be allowed for nonideal EOS");
    }
  }
  auto eos = pele::physics::PhysicsType::eos();
  if (!allow_small_energy) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke = 0.5 * (Up * Up + Vp * Vp + Wp * Wp);
    const amrex::Real eden = S(i, j, k, UEDEN) * rhoInv;
    const amrex::Real eos_state_T = constants::very_small_num();
    amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
    amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
    eos.T2Ei(eos_state_T, eos_state_ei);
    amrex::Real eos_state_e = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
      eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
    }
    const amrex::Real small_e = eos_state_e;
    if (eden < small_e) {
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc = amrex::max<amrex::Real>(
          S(i, j, k, UTEMP), std::numeric_limits<amrex::Real>::min());
        eos.T2Ei(eos_state_T_loc, eos_state_ei);
        eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
      S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
    } else {
      const amrex::Real rho_eint = S(i, j, k, UEDEN) - S(i, j, k, URHO) * ke;
      if (
        (rho_eint >
         (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) &&
        ((rho_eint -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) /
         (S(i, j, k, UEDEN) -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) >
          dual_energy_eta2) {
        S(i, j, k, UEINT) = rho_eint;
      }
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc = amrex::max<amrex::Real>(
          S(i, j, k, UTEMP), std::numeric_limits<amrex::Real>::min());
        eos.T2Ei(eos_state_T_loc, eos_state_ei);
        eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        if (dual_energy_update_E_from_e) {
          S(i, j, k, UEDEN) =
            S(i, j, k, UEDEN) +
            (S(i, j, k, URHO) * eos_state_e - S(i, j, k, UEINT));
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
    }
  } else if (!allow_negative_energy) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke = 0.5 * (Up * Up + Vp * Vp + Wp * Wp);
    if (
      S(i, j, k, UEDEN) <
      (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
      if (
        S(i, j, k, UEINT) <
        (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
        const amrex::Real eos_state_T = constants::very_small_num();
        amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
        amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
        eos.T2Ei(eos_state_T, eos_state_ei);
        amrex::Real eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
      S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
    } else {
      const amrex::Real rho_eint = S(i, j, k, UEDEN) - S(i, j, k, URHO) * ke;
      if (
        (rho_eint >
         (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) &&
        ((rho_eint -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) /
         (S(i, j, k, UEDEN) -
          S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) >
          dual_energy_eta2) {
        S(i, j, k, UEINT) = rho_eint;
      } else if (
        S(i, j, k, UEINT) >
          (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min()) &&
        (dual_energy_update_E_from_e)) {
        S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
      } else if (
        S(i, j, k, UEINT) <=
        (S(i, j, k, URHO) * std::numeric_limits<amrex::Real>::min())) {
        const amrex::Real eos_state_T = constants::very_small_num();
        amrex::Real eos_state_massfrac[NUM_SPECIES] = {0.0};
        amrex::Real eos_state_ei[NUM_SPECIES] = {0.0};
        eos.T2Ei(eos_state_T, eos_state_ei);
        amrex::Real eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        const amrex::Real eint_new = eos_state_e;
        if (dual_energy_update_E_from_e) {
          S(i, j, k, UEDEN) = S(i, j, k, UEDEN) +
                              (S(i, j, k, URHO) * eint_new - S(i, j, k, UEINT));
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eint_new;
      }
    }
  } else {
    const amrex::Real rho = S(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke = 0.5 * (Up * Up + Vp * Vp + Wp * Wp);
    S(i, j, k, UEINT) = S(i, j, k, UEDEN) - rho * ke;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& qa,
  PassMap const& pmap)
{
  auto eos = pele::physics::PhysicsType::eos();
  const amrex::Real rho = u(i, j, k, URHO);
  const amrex::Real rhoinv = 1.0 / rho;
  const amrex::Real vx = u(i, j, k, UMX) * rhoinv;
  const amrex::Real vy = u(i, j, k, UMY) * rhoinv;
  const amrex::Real vz = u(i, j, k, UMZ) * rhoinv;
  const amrex::Real kineng = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
  q(i, j, k, QRHO) = rho;
  q(i, j, k, QU) = vx;
  q(i, j, k, QV) = vy;
  q(i, j, k, QW) = vz;
  for (int ipassive = 0; ipassive < NPASSIVE; ++ipassive) {
    const int n = pmap.upassMap[ipassive];
    const int nq = pmap.qpassMap[ipassive];
    q(i, j, k, nq) = u(i, j, k, n) / rho;
  }

  const amrex::Real e = (u(i, j, k, UEDEN) - kineng) * rhoinv;
  amrex::Real T = u(i, j, k, UTEMP);
  amrex::Real massfrac[NUM_SPECIES];
  for (int sp = 0; sp < NUM_SPECIES; ++sp) {
    if (
      (-1e-4 * std::numeric_limits<amrex::Real>::epsilon() <
       q(i, j, k, sp + QFS)) &&
      (q(i, j, k, sp + QFS) <
       1e-4 * std::numeric_limits<amrex::Real>::epsilon())) {
      q(i, j, k, sp + QFS) = 0.0;
    }
    massfrac[sp] = q(i, j, k, sp + QFS);
  }

  amrex::Real dpdr_e, dpde, gam1, cs, wbar, p;

  eos.Y2WBAR(massfrac, wbar);
  eos.REY2T(rho, e, massfrac, T);
  eos.RTY2P(rho, T, massfrac, p);
  eos.RTY2Cs(rho, T, massfrac, cs);
  eos.RTY2G(rho, T, massfrac, gam1);
  eos.RTY2dpde_dpdre(rho, T, massfrac, dpde, dpdr_e);

  q(i, j, k, QTEMP) = T;
  q(i, j, k, QREINT) = e * rho;
  q(i, j, k, QPRES) = p;
  q(i, j, k, QGAME) = p / (e * rho) + 1.0;

  // Auxilary Fab
  qa(i, j, k, QDPDR) = dpdr_e;
  qa(i, j, k, QDPDE) = dpde;
  qa(i, j, k, QGAMC) = gam1;
  qa(i, j, k, QC) = cs;
  qa(i, j, k, QCSML) = amrex::max<amrex::Real>(
    constants::small_num(), constants::small_num() * cs);
  qa(i, j, k, QRSPEC) = pele::physics::Constants::RU / wbar;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_move_transcoefs_to_ec(
  const int i,
  const int j,
  const int k,
  const int n,
  const amrex::Array4<const amrex::Real>& carr,
  amrex::Real* earr,
  const int dir,
  const int do_harmonic)

{
  if (do_harmonic == 0) {
    if (dir == 0) {
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i - 1, j, k, n));
    } else if (dir == 1) {
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j - 1, k, n));
    } else {
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j, k - 1, n));
    }
  } else {
    amrex::Real a = carr(i, j, k, n);
    if (dir == 0) {
      amrex::Real b = carr(i - 1, j, k, n);
      if (a * b > 0.0) {
        earr[n] = 2.0 * (a * b) / (a + b);
      } else {
        earr[n] = 0.0;
      }
    } else if (dir == 1) {
      amrex::Real b = carr(i, j - 1, k, n);
      if (a * b > 0.0) {
        earr[n] = 2.0 * (a * b) / (a + b);
      } else {
        earr[n] = 0.0;
      }
    } else {
      amrex::Real b = carr(i, j, k - 1, n);
      if (a * b > 0.0) {
        earr[n] = 2.0 * (a * b) / (a + b);
      } else {
        earr[n] = 0.0;
      }
    }
  }
}

AMREX_FORCE_INLINE
void
setV(
  const amrex::Box& box,
  const int Ncomp,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      array(i, j, k, n) = val;
    });
}

AMREX_FORCE_INLINE
void
setC(
  const amrex::Box& box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++) {
      array(i, j, k, n) = val;
    }
  });
}

AMREX_FORCE_INLINE
void
copy_array4(
  const amrex::Box& box,
  const int Ncomp,
  const amrex::Array4<const amrex::Real>& in,
  const amrex::Array4<amrex::Real>& out)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      out(i, j, k, n) = in(i, j, k, n);
    });
}

AMREX_FORCE_INLINE
void
lincomb_array4(
  const amrex::Box& box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real alpha,
  const amrex::Real beta,
  const amrex::Array4<amrex::Real>& d)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++) {
      d(i, j, k, n) = alpha * f1(i, j, k, n) + beta * f2(i, j, k, n);
    }
  });
}

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<amrex::Real>& data);

void read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data);

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

template <typename T>
void
sort(T& vec)
{
  typename T::value_type* d_vec = vec.data();
  const int vec_size = vec.size();
  // Serial loop on the GPU
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) {
    for (int i = 0; i < vec_size - 1; i++) {
      for (int j = 0; j < vec_size - i - 1; j++) {
        if (d_vec[j + 1] < d_vec[j]) {
          // swap<typename T::value_type>(&d_vec[j], &d_vec[j + 1]);
          typename T::value_type temp = d_vec[j];
          d_vec[j] = d_vec[j + 1];
          d_vec[j + 1] = temp;
        }
      }
    }
  });
}

// Return vector of unique elements in input. Assume input is sorted. This is
// device safe code
template <typename T>
T
unique(T input)
{
  typename T::value_type* d_input = input.data();
  const int input_size = input.size();
  if (input_size == 0) {
    T output(input_size);
    return output;
  }

  // count the number of uniques
  const int Nunique = amrex::Reduce::Sum<int>(
    input.size() - 1,
    [=] AMREX_GPU_DEVICE(int i) noexcept -> int {
      return (d_input[i] != d_input[i + 1]);
    },
    1);

  // allocate the memory
  T output(Nunique);
  if (output.empty()) {
    return output;
  }
  typename T::value_type* d_output = output.data();

  // get the uniques
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*dummy*/) {
    int cnt = 1;
    d_output[0] = d_input[0];
    for (int i = 0; i < input_size - 1; i++) {
      if (d_input[i] != d_input[i + 1]) {
        d_output[cnt] = d_input[i + 1];
        cnt++;
      }
    }
  });

  return output;
}

// Find position of element in vector
template <typename T>
int
find_position(const std::vector<T>& vec, const T& element)
{
  auto it = std::find(vec.begin(), vec.end(), element);
  if (it != vec.end()) {
    return distance(vec.begin(), it);
  }
  return -1;
}

std::string convertIntGG(int number);

// Clean the mass fractions on state, given a mask
void clean_massfrac(
  const amrex::Box& /*bx*/,
  const amrex::Real /*threshold*/,
  amrex::Array4<const int> const& /*mask*/,
  amrex::Array4<amrex::Real> const& /*state*/);

#endif
