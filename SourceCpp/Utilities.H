#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <AMReX_FArrayBox.H>
#include "Constants.H"
#include "IndexDefines.H"
#include "EOS.H"

AMREX_GPU_DEVICE
void pc_cmpTemp(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& S);

AMREX_GPU_DEVICE
void pc_rst_int_e(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& S,
  const int allow_small_energy,
  const int allow_negative_energy,
  const int dual_energy_update_E_from_e,
  const int verbose);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_ctoprim(
  const int i,
  const int j,
  const int k,
  amrex::Array4<const amrex::Real> const& u,
  amrex::Array4<amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& qa)
{
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
    const int n = indxmap::upass_map(ipassive);
    const int nq = indxmap::qpass_map(ipassive);
    q(i, j, k, nq) = u(i, j, k, n) / rho;
  }

  const amrex::Real e = (u(i, j, k, UEDEN) - kineng) * rhoinv;
  amrex::Real T = u(i, j, k, UTEMP);
  amrex::Real massfrac[NUM_SPECIES];
  //    amrex::Real aux[NUM_AUX];
  for (int sp = 0; sp < NUM_SPECIES; ++sp)
    massfrac[sp] = q(i, j, k, sp + QFS);
  //    for(int ax = 0; ax < NUM_AUX; ++ax) aux[ax] = q(i,j,k,ax+QFX);
  amrex::Real dpdr_e, dpde, gam1, cs, wbar, p;

  // Are all these EOS calls needed? Seems fairly convoluted.
  EOS::Y2WBAR(massfrac, wbar);
  EOS::EY2T(e, massfrac, T);
  EOS::RTY2P(rho, T, massfrac, p);
  EOS::RTY2Cs(rho, T, massfrac, cs);
  EOS::TY2G(T, massfrac, gam1);
  EOS::RPE2dpdr_e(rho, p, e, dpdr_e);
  EOS::RG2dpde(rho, gam1, dpde);

  q(i, j, k, QTEMP) = T;
  q(i, j, k, QREINT) = e * rho;
  q(i, j, k, QPRES) = p;
  q(i, j, k, QGAME) = p / (e * rho) + 1.0;

  // Auxilary Fab
  qa(i, j, k, QDPDR) = dpdr_e;
  qa(i, j, k, QDPDE) = dpde;
  qa(i, j, k, QGAMC) = gam1;
  qa(i, j, k, QC) = cs;
  qa(i, j, k, QCSML) = amrex::max(SMALL, SMALL * cs);
  qa(i, j, k, QRSPEC) = EOS::RU / wbar;
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
    if (dir == 0)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i - 1, j, k, n));
    else if (dir == 1)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j - 1, k, n));
    else
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j, k - 1, n));
  } else {
    amrex::Real a = carr(i, j, k, n);
    if (dir == 0) {
      amrex::Real b = carr(i - 1, j, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else if (dir == 1) {
      amrex::Real b = carr(i, j - 1, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else {
      amrex::Real b = carr(i, j, k - 1, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    }
  }
}

AMREX_FORCE_INLINE
void
setV(
  const amrex::Box box,
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
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      array(i, j, k, n) = val;
  });
}

AMREX_FORCE_INLINE
void
copy_array4(
  const amrex::Box box,
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
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real alpha,
  const amrex::Real beta,
  const amrex::Array4<amrex::Real>& d)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      d(i, j, k, n) = alpha * f1(i, j, k, n) + beta * f2(i, j, k, n);
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
  const std::string iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<amrex::Real>& data);

void read_csv(
  const std::string iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data);

AMREX_GPU_HOST_DEVICE
void locate(const amrex::Real* xtable, const int n, amrex::Real& x, int& idxlo);

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
swap(T* x, T* y)
{
  T temp = *x;
  *x = *y;
  *y = temp;
}

template <typename T>
void
sort(T& vec)
{
  typename T::value_type* d_vec = vec.data();
  const int vec_size = vec.size();
  // Serial loop on the GPU
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int dummy) {
    for (int i = 0; i < vec_size - 1; i++) {
      for (int j = 0; j < vec_size - i - 1; j++) {
        if (d_vec[j + 1] < d_vec[j]) {
          swap<typename T::value_type>(&d_vec[j], &d_vec[j + 1]);
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
  int init_val = 0;
  amrex::Gpu::LaunchSafeGuard lsg(true);
  amrex::Gpu::DeviceScalar<int> ds(init_val);
  int* dp = ds.dataPtr();
  amrex::VecReduce(
    input.size() - 1, init_val,
    [=] AMREX_GPU_DEVICE(int i, int* r) noexcept {
      if (d_input[i] != d_input[i + 1]) {
        *r += 1;
      }
    },
    [=] AMREX_GPU_DEVICE(int const& r) noexcept {
      amrex::Gpu::deviceReduceSum(dp, r);
    });
  const int Nunique = ds.dataValue() + 1;

  // allocate the memory
  T output(Nunique);
  if (output.size() == 0) {
    return output;
  }
  typename T::value_type* d_output = output.data();

  // get the uniques
  amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int dummy) {
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

#endif
