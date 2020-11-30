#include "Utilities.H"

AMREX_GPU_DEVICE
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
  EOS::EY2T(e, massfrac, T);
  S(i, j, k, UTEMP) = T;
}

AMREX_GPU_DEVICE
void
pc_rst_int_e(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::Real> const& S,
  const int allow_small_energy,
  const int allow_negative_energy,
  const int dual_energy_update_E_from_e,
  const int /*verbose*/)
{
  if (allow_small_energy == 0) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke = 0.5 * (Up * Up + Vp * Vp + Wp * Wp);
    const amrex::Real eden = S(i, j, k, UEDEN) * rhoInv;
    // const amrex::Real eos_state_rho = S(i, j, k, URHO);
    const amrex::Real eos_state_T = SMALL_TEMP;
    amrex::Real eos_state_massfrac[NUM_SPECIES];
    amrex::Real eos_state_ei[NUM_SPECIES];
    EOS::T2Ei(eos_state_T, eos_state_ei);
    amrex::Real eos_state_e = 0.0;
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
      eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
    }
    const amrex::Real small_e = eos_state_e;
    if (eden < small_e) {
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc =
          amrex::max(S(i, j, k, UTEMP), SMALL_TEMP);
        EOS::T2Ei(eos_state_T_loc, eos_state_ei);
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
        (rho_eint > (S(i, j, k, URHO) * SMALL_E)) and
        ((rho_eint - S(i, j, k, URHO) * SMALL_E) /
         (S(i, j, k, UEDEN) - S(i, j, k, URHO) * SMALL_E)) > DUAL_ENERGY_ETA2) {
        S(i, j, k, UEINT) = rho_eint;
      }
      if (S(i, j, k, UEINT) * rhoInv < small_e) {
        const amrex::Real eos_state_T_loc =
          amrex::max(S(i, j, k, UTEMP), SMALL_TEMP);
        EOS::T2Ei(eos_state_T_loc, eos_state_ei);
        eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        if (dual_energy_update_E_from_e == 1) {
          S(i, j, k, UEDEN) =
            S(i, j, k, UEDEN) +
            (S(i, j, k, URHO) * eos_state_e - S(i, j, k, UEINT));
        }
        S(i, j, k, UEINT) = S(i, j, k, URHO) * eos_state_e;
      }
    }
  } else if (allow_negative_energy == 0) {
    const amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
    const amrex::Real Up = S(i, j, k, UMX) * rhoInv;
    const amrex::Real Vp = S(i, j, k, UMY) * rhoInv;
    const amrex::Real Wp = S(i, j, k, UMZ) * rhoInv;
    const amrex::Real ke = 0.5 * (Up * Up + Vp * Vp + Wp * Wp);
    if (S(i, j, k, UEDEN) < (S(i, j, k, URHO) * SMALL_E)) {
      if (S(i, j, k, UEINT) < (S(i, j, k, URHO) * SMALL_E)) {
        // const amrex::Real eos_state_rho = S(i, j, k, URHO);
        const amrex::Real eos_state_T = SMALL_TEMP;
        amrex::Real eos_state_massfrac[NUM_SPECIES];
        amrex::Real eos_state_ei[NUM_SPECIES];
        EOS::T2Ei(eos_state_T, eos_state_ei);
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
        (rho_eint > (S(i, j, k, URHO) * SMALL_E)) and
        ((rho_eint - S(i, j, k, URHO) * SMALL_E) /
         (S(i, j, k, UEDEN) - S(i, j, k, URHO) * SMALL_E)) > DUAL_ENERGY_ETA2) {
        S(i, j, k, UEINT) = rho_eint;
      } else if (
        S(i, j, k, UEINT) > (S(i, j, k, URHO) * SMALL_E) and
        (dual_energy_update_E_from_e == 1)) {
        S(i, j, k, UEDEN) = S(i, j, k, UEINT) + S(i, j, k, URHO) * ke;
      } else if (S(i, j, k, UEINT) <= (S(i, j, k, URHO) * SMALL_E)) {
        // const amrex::Real eos_state_rho = S(i, j, k, URHO);
        const amrex::Real eos_state_T = SMALL_TEMP;
        amrex::Real eos_state_massfrac[NUM_SPECIES];
        amrex::Real eos_state_ei[NUM_SPECIES];
        EOS::T2Ei(eos_state_T, eos_state_ei);
        amrex::Real eos_state_e = 0.0;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          eos_state_massfrac[sp] = S(i, j, k, sp + UFS) * rhoInv;
          eos_state_e += eos_state_massfrac[sp] * eos_state_ei[sp];
        }
        const amrex::Real eint_new = eos_state_e;
        /* Avoiding this for now due to GPU restricitions
        if (verbose) {
          amrex::Print() << std::endl;
          amrex::Print()
            << ">>> Warning: PeleC_util.F90::reset_internal_energy "
            << " " << i << ", " << j << ", " << k << std::endl;
          amrex::Print() << ">>> ... resetting neg. e from EOS using SMALL_TEMP"
                         << std::endl;
          amrex::Print() << ">>> ... from ',S(i,j,k,UEINT)/S(i,j,k,URHO),' to "
                         << std::endl;
          amrex::Print() << eint_new << std::endl;
        }*/
        if (dual_energy_update_E_from_e == 1) {
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

// -----------------------------------------------------------
// Read a binary file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<double>& data /*needs to be double*/)
{

  std::ifstream infile(iname, std::ios::in | std::ios::binary);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }

  for (size_t i = 0; i < nx * ny * nz * ncol; i++) {
    infile.read(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
  }
  infile.close();
}

// -----------------------------------------------------------
// Read a csv file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  size_t nlines = 0;
  std::string firstline, line;
  std::getline(iss, firstline); // skip header
  while (getline(iss, line))
    ++nlines;

  // Quick sanity check
  if (nlines != nx * ny * nz)
    amrex::Abort(
      "Number of lines in the input file (= " + std::to_string(nlines) +
      ") does not match the input resolution (=" + std::to_string(nx) + ")");

  // Read the data from the file
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline); // skip header
  int cnt = 0;
  while (std::getline(iss, line)) {
    std::istringstream linestream(line);
    std::string value;
    while (getline(linestream, value, ',')) {
      std::istringstream sinput(value);
      sinput >> data[cnt];
      cnt++;
    }
  }
}

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_HOST_DEVICE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  } else if (x <= xtable[0]) {
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
