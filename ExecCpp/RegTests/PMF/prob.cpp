#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real pamb = 1013250.0 * 100.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real phi_in = -0.2;
AMREX_GPU_DEVICE_MANAGED amrex::Real T_in = 298.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real vn_in = 0.2;
AMREX_GPU_DEVICE_MANAGED amrex::Real pertmag = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {1.0};
AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N = 0;
AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M = 0;
AMREX_GPU_DEVICE_MANAGED bool pmf_do_average = false;

amrex::Gpu::ManagedVector<amrex::Real>* pmf_X = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* pmf_Y = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* fuel_state = nullptr;

AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_X = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_Y = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_fuel_state = nullptr;

std::string pmf_datafile = "";
amrex::Vector<std::string> pmf_names;
} // namespace ProbParm

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"')
      count++;
  }
  if ((count % 2) == 0)
    return true;
  else
    return false;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline, secondline, remaininglines;
  unsigned int pos1, pos2;
  int variable_count, line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline))
    amrex::Abort("PMF file variable quotes unbalanced");
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  ProbParm::pmf_names.resize(variable_count);
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    ProbParm::pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << ProbParm::pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  ProbParm::pmf_N = line_count;
  ProbParm::pmf_M = variable_count - 1;
  ProbParm::pmf_X->resize(ProbParm::pmf_N);
  ProbParm::pmf_Y->resize(ProbParm::pmf_N * ProbParm::pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < ProbParm::pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> (*ProbParm::pmf_X)[i];
    for (unsigned int j = 0; j < ProbParm::pmf_M; j++) {
      sinput >> (*ProbParm::pmf_Y)[j * ProbParm::pmf_N + i];
    }
  }
  ProbParm::d_pmf_X = ProbParm::pmf_X->dataPtr();
  ProbParm::d_pmf_Y = ProbParm::pmf_Y->dataPtr();
}

void
init_bc()
{
  amrex::Real vt, ek, T, rho, e;
  amrex::Real molefrac[NUM_SPECIES], massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};

  if (ProbParm::phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    pmf(yl, yr, pmf_vals);
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    ProbParm::vn_in = pmf_vals[1];
  } else {
    const amrex::Real a = 0.5;
    for (int n = 0; n < NUM_SPECIES; n++)
      molefrac[n] = 0.0;
    molefrac[O2_ID] = 1.0 / (1.0 + ProbParm::phi_in / a + 0.79 / 0.21);
    molefrac[H2_ID] = ProbParm::phi_in * molefrac[O2_ID] / a;
    molefrac[N2_ID] = 1.0 - molefrac[H2_ID] - molefrac[O2_ID];
    T = ProbParm::T_in;
  }
  const int p = ProbParm::pamb;

  EOS::X2Y(molefrac, massfrac);
  EOS::PYT2RE(p, massfrac, T, rho, e);

  vt = ProbParm::vn_in;
  ek = 0.5 * (vt * vt);

  (*ProbParm::fuel_state)[URHO] = rho;
  (*ProbParm::fuel_state)[UMX] = 0.0;
  (*ProbParm::fuel_state)[UMY] = rho * vt;
  (*ProbParm::fuel_state)[UMZ] = 0.0;
  (*ProbParm::fuel_state)[UEINT] = rho * e;
  (*ProbParm::fuel_state)[UEDEN] = rho * (e + ek);
  (*ProbParm::fuel_state)[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    (*ProbParm::fuel_state)[UFS + n - 1] = rho * massfrac[n];
  }
  ProbParm::d_fuel_state = ProbParm::fuel_state->dataPtr();
}

void
pc_prob_close()
{
  delete ProbParm::pmf_X;
  delete ProbParm::pmf_Y;
  delete ProbParm::fuel_state;

  ProbParm::pmf_X = nullptr;
  ProbParm::pmf_Y = nullptr;
  ProbParm::fuel_state = nullptr;
  ProbParm::d_pmf_X = nullptr;
  ProbParm::d_pmf_Y = nullptr;
  ProbParm::d_fuel_state = nullptr;
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
  amrex::ParmParse pp("prob");
  pp.query("pamb", ProbParm::pamb);
  pp.query("phi_in", ProbParm::phi_in);
  pp.query("T_in", ProbParm::T_in);
  pp.query("vn_in", ProbParm::vn_in);
  pp.query("pertmag", ProbParm::pertmag);
  pp.query("pmf_datafile", ProbParm::pmf_datafile);

  ProbParm::L[0] = probhi[0] - problo[0];
  ProbParm::L[1] = probhi[1] - problo[1];
  ProbParm::L[2] = probhi[2] - problo[2];

  ProbParm::pmf_X = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::pmf_Y = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::fuel_state = new amrex::Gpu::ManagedVector<amrex::Real>;
  ProbParm::fuel_state->resize(NVAR);

  read_pmf(ProbParm::pmf_datafile);

  init_bc();
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
