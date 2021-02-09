#include "prob.H"

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
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
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

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::prob_parm_device->pmf_N = line_count;
  PeleC::prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    PeleC::prob_parm_device->pmf_N * PeleC::prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    PeleC::prob_parm_device->pmf_N * PeleC::prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < PeleC::prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (unsigned int j = 0; j < PeleC::prob_parm_device->pmf_M; j++) {
      sinput >>
        PeleC::prob_parm_host->h_pmf_Y[j * PeleC::prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real molefrac[NUM_SPECIES];
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};

  if (PeleC::prob_parm_device->phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    pmf(yl, yr, pmf_vals, *PeleC::prob_parm_device);
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    PeleC::prob_parm_device->vn_in = pmf_vals[1];
  } else {
    const amrex::Real a = 0.5;
    for (amrex::Real& n : molefrac) {
      n = 0.0;
    }
    molefrac[O2_ID] =
      1.0 / (1.0 + PeleC::prob_parm_device->phi_in / a + 0.79 / 0.21);
    molefrac[H2_ID] = PeleC::prob_parm_device->phi_in * molefrac[O2_ID] / a;
    molefrac[N2_ID] = 1.0 - molefrac[H2_ID] - molefrac[O2_ID];
    T = PeleC::prob_parm_device->T_in;
  }
  const amrex::Real p = PeleC::prob_parm_device->pamb;

  EOS::X2Y(molefrac, massfrac);
  EOS::PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::prob_parm_device->fuel_state[URHO] = rho;
  PeleC::prob_parm_device->fuel_state[UMX] = 0.0;
  PeleC::prob_parm_device->fuel_state[UMY] = rho * vt;
  PeleC::prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::prob_parm_device->fuel_state[UFS + n - 1] = rho * massfrac[n];
  }
}

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
  std::string pmf_datafile;

  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::prob_parm_device->pamb);
  pp.query("phi_in", PeleC::prob_parm_device->phi_in);
  pp.query("T_in", PeleC::prob_parm_device->T_in);
  pp.query("vn_in", PeleC::prob_parm_device->vn_in);
  pp.query("pertmag", PeleC::prob_parm_device->pertmag);
  pp.query("pmf_datafile", pmf_datafile);

  PeleC::prob_parm_device->L[0] = probhi[0] - problo[0];
  PeleC::prob_parm_device->L[1] = probhi[1] - problo[1];
  PeleC::prob_parm_device->L[2] = probhi[2] - problo[2];

  read_pmf(pmf_datafile);

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
