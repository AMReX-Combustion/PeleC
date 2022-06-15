#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>
#include <memory>

#ifdef PELEC_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "Transport.H"
#include "mechanism.H"
#include "PeleC.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"

ProbParmDevice* PeleC::d_prob_parm_device = nullptr;
ProbParmDevice* PeleC::h_prob_parm_device = nullptr;
ProbParmHost* PeleC::prob_parm_host = nullptr;
TaggingParm* PeleC::tagging_parm = nullptr;

// Components are:
// Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall, UserBC
static int scalar_bc[] = {INT_DIR,      EXT_DIR,      FOEXTRAP, REFLECT_EVEN,
                          REFLECT_EVEN, REFLECT_EVEN, EXT_DIR};

static int norm_vel_bc[] = {INT_DIR,     EXT_DIR,     FOEXTRAP, REFLECT_ODD,
                            REFLECT_ODD, REFLECT_ODD, EXT_DIR};

static int tang_vel_bc[] = {INT_DIR,      EXT_DIR,     FOEXTRAP, REFLECT_EVEN,
                            REFLECT_EVEN, REFLECT_ODD, EXT_DIR};

static int react_src_bc[] = {INT_DIR,      REFLECT_EVEN, REFLECT_EVEN,
                             REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN,
                             REFLECT_EVEN};

static void
set_scalar_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, scalar_bc[lo_bc[dir]]);
    bc.setHi(dir, scalar_bc[hi_bc[dir]]);
  }
}

static void
set_x_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, norm_vel_bc[lo_bc[0]]); bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_y_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, norm_vel_bc[lo_bc[1]]); bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_react_src_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, react_src_bc[lo_bc[dir]]);
    bc.setHi(dir, react_src_bc[hi_bc[dir]]);
  }
}

static void
set_z_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, norm_vel_bc[lo_bc[2]]); bc.setHi(2, norm_vel_bc[hi_bc[2]]););
}

void
PeleC::variableSetUp()
{
  // PeleC::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job

  // initialize the start time for our CPU-time tracker
  startCPUTime = amrex::ParallelDescriptor::second();

  // Output the git commit hashes used to build the executable.

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const char* pelec_hash = amrex::buildInfoGetGitHash(1);
    const char* amrex_hash = amrex::buildInfoGetGitHash(2);
    const char* pelephysics_hash = amrex::buildInfoGetGitHash(3);
    const char* amrexhydro_hash = amrex::buildInfoGetGitHash(4);
    const char* sundials_hash = amrex::buildInfoGetGitHash(5);
    const char* buildgithash = amrex::buildInfoGetBuildGitHash();
    const char* buildgitname = amrex::buildInfoGetBuildGitName();

    if (strlen(pelec_hash) > 0) {
      amrex::Print() << "\n"
                     << "PeleC git hash: " << pelec_hash << "\n";
    }
    if (strlen(amrex_hash) > 0) {
      amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
    }
    if (strlen(pelephysics_hash) > 0) {
      amrex::Print() << "PelePhysics git hash: " << pelephysics_hash << "\n";
    }
    if (strlen(amrexhydro_hash) > 0) {
      amrex::Print() << "AMReX-Hydro git hash: " << amrexhydro_hash << "\n";
    }
    if (strlen(sundials_hash) > 0) {
      amrex::Print() << "SUNDIALS git hash: " << sundials_hash << "\n";
    }
    if (strlen(buildgithash) > 0) {
      amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
    }

    amrex::Print() << "\n";
  }

  AMREX_ASSERT(desc_lst.size() == 0);

  prob_parm_host = new ProbParmHost{};
  h_prob_parm_device = new ProbParmDevice{};
  tagging_parm = new TaggingParm{};
  d_prob_parm_device = static_cast<ProbParmDevice*>(
    amrex::The_Arena()->alloc(sizeof(ProbParmDevice)));
  trans_parms.allocate();
  turb_inflow.init(amrex::DefaultGeometry());

  // Get options, set phys_bc
  eb_in_domain = ebInDomain();
  read_params();

#ifdef PELEC_USE_MASA
  if (do_mms) {
    init_mms();
  }
#endif

  // Set number of state variables and pointers to components

  int cnt = 0;
  Density = cnt++;
  Xmom = cnt++;
  Ymom = cnt++;
  Zmom = cnt++;
  Eden = cnt++;
  Eint = cnt++;
  Temp = cnt++;

  if (NUM_ADV > 0) {
    FirstAdv = cnt;
    cnt += NUM_ADV;
  }

  // int dm = AMREX_SPACEDIM;

  if (NUM_SPECIES > 0) {
    FirstSpec = cnt;
    cnt += NUM_SPECIES; // NOLINT
  }

  if (NUM_AUX > 0) {
    FirstAux = cnt;
    cnt += NUM_AUX;
  }

  if (NUM_LIN > 0) {
    FirstLin = cnt;
    cnt += NUM_LIN;
  }

  // NUM_LIN variables are will be added by the specific models
  // NVAR = cnt;

  // const amrex::Real run_strt = amrex::ParallelDescriptor::second() ;
  // Real run_stop = ParallelDescriptor::second() - run_strt;
  // ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());

  // if (ParallelDescriptor::IOProcessor())
  //    amrex::Print() << "\nTime in set_method_params: " << run_stop << '\n'
  //    ;

  // if (nscbc_adv == 1 && amrex::ParallelDescriptor::IOProcessor()) {
  //  amrex::Print() << "Using Ghost-Cells Navier-Stokes Characteristic BCs
  //  for
  //  "
  //                    "advection: nscbc_adv = "
  //                 << nscbc_adv << '\n'
  //                 << '\n';
  //}

  // if (nscbc_diff == 1 && amrex::ParallelDescriptor::IOProcessor()) {
  //  amrex::Print() << "Using Ghost-Cells Navier-Stokes Characteristic BCs
  //  for
  //  "
  //                    "diffusion: nscbc_diff = "
  //                 << nscbc_diff << '\n'
  //                 << '\n';
  //}

  // int coord_type = amrex::DefaultGeometry().Coord();

  amrex::MFInterpolater* interp;

  if (state_interp_order == 0) {
    interp = &amrex::mf_pc_interp;
  } else {
    if (lin_limit_state_interp == 1) {
      interp = &amrex::mf_lincc_interp;
    } else {
      interp = &amrex::mf_cell_cons_interp;
    }
  }

  // Decide from input whether eb interp is needed for FillPatch operations
  if (eb_in_domain && (state_interp_order != 0)) {
    if (lin_limit_state_interp == 1) {
      interp = &amrex::eb_mf_lincc_interp;
    } else {
      interp = &amrex::eb_mf_cell_cons_interp;
    }
  }

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint = true;

  int ngrow_state = state_nghost;
  AMREX_ASSERT(ngrow_state >= 0);

  desc_lst.addDescriptor(
    State_Type, amrex::IndexType::TheCellType(), amrex::StateDescriptor::Point,
    ngrow_state, NVAR, interp, state_data_extrap, store_in_checkpoint);

  // Components 0:Numspec-1 are rho.omega_i
  // Component NUM_SPECIES is rho.edot = (rho.eout-rho.ein)
  store_in_checkpoint = do_react;
  desc_lst.addDescriptor(
    Reactions_Type, amrex::IndexType::TheCellType(),
    amrex::StateDescriptor::Point, 0, NUM_SPECIES + 2, interp,
    state_data_extrap, store_in_checkpoint);
  amrex::Vector<amrex::BCRec> bcs(NVAR);
  amrex::Vector<std::string> name(NVAR);
  amrex::Vector<amrex::BCRec> react_bcs(NUM_SPECIES + 2);
  amrex::Vector<std::string> react_name(NUM_SPECIES + 2);

  amrex::BCRec bc;
  cnt = 0;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "density";
  cnt++;
  set_x_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "xmom";
  cnt++;
  set_y_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "ymom";
  cnt++;
  set_z_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "zmom";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "rho_E";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "rho_e";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "Temp";

  for (int i = 0; i < NUM_ADV; ++i) {
    char buf[64];
    sprintf(buf, "rho_adv_%d", i);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = std::string(buf);
  }

  // Get the species names from the network model
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << NUM_SPECIES << " Species: " << std::endl;
    for (int i = 0; i < NUM_SPECIES; i++) {
      amrex::Print() << spec_names[i] << ' ' << ' ';
    }
    amrex::Print() << std::endl;
  }

  for (int i = 0; i < NUM_SPECIES; ++i) {
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "rho_" + spec_names[i];
  }
  // Get the auxiliary names from the network model.
  amrex::Vector<std::string> aux_names;
  for (int i = 0; i < NUM_AUX; i++) {
    int len = 20;
    amrex::Vector<int> int_aux_names(len);

    // Disabling for the GPU at the moment. Look at the species names to see
    // how to do this in C++. AUX stuff is usually 0 anyway. This call returns
    // the actual length of each string in "len"
    // get_aux_names(int_aux_names.dataPtr(),&i,&len);

    char* char_aux_names = new char[len + 1];
    for (int j = 0; j < len; j++) {
      char_aux_names[j] = int_aux_names[j];
    }
    char_aux_names[len] = '\0';
    aux_names.push_back(std::string(char_aux_names));
    delete[] char_aux_names;
  }
  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << NUM_AUX << " Auxiliary Variables: " << std::endl;
    for (int i = 0; i < NUM_AUX; i++) {
      amrex::Print() << aux_names[i] << ' ' << ' ';
    }
    amrex::Print() << std::endl;
  }

  for (int i = 0; i < NUM_AUX; ++i) {
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "rho_" + aux_names[i];
  }

  amrex::StateDescriptor::BndryFunc bndryfunc1(pc_bcfill_hyp);
  bndryfunc1.setRunOnGPU(true);

  desc_lst.setComponent(State_Type, Density, name, bcs, bndryfunc1);

  for (int i = 0; i < NUM_SPECIES; ++i) {
    set_react_src_bc(bc, phys_bc);
    react_bcs[i] = bc;
    react_name[i] = "rho_omega_" + spec_names[i];
  }
  set_react_src_bc(bc, phys_bc);
  react_bcs[NUM_SPECIES] = bc;
  react_name[NUM_SPECIES] = "rhoe_dot";
  react_bcs[NUM_SPECIES + 1] = bc;
  react_name[NUM_SPECIES + 1] = "heatRelease";

  amrex::StateDescriptor::BndryFunc bndryfunc2(pc_reactfill_hyp);
  bndryfunc2.setRunOnGPU(true);

  desc_lst.setComponent(Reactions_Type, 0, react_name, react_bcs, bndryfunc2);

  const bool workest_store_in_checkpoint = false;
  const bool workest_data_extrap = false;
  desc_lst.addDescriptor(
    Work_Estimate_Type, amrex::IndexType::TheCellType(),
    amrex::StateDescriptor::Point, 0, 1, &amrex::pc_interp, workest_data_extrap,
    workest_store_in_checkpoint);
  // Because we use piecewise constant interpolation, we do not use bc and
  // BndryFunc.
  desc_lst.setComponent(
    Work_Estimate_Type, 0, "WorkEstimate", bc,
    amrex::StateDescriptor::BndryFunc(pc_nullfill));

  num_state_type = desc_lst.size();

  // DEFINE DERIVED QUANTITIES

  // Pressure
  derive_lst.add(
    "pressure", amrex::IndexType::TheCellType(), 1, pc_derpres,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("pressure", desc_lst, State_Type, Density, NVAR);

  // Kinetic energy
  derive_lst.add(
    "kineng", amrex::IndexType::TheCellType(), 1, pc_derkineng,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("kineng", desc_lst, State_Type, Density, NVAR);

  // Enstrophy
  derive_lst.add(
    "enstrophy", amrex::IndexType::TheCellType(), 1, pc_derenstrophy,
    amrex::DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("enstrophy", desc_lst, State_Type, Density, NVAR);

  // Sound speed (c)
  derive_lst.add(
    "soundspeed", amrex::IndexType::TheCellType(), 1, pc_dersoundspeed,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density, NVAR);

  // Mach number(M)
  derive_lst.add(
    "MachNumber", amrex::IndexType::TheCellType(), 1, pc_dermachnumber,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density, NVAR);

  // Entropy (S)
  derive_lst.add(
    "entropy", amrex::IndexType::TheCellType(), 1, pc_derentropy,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("entropy", desc_lst, State_Type, Density, NVAR);

  // Vorticity
  derive_lst.add(
    "magvort", amrex::IndexType::TheCellType(), 1, pc_dermagvort,
    amrex::DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("magvort", desc_lst, State_Type, Density, NVAR);

  // Div(u)
  derive_lst.add(
    "divu", amrex::IndexType::TheCellType(), 1, pc_derdivu,
    amrex::DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("divu", desc_lst, State_Type, Density, NVAR);

  // Internal energy as derived from rho*E, part of the state
  derive_lst.add(
    "eint_E", amrex::IndexType::TheCellType(), 1, pc_dereint1,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("eint_E", desc_lst, State_Type, Density, NVAR);

  // Internal energy as derived from rho*e, part of the state
  derive_lst.add(
    "eint_e", amrex::IndexType::TheCellType(), 1, pc_dereint2,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("eint_e", desc_lst, State_Type, Density, NVAR);

  // Log(density)
  derive_lst.add(
    "logden", amrex::IndexType::TheCellType(), 1, pc_derlogden,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("logden", desc_lst, State_Type, Density, NVAR);

  // Y from rhoY
  amrex::Vector<std::string> var_names_massfrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    var_names_massfrac[i] = "Y(" + spec_names[i] + ")";
  }

  derive_lst.add(
    "massfrac", amrex::IndexType::TheCellType(), NUM_SPECIES,
    var_names_massfrac, pc_derspec, amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("massfrac", desc_lst, State_Type, Density, NVAR);

  // adv from rho_adv
  amrex::Vector<std::string> var_names_adv(NUM_ADV);
  for (int i = 0; i < NUM_ADV; i++) {
    var_names_adv[i] = "adv_" + std::to_string(i);
  }

  derive_lst.add(
    "adv", amrex::IndexType::TheCellType(), NUM_ADV, var_names_adv, pc_deradv,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("adv", desc_lst, State_Type, Density, NVAR);

  // Species mole fractions
  amrex::Vector<std::string> var_names_molefrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    var_names_molefrac[i] = "X(" + spec_names[i] + ")";
  }

  derive_lst.add(
    "molefrac", amrex::IndexType::TheCellType(), NUM_SPECIES,
    var_names_molefrac, pc_dermolefrac, amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("molefrac", desc_lst, State_Type, Density, NVAR);

  // Velocities
  derive_lst.add(
    "x_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelx,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "y_velocity", amrex::IndexType::TheCellType(), 1, pc_dervely,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "z_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelz,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "magvel", amrex::IndexType::TheCellType(), 1, pc_dermagvel,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("magvel", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "radvel", amrex::IndexType::TheCellType(), 1, pc_derradialvel,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("radvel", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "magmom", amrex::IndexType::TheCellType(), 1, pc_dermagmom,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("magmom", desc_lst, State_Type, Density, NVAR);

#ifdef AMREX_PARTICLES
  // We want a derived type that corresponds to the number of particles
  // in each cell.  We only intend to use it in plotfiles for debugging
  // purposes. We'll actually set the values in writePlotFile().
  derive_lst.add(
    "particle_count", amrex::IndexType::TheCellType(), 1, pc_dernull,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("particle_count", desc_lst, State_Type, Density, 1);

  derive_lst.add(
    "total_particle_count", amrex::IndexType::TheCellType(), 1, pc_dernull,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent(
    "total_particle_count", desc_lst, State_Type, Density, 1);

  derive_lst.add(
    "particle_density", amrex::IndexType::TheCellType(), 1, pc_dernull,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("particle_density", desc_lst, State_Type, Density, 1);
#endif

  derive_lst.add(
    "cp", amrex::IndexType::TheCellType(), 1, pc_dercp,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("cp", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "cv", amrex::IndexType::TheCellType(), 1, pc_dercv,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("cv", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "viscosity", amrex::IndexType::TheCellType(), 1, pc_derviscosity,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("viscosity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "bulk_viscosity", amrex::IndexType::TheCellType(), 1, pc_derbulkviscosity,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent(
    "bulk_viscosity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "conductivity", amrex::IndexType::TheCellType(), 1, pc_derconductivity,
    amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("conductivity", desc_lst, State_Type, Density, NVAR);

  amrex::Vector<std::string> var_names_diffusivity(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    var_names_diffusivity[i] = "D(" + spec_names[i] + ")";
  }
  derive_lst.add(
    "diffusivity", amrex::IndexType::TheCellType(), NUM_SPECIES,
    var_names_diffusivity, pc_derdiffusivity, amrex::DeriveRec::TheSameBox);
  derive_lst.addComponent("diffusivity", desc_lst, State_Type, Density, NVAR);

  // LES coefficients
  if (do_les) {
    derive_lst.add(
      "C_s2", amrex::IndexType::TheCellType(), 1, pc_dernull,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("C_s2", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "C_I", amrex::IndexType::TheCellType(), 1, pc_dernull,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("C_I", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "Pr_T", amrex::IndexType::TheCellType(), 1, pc_dernull,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("Pr_T", desc_lst, State_Type, Density, 1);
  }

  // MMS derives
#ifdef PELEC_USE_MASA
  if (do_mms) {
    derive_lst.add(
      "rhommserror", amrex::IndexType::TheCellType(), 1, pc_derrhommserror,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("rhommserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "ummserror", amrex::IndexType::TheCellType(), 1, pc_derummserror,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("ummserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "vmmserror", amrex::IndexType::TheCellType(), 1, pc_dervmmserror,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("vmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "wmmserror", amrex::IndexType::TheCellType(), 1, pc_derwmmserror,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("wmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "pmmserror", amrex::IndexType::TheCellType(), 1, pc_derpmmserror,
      amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("pmmserror", desc_lst, State_Type, Density, NVAR);
  }
#endif

  // Problem-specific derives
  add_problem_derives<ProblemDerives>(derive_lst, desc_lst);

  // Set list of active sources
  set_active_sources();
#ifdef AMREX_PARTICLES
  defineParticles();
#endif
}

void
PeleC::variableCleanUp()
{
  derive_lst.clear();

  desc_lst.clear();

  clear_prob();

  eb_initialized = false;

  delete prob_parm_host;
  delete tagging_parm;
  delete h_prob_parm_device;
  amrex::The_Arena()->free(d_prob_parm_device);
  trans_parms.deallocate();
}

void
PeleC::set_active_sources()
{
  if (do_diffuse && !do_mol) {
    src_list.push_back(diff_src);
  }

  // optional external source
  if (add_ext_src) {
    src_list.push_back(ext_src);
  }

  // optional forcing source
  if (add_forcing_src) {
    src_list.push_back(forcing_src);
  }

  if (do_spray_particles) {
    src_list.push_back(spray_src);
  }

  // optional LES source
  if (do_les) {
    src_list.push_back(les_src);
  }

#ifdef PELEC_USE_MASA
  // optional MMS source
  if (do_mms) {
    src_list.push_back(mms_src);
  }
#endif
}
