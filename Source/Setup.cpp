#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>
#include <memory>

#ifdef PELEC_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#if defined(PELEC_USE_REACTIONS) && defined(AMREX_USE_GPU) && \
  defined(USE_SUNDIALS_PP)
#include <AMReX_SUNMemory.H>
#endif

#include "Transport.H"
#include "mechanism.h"
#include "PeleC.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"
#include "chemistry_file.H"

ProbParmDevice* PeleC::d_prob_parm_device = nullptr;
ProbParmDevice* PeleC::h_prob_parm_device = nullptr;
ProbParmHost* PeleC::prob_parm_host = nullptr;
TaggingParm* PeleC::tagging_parm = nullptr;
PassMap* PeleC::d_pass_map = nullptr;
PassMap* PeleC::h_pass_map = nullptr;

// Components are:
// Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall, UserBC
static int scalar_bc[] = {INT_DIR,      EXT_DIR,      FOEXTRAP, REFLECT_EVEN,
                          REFLECT_EVEN, REFLECT_EVEN, EXT_DIR};

static int norm_vel_bc[] = {INT_DIR,     EXT_DIR,     FOEXTRAP, REFLECT_ODD,
                            REFLECT_ODD, REFLECT_ODD, EXT_DIR};

static int tang_vel_bc[] = {INT_DIR,      EXT_DIR,     FOEXTRAP, REFLECT_EVEN,
                            REFLECT_EVEN, REFLECT_ODD, EXT_DIR};

#ifdef PELEC_USE_REACTIONS
static int react_src_bc[] = {INT_DIR,      REFLECT_EVEN, REFLECT_EVEN,
                             REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN,
                             REFLECT_EVEN};
#endif

static amrex::Box
the_same_box(const amrex::Box& b)
{
  return b;
}
static amrex::Box
grow_box_by_one(const amrex::Box& b)
{
  return amrex::grow(b, 1);
}

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

#ifdef PELEC_USE_REACTIONS
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
#endif

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
    if (strlen(buildgithash) > 0) {
      amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
    }

    amrex::Print() << "\n";
  }

  AMREX_ASSERT(desc_lst.size() == 0);

  prob_parm_host = new ProbParmHost{};
  h_prob_parm_device = new ProbParmDevice{};
  tagging_parm = new TaggingParm{};
  h_pass_map = new PassMap{};
  d_prob_parm_device = static_cast<ProbParmDevice*>(
    amrex::The_Arena()->alloc(sizeof(ProbParmDevice)));
  d_pass_map =
    static_cast<PassMap*>(amrex::The_Arena()->alloc(sizeof(PassMap)));

  // Get options, set phys_bc
  read_params();

  pele::physics::transport::InitTransport<
    pele::physics::PhysicsType::eos_type>()();

#ifdef PELEC_USE_REACTIONS
#if defined(AMREX_USE_GPU) && defined(USE_SUNDIALS_PP)
  amrex::sundials::MemoryHelper::Initialize();
#endif

  if (chem_integrator == 1) {
    amrex::Print() << "Using built-in RK64 chemistry integrator\n";
  }
#ifdef USE_SUNDIALS_PP
  else if (chem_integrator == 2) {
    amrex::Print()
      << "Using sundials chemistry integrator with flattened arrays\n";
  } else if (chem_integrator == 3) {
    amrex::Print() << "Using sundials chemistry integrator with boxes\n";
  }
#endif
  else {
    amrex::Abort("Invalid chem_integrator choice.");
  }

  // Initialize the reactor
  if (do_react == 1) {
    init_reactor();
  }
#endif

  init_pass_map(h_pass_map);

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, h_pass_map, h_pass_map + 1, d_pass_map);

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

#ifdef NUM_ADV
  NumAdv = NUM_ADV;
#else
  NumAdv = 0;
#endif

  if (NumAdv > 0) {
    FirstAdv = cnt;
    cnt += NumAdv;
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

  // NVAR = cnt;

#ifdef AMREX_PARTICLES
  // Set index locations for particle state vector
  pstateVel = 0;
  pstateT = pstateVel + AMREX_SPACEDIM;
  pstateDia = pstateT + 1;
  pstateRho = pstateDia + 1;
  pstateY = pstateRho + 1;
  pstateNum = pstateY + SPRAY_FUEL_NUM;
#endif

  // const amrex::Real run_strt = amrex::ParallelDescriptor::second() ;
  // Real run_stop = ParallelDescriptor::second() - run_strt;
  // ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());

  // if (ParallelDescriptor::IOProcessor())
  //    amrex::Print() << "\nTime in set_method_params: " << run_stop << '\n' ;

  // if (nscbc_adv == 1 && amrex::ParallelDescriptor::IOProcessor()) {
  //  amrex::Print() << "Using Ghost-Cells Navier-Stokes Characteristic BCs for
  //  "
  //                    "advection: nscbc_adv = "
  //                 << nscbc_adv << '\n'
  //                 << '\n';
  //}

  // if (nscbc_diff == 1 && amrex::ParallelDescriptor::IOProcessor()) {
  //  amrex::Print() << "Using Ghost-Cells Navier-Stokes Characteristic BCs for
  //  "
  //                    "diffusion: nscbc_diff = "
  //                 << nscbc_diff << '\n'
  //                 << '\n';
  //}

  // int coord_type = amrex::DefaultGeometry().Coord();

  amrex::Vector<amrex::Real> center(AMREX_SPACEDIM, 0.0);
  amrex::ParmParse ppc("pelec");
  ppc.queryarr("center", center, 0, AMREX_SPACEDIM);

  amrex::Interpolater* interp;

  if (state_interp_order == 0) {
    interp = &amrex::pc_interp;
  } else {
    if (lin_limit_state_interp == 1) {
      interp = &amrex::lincc_interp;
    } else {
      interp = &amrex::cell_cons_interp;
    }
  }

#ifdef PELEC_USE_EB
  // Decide from input whether eb interp is needed for FillPatch operations
  if ((eb_in_domain) and (state_interp_order != 0)) {
    interp = &amrex::eb_cell_cons_interp;
  }
#endif

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = state_nghost;
  AMREX_ASSERT(ngrow_state >= 0);

  store_in_checkpoint = true;
  desc_lst.addDescriptor(
    State_Type, amrex::IndexType::TheCellType(), amrex::StateDescriptor::Point,
    ngrow_state, NVAR, interp, state_data_extrap, store_in_checkpoint);

  // Components 0:Numspec-1 are rho.omega_i
  // Component NUM_SPECIES is rho.edot = (rho.eout-rho.ein)
#ifdef PELEC_USE_REACTIONS
  store_in_checkpoint = true;
  desc_lst.addDescriptor(
    Reactions_Type, amrex::IndexType::TheCellType(),
    amrex::StateDescriptor::Point, 0, NUM_SPECIES + 2, interp,
    state_data_extrap, store_in_checkpoint);
#endif

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

  for (int i = 0; i < NumAdv; ++i) {
    char buf[64];
    sprintf(buf, "adv_%d", i);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = std::string(buf);
  }

  // Get the species names from the network model
  // Set it for Null mechanism let it be overwritten for others
  spec_names.resize(1);
  spec_names[0] = "Null";
  CKSYMS_STR(spec_names);

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

    // Disabling for the GPU at the moment. Look at the species names to see how
    // to do this in C++. AUX stuff is usually 0 anyway. This call returns the
    // actual length of each string in "len"
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

#ifdef PELEC_USE_REACTIONS
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
#endif

  if (do_react_load_balance || do_mol_load_balance) {
    desc_lst.addDescriptor(
      Work_Estimate_Type, amrex::IndexType::TheCellType(),
      amrex::StateDescriptor::Point, 0, 1, &amrex::pc_interp);
    // Because we use piecewise constant interpolation, we do not use bc and
    // BndryFunc.
    desc_lst.setComponent(
      Work_Estimate_Type, 0, "WorkEstimate", bc,
      amrex::StateDescriptor::BndryFunc(pc_nullfill));
  }

  num_state_type = desc_lst.size();

  // DEFINE DERIVED QUANTITIES

  // Pressure
  derive_lst.add(
    "pressure", amrex::IndexType::TheCellType(), 1, pc_derpres, the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, Density, NVAR);

  // Kinetic energy
  derive_lst.add(
    "kineng", amrex::IndexType::TheCellType(), 1, pc_derkineng, the_same_box);
  derive_lst.addComponent("kineng", desc_lst, State_Type, Density, NVAR);

  // Enstrophy
  derive_lst.add(
    "enstrophy", amrex::IndexType::TheCellType(), 1, pc_derenstrophy,
    grow_box_by_one);
  derive_lst.addComponent("enstrophy", desc_lst, State_Type, Density, NVAR);

  // Sound speed (c)
  derive_lst.add(
    "soundspeed", amrex::IndexType::TheCellType(), 1, pc_dersoundspeed,
    the_same_box);
  derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density, NVAR);

  // Mach number(M)
  derive_lst.add(
    "MachNumber", amrex::IndexType::TheCellType(), 1, pc_dermachnumber,
    the_same_box);
  derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density, NVAR);

  // Entropy (S)
  derive_lst.add(
    "entropy", amrex::IndexType::TheCellType(), 1, pc_derentropy, the_same_box);
  derive_lst.addComponent("entropy", desc_lst, State_Type, Density, NVAR);

  // Vorticity
  derive_lst.add(
    "magvort", amrex::IndexType::TheCellType(), 1, pc_dermagvort,
    grow_box_by_one);
  derive_lst.addComponent("magvort", desc_lst, State_Type, Density, NVAR);

  // Div(u)
  derive_lst.add(
    "divu", amrex::IndexType::TheCellType(), 1, pc_derdivu, grow_box_by_one);
  derive_lst.addComponent("divu", desc_lst, State_Type, Density, NVAR);

  // Internal energy as derived from rho*E, part of the state
  derive_lst.add(
    "eint_E", amrex::IndexType::TheCellType(), 1, pc_dereint1, the_same_box);
  derive_lst.addComponent("eint_E", desc_lst, State_Type, Density, NVAR);

  // Internal energy as derived from rho*e, part of the state
  derive_lst.add(
    "eint_e", amrex::IndexType::TheCellType(), 1, pc_dereint2, the_same_box);
  derive_lst.addComponent("eint_e", desc_lst, State_Type, Density, NVAR);

  // Log(density)
  derive_lst.add(
    "logden", amrex::IndexType::TheCellType(), 1, pc_derlogden, the_same_box);
  derive_lst.addComponent("logden", desc_lst, State_Type, Density, NVAR);

  // Y from rhoY
  amrex::Vector<std::string> var_names_massfrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    var_names_massfrac[i] = "Y(" + spec_names[i] + ")";
  }

  derive_lst.add(
    "massfrac", amrex::IndexType::TheCellType(), NUM_SPECIES,
    var_names_massfrac, pc_derspec, the_same_box);
  derive_lst.addComponent("massfrac", desc_lst, State_Type, Density, NVAR);

  // Species mole fractions
  amrex::Vector<std::string> var_names_molefrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    var_names_molefrac[i] = "X(" + spec_names[i] + ")";
  }

  derive_lst.add(
    "molefrac", amrex::IndexType::TheCellType(), NUM_SPECIES,
    var_names_molefrac, pc_dermolefrac, the_same_box);
  derive_lst.addComponent("molefrac", desc_lst, State_Type, Density, NVAR);

  // Velocities
  derive_lst.add(
    "x_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelx, the_same_box);
  derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "y_velocity", amrex::IndexType::TheCellType(), 1, pc_dervely, the_same_box);
  derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "z_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelz, the_same_box);
  derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "magvel", amrex::IndexType::TheCellType(), 1, pc_dermagvel, the_same_box);
  derive_lst.addComponent("magvel", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "radvel", amrex::IndexType::TheCellType(), 1, pc_derradialvel,
    the_same_box);
  derive_lst.addComponent("radvel", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "magmom", amrex::IndexType::TheCellType(), 1, pc_dermagmom, the_same_box);
  derive_lst.addComponent("magmom", desc_lst, State_Type, Density, NVAR);

#ifdef PELEC_USE_EB
  // A dummy
  derive_lst.add(
    "vfrac", amrex::IndexType::TheCellType(), 1, pc_dermagvel, the_same_box);
#endif

#ifdef AMREX_PARTICLES
  // We want a derived type that corresponds to the number of particles
  // in each cell.  We only intend to use it in plotfiles for debugging
  // purposes. We'll actually set the values in writePlotFile().
  derive_lst.add(
    "particle_count", amrex::IndexType::TheCellType(), 1, pc_dernull,
    the_same_box);
  derive_lst.addComponent("particle_count", desc_lst, State_Type, Density, 1);

  derive_lst.add(
    "total_particle_count", amrex::IndexType::TheCellType(), 1, pc_dernull,
    the_same_box);
  derive_lst.addComponent(
    "total_particle_count", desc_lst, State_Type, Density, 1);

  derive_lst.add(
    "particle_density", amrex::IndexType::TheCellType(), 1, pc_dernull,
    the_same_box);
  derive_lst.addComponent("particle_density", desc_lst, State_Type, Density, 1);
#endif

  // LES coefficients
  if (do_les) {
    derive_lst.add(
      "C_s2", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("C_s2", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "C_I", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("C_I", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "Pr_T", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("Pr_T", desc_lst, State_Type, Density, 1);
  }

  // MMS derives
#ifdef PELEC_USE_MASA
  if (do_mms) {
    derive_lst.add(
      "rhommserror", amrex::IndexType::TheCellType(), 1, pc_derrhommserror,
      the_same_box);
    derive_lst.addComponent("rhommserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "ummserror", amrex::IndexType::TheCellType(), 1, pc_derummserror,
      the_same_box);
    derive_lst.addComponent("ummserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "vmmserror", amrex::IndexType::TheCellType(), 1, pc_dervmmserror,
      the_same_box);
    derive_lst.addComponent("vmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "wmmserror", amrex::IndexType::TheCellType(), 1, pc_derwmmserror,
      the_same_box);
    derive_lst.addComponent("wmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "pmmserror", amrex::IndexType::TheCellType(), 1, pc_derpmmserror,
      the_same_box);
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

  pele::physics::transport::CloseTransport<
    pele::physics::PhysicsType::eos_type>()();

#ifdef PELEC_USE_REACTIONS
  if (do_react == 1) {
    close_reactor();
  }
#endif

  clear_prob();

#ifdef PELEC_USE_EB
  eb_initialized = false;
#endif

  delete prob_parm_host;
  delete tagging_parm;
  delete h_prob_parm_device;
  delete h_pass_map;
  amrex::The_Arena()->free(d_prob_parm_device);
  amrex::The_Arena()->free(d_pass_map);
}

void
PeleC::set_active_sources()
{
  if (do_diffuse && !do_mol) {
    src_list.push_back(diff_src);
  }

  // optional external source
  if (add_ext_src == 1) {
    src_list.push_back(ext_src);
  }

  // optional forcing source
  if (add_forcing_src == 1) {
    src_list.push_back(forcing_src);
  }

#ifdef AMREX_PARTICLES
  if (do_spray_particles) {
    src_list.push_back(spray_src);
  }
#endif

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
