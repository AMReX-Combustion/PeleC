#include <cstdio>

#include "AMReX_LevelBld.H"
#include <AMReX_ParmParse.H>
#include "PeleC.H"
#include "PeleC_F.H"
#include <Derive_F.H>
#include "AMReX_buildInfo.H"

#ifdef USE_MASA
#include <masa.h>
using namespace MASA;
#endif

using std::string;
using namespace amrex;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

typedef StateDescriptor::BndryFunc BndryFunc;


//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall, UserBC
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD, EXT_DIR
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_ODD, EXT_DIR
};

static int react_src_bc[] =
{
    INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	bc.setLo(i,scalar_bc[lo_bc[i]]);
	bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)    
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_react_src_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	bc.setLo(i,react_src_bc[lo_bc[i]]);
	bc.setHi(i,react_src_bc[hi_bc[i]]);
    }
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

void
PeleC::variableSetUp ()
{

    // PeleC::variableSetUp is called in the constructor of Amr.cpp, so
    // it should get called every time we start or restart a job


    // initialize the start time for our CPU-time tracker
    startCPUTime = ParallelDescriptor::second();


    // Output the git commit hashes used to build the executable.

    if (ParallelDescriptor::IOProcessor()) {

	const char* pelec_hash  = buildInfoGetGitHash(1);
	const char* amrex_hash  = buildInfoGetGitHash(2);
	const char* pelephysics_hash = buildInfoGetGitHash(3);
	const char* buildgithash = buildInfoGetBuildGitHash();
	const char* buildgitname = buildInfoGetBuildGitName();

	if (strlen(pelec_hash) > 0) {
	    std::cout << "\n" << "PeleC git hash: " << pelec_hash << "\n";
	}
	if (strlen(amrex_hash) > 0) {
	    std::cout << "AMReX git hash: " << amrex_hash << "\n";
	}
	if (strlen(pelephysics_hash) > 0) {
	    std::cout << "PelePhysics git hash: " << pelephysics_hash << "\n";
	}
	if (strlen(buildgithash) > 0){
	    std::cout << buildgitname << " git hash: " << buildgithash << "\n";
	}
    
	std::cout << "\n";
    }
  
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    // Initialize the runtime parameters for any of the external code
    init_extern();

    // Initialize the network
    init_network();

#ifdef REACTIONS
    // Initialize the reactor
    if (do_react == 1) {
	init_reactor();
    }
#endif

    init_transport();

#ifdef USE_MASA
    if (do_mms) {
      init_mms();
    }
#endif

    //
    // Set number of state variables and pointers to components
    //

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
    
    if (NumAdv > 0)
    {
	FirstAdv = cnt;
	cnt += NumAdv;
    }

    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    get_num_spec(&NumSpec);
  
    if (NumSpec > 0)
    {
	FirstSpec = cnt;
	cnt += NumSpec;
    }

    // Get the number of auxiliary quantities from the network model.
    get_num_aux(&NumAux);
  
    if (NumAux > 0)
    {
	FirstAux = cnt;
	cnt += NumAux;
    }

    NUM_STATE = cnt;


#ifdef AMREX_PARTICLES
    // Set index locations for particle state vector and storage of field variables that 
    // get computed first and then interpolated to particle position
    pstate_loc = 0;
    pstate_vel = pstate_loc + BL_SPACEDIM;
    pstate_T = pstate_vel + BL_SPACEDIM;
    pstate_dia = pstate_T + 1;
    pstate_rho = pstate_dia + 1;
    pstate_spc = pstate_rho + 1;
    n_pstate = pstate_spc + 1;

    pfld_vel = 0;
    pfld_rho = pfld_vel + BL_SPACEDIM;
    pfld_T = pfld_rho + 1;
    pfld_p = pfld_T + 1;
    pfld_spc = pfld_p + 1;
    n_pfld = pfld_spc + NumSpec; // increase this for multi-component droplets

//  std::cout << "n_pstate = " << n_pstate << std::endl << "declared components: " << SPRAY_COMPONENTS << std::endl;

#endif
    

    const Real run_strt = ParallelDescriptor::second() ; 

    // Read in the input values to Fortran.

    set_pelec_method_params();

    set_method_params(dm, Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux, 
		      NumAdv,
		      diffuse_cutoff_density,
              pstate_loc, pstate_vel, pstate_T, pstate_dia, pstate_rho, pstate_spc,
              pfld_vel, pfld_rho, pfld_T, pfld_p, pfld_spc);

    // Get various values from Fortran
    get_method_params(&NUM_GROW,&QTHERM,&QVAR,&cQRHO,&cQU,&cQV,&cQW,&cQGAME,&cQPRES,
                      &cQREINT,&cQTEMP,&cQFA,&cQFS,&cQFX,&NQAUX,&cQGAMC,&cQC,&cQCSML,
                      &cQDPDR,&cQDPDE,&cQRSPEC);

    Real run_stop = ParallelDescriptor::second() - run_strt;
 
    ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
 
    if (ParallelDescriptor::IOProcessor())
	std::cout << "\nTime in set_method_params: " << run_stop << '\n' ;
  
    if (nscbc_adv == 1 && ParallelDescriptor::IOProcessor())
	std::cout << '\n' << "Using Ghost-Cells Navier-Stokes Characteristic BCs for advection: nscbc_adv = " << nscbc_adv << '\n' << '\n' ;

    if (nscbc_diff == 1 && ParallelDescriptor::IOProcessor())
	std::cout << "Using Ghost-Cells Navier-Stokes Characteristic BCs for diffusion: nscbc_diff = " << nscbc_diff << '\n' << '\n' ;

  
    int coord_type = DefaultGeometry().Coord();

    // Get the center variable from the inputs and pass it directly to Fortran.
    Vector<Real> center(BL_SPACEDIM, 0.0);
    ParmParse ppc("pelec");
    ppc.queryarr("center",center,0,BL_SPACEDIM);
  
    set_problem_params(dm,phys_bc.lo(),phys_bc.hi(),
		       Interior,UserBC,Inflow,Outflow,Symmetry,SlipWall,NoSlipWall,coord_type,
		       DefaultGeometry().ProbLo(),DefaultGeometry().ProbHi(),center.dataPtr());
  
    // Read in the parameters for the tagging criteria
    // and store them in the Fortran module.
  
    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);
  
    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];
  
    get_tagging_params(probin_file_name.dataPtr(),&probin_file_length);
  
    Interpolater* interp;


    if (state_interp_order == 0)
    {
        interp = &pc_interp;
    }
    else
    {
        if (lin_limit_state_interp == 1)
        {
            interp = &lincc_interp;
        }
        else {
            interp = &cell_cons_interp;
        }
    }

#ifdef AMREX_USE_EB
    // Decide from input whether eb interp is needed for FillPatch oeprations
    ParmParse pp;
    std::string geom_type("all_regular");
    pp.query("geom_type", geom_type);
    if (geom_type != "all_regular") {
        if (state_interp_order != 0) {
            interp = &eb_cell_cons_interp;
        }
    }
#endif

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    int ngrow_state = state_nghost;
    BL_ASSERT(ngrow_state >= 0);

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			   StateDescriptor::Point,ngrow_state,NUM_STATE,
			   interp,state_data_extrap,store_in_checkpoint);

    // Components 0:Numspec-1         are      rho.omega_i
    // Component    NumSpec            is      rho.edot = (rho.eout-rho.ein)
#ifdef REACTIONS
    store_in_checkpoint = true;
    desc_lst.addDescriptor(Reactions_Type,IndexType::TheCellType(),
			   StateDescriptor::Point,0,NumSpec+1,
			   &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

    Vector<BCRec>       bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);
    Vector<BCRec>       react_bcs(NumSpec+1);
    Vector<std::string> react_name(NumSpec+1);

    BCRec bc;
    cnt = 0;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";

    for (int i=0; i<NumAdv; ++i)
    {
	char buf[64];
	sprintf(buf, "adv_%d", i);
	cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = string(buf);
    }

    // Get the species names from the network model.
    for (int i = 0; i < NumSpec; i++) {
	int len = 20;
	Vector<int> int_spec_names(len);
	// This call return the actual length of each string in "len" 
	get_spec_names(int_spec_names.dataPtr(),&i,&len);
	char char_spec_names[len+1];
	for (int j = 0; j < len; j++) 
	    char_spec_names[j] = int_spec_names[j];
	char_spec_names[len] = '\0';
	spec_names.push_back(std::string(char_spec_names));
    }

    if ( ParallelDescriptor::IOProcessor())
    {
	std::cout << NumSpec << " Species: " << std::endl;
	for (int i = 0; i < NumSpec; i++)  
	    std::cout << spec_names[i] << ' ' << ' ';
	std::cout << std::endl;
    } 

    for (int i=0; i<NumSpec; ++i)
    {
	cnt++; 
	set_scalar_bc(bc,phys_bc); 
	bcs[cnt] = bc; 
	name[cnt] = "rho_" + spec_names[i];
    }

    // Get the auxiliary names from the network model.
    std::vector<std::string> aux_names;
    for (int i = 0; i < NumAux; i++) {
	int len = 20;
	Vector<int> int_aux_names(len);
	// This call return the actual length of each string in "len"
	get_aux_names(int_aux_names.dataPtr(),&i,&len);
	char char_aux_names[len+1];
	for (int j = 0; j < len; j++)
	    char_aux_names[j] = int_aux_names[j];
	char_aux_names[len] = '\0';
	aux_names.push_back(std::string(char_aux_names));
    }

    if ( ParallelDescriptor::IOProcessor())
    {
	std::cout << NumAux << " Auxiliary Variables: " << std::endl;
	for (int i = 0; i < NumAux; i++)
	    std::cout << aux_names[i] << ' ' << ' ';
	std::cout << std::endl;
    }

    for (int i=0; i<NumAux; ++i)
    {
	cnt++;
	set_scalar_bc(bc,phys_bc);
	bcs[cnt] = bc;
	name[cnt] = "rho_" + aux_names[i];
    }

  desc_lst.setComponent(State_Type,
                        Density,
                        name,
                        bcs,
                        StateDescriptor::BndryFunc(pc_bcfill_hyp));

#ifdef REACTIONS
    for (int i=0; i<NumSpec; ++i) {
      set_react_src_bc(bc, phys_bc);
      react_bcs[i] = bc;
      react_name[i] = "rho_omega_" + spec_names[i];
    }
   set_react_src_bc(bc, phys_bc);
   react_bcs[NumSpec] = bc;
   react_name[NumSpec] = "rhoe_dot";

   desc_lst.setComponent(Reactions_Type,
                         0,
                         react_name,
                         react_bcs,
                         StateDescriptor::BndryFunc(pc_reactfill_hyp));
#endif


    if (do_react_load_balance || do_mol_load_balance) {
      desc_lst.addDescriptor(Work_Estimate_Type, IndexType::TheCellType(),
                             StateDescriptor::Point, 0, 1, &pc_interp);
      // Because we use piecewise constant interpolation, we do not use bc and BndryFunc.
      desc_lst.setComponent(Work_Estimate_Type, 0, "WorkEstimate",
                            bc, BndryFunc(pc_nullfill));
    }

    num_state_type = desc_lst.size();

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure",IndexType::TheCellType(),1,pc_derpres,the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng",IndexType::TheCellType(),1,pc_derkineng,the_same_box);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Enstrophy
    //
    derive_lst.add("enstrophy",IndexType::TheCellType(),1,pc_derenstrophy,grow_box_by_one);
    derive_lst.addComponent("enstrophy",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed",IndexType::TheCellType(),1,pc_dersoundspeed,the_same_box);
    derive_lst.addComponent("soundspeed",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber",IndexType::TheCellType(),1,pc_dermachnumber,the_same_box);
    derive_lst.addComponent("MachNumber",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM == 1)
    //
    // Wave speed u+c
    //
    derive_lst.add("uplusc",IndexType::TheCellType(),1,pc_deruplusc,the_same_box);
    derive_lst.addComponent("uplusc",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Wave speed u-c
    //
    derive_lst.add("uminusc",IndexType::TheCellType(),1,pc_deruminusc,the_same_box);
    derive_lst.addComponent("uminusc",desc_lst,State_Type,Density,NUM_STATE);
#endif

    //
    // Entropy (S)
    //
    derive_lst.add("entropy",IndexType::TheCellType(),1,pc_derentropy,the_same_box);
    derive_lst.addComponent("entropy",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Vorticity
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,pc_dermagvort,grow_box_by_one);
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Div(u)
    //
    derive_lst.add("divu",IndexType::TheCellType(),1,pc_derdivu,grow_box_by_one);
    derive_lst.addComponent("divu",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E",IndexType::TheCellType(),1,pc_dereint1,the_same_box);
    derive_lst.addComponent("eint_E",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e",IndexType::TheCellType(),1,pc_dereint2,the_same_box);
    derive_lst.addComponent("eint_e",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden",IndexType::TheCellType(),1,pc_derlogden,the_same_box);
    derive_lst.addComponent("logden",desc_lst,State_Type,Density,NUM_STATE);
    
#ifdef DO_HIT_FORCE
    //
    // forcing - used to calculate the rate of injection of energy in probtype 14 (HIT)
    //

    derive_lst.add("forcing",IndexType::TheCellType(),1,pc_derforcing,the_same_box);
    derive_lst.addComponent("forcing",desc_lst,State_Type,Density,NUM_STATE);
    //
    // forcex - used to put the forcing term in the plot file
    //
    derive_lst.add("forcex",IndexType::TheCellType(),1,pc_derforcex,the_same_box);
    derive_lst.addComponent("forcex",desc_lst,State_Type,Density,NUM_STATE);
    //
    // forcey - used to put the forcing term in the plot file
    //
    derive_lst.add("forcey",IndexType::TheCellType(),1,pc_derforcey,the_same_box);
    derive_lst.addComponent("forcey",desc_lst,State_Type,Density,NUM_STATE);
    //
    // forcez - used to put the forcing term in the plot file
    //
    derive_lst.add("forcez",IndexType::TheCellType(),1,pc_derforcez,the_same_box);
    derive_lst.addComponent("forcez",desc_lst,State_Type,Density,NUM_STATE);
#endif

    //
    // Y from rhoY
    //
    Vector<std::string> var_names_massfrac(NumSpec);
    for (int i = 0; i < NumSpec; i++){
      var_names_massfrac[i] = "Y("+spec_names[i]+")";
    }

    derive_lst.add("massfrac",IndexType::TheCellType(),NumSpec,var_names_massfrac,
                   pc_derspec,the_same_box);
    derive_lst.addComponent("massfrac",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Species mole fractions
    //
    
    Vector<std::string> var_names_molefrac(NumSpec);
    for (int i = 0; i < NumSpec; i++){
      var_names_molefrac[i] = "X("+spec_names[i]+")";
    }  
     
    derive_lst.add("molefrac",IndexType::TheCellType(),NumSpec,var_names_molefrac,
                   pc_dermolefrac,the_same_box);
    derive_lst.addComponent("molefrac",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Velocities
    //
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,pc_dervelx,the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("y_velocity",IndexType::TheCellType(),1,pc_dervely,the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("z_velocity",IndexType::TheCellType(),1,pc_dervelz,the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("magvel",IndexType::TheCellType(),1,pc_dermagvel,the_same_box);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("radvel",IndexType::TheCellType(),1,pc_derradialvel,the_same_box);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Density,NUM_STATE);

    derive_lst.add("magmom",IndexType::TheCellType(),1,pc_dermagmom,the_same_box);
    derive_lst.addComponent("magmom",desc_lst,State_Type,Density,NUM_STATE);

#ifdef AMREX_USE_EB
    derive_lst.add("vfrac",IndexType::TheCellType(),1,pc_dermagvel,the_same_box); // A dummy
#endif
    
#ifdef AMREX_PARTICLES
    //
    // We want a derived type that corresponds to the number of particles
    // in each cell.  We only intend to use it in plotfiles for debugging
    // purposes.  We'll just use the DERNULL since don't do anything in
    // fortran for now.  We'll actually set the values in writePlotFile().
    //
    derive_lst.add("particle_count",IndexType::TheCellType(),1,pc_dernull,the_same_box);
    derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);

    derive_lst.add("total_particle_count",IndexType::TheCellType(),1,pc_dernull,the_same_box);
    derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);

    derive_lst.add("particle_density",IndexType::TheCellType(),1,pc_dernull,the_same_box);
    derive_lst.addComponent("particle_density",desc_lst,State_Type,Density,1);
#endif

#if 0
    //
    // A derived quantity equal to all the state variables.
    //
    derive_lst.add("FULLSTATE",IndexType::TheCellType(),NUM_STATE,FORT_DERCOPY,the_same_box);
    derive_lst.addComponent("FULLSTATE",desc_lst,State_Type,Density,NUM_STATE);

#endif


    //
    // LES coefficients
    //
    if(do_les){
      derive_lst.add("C_s2",IndexType::TheCellType(),1,pc_dernull,the_same_box);
      derive_lst.addComponent("C_s2",desc_lst,State_Type,Density,1);

      derive_lst.add("C_I",IndexType::TheCellType(),1,pc_dernull,the_same_box);
      derive_lst.addComponent("C_I",desc_lst,State_Type,Density,1);

      derive_lst.add("Pr_T",IndexType::TheCellType(),1,pc_dernull,the_same_box);
      derive_lst.addComponent("Pr_T",desc_lst,State_Type,Density,1);
    }

    // 
    // Problem-specific adds
#include <Problem_Derives.H>

    // Set list of active sources
    set_active_sources();
}

void
PeleC::set_active_sources()
{
    if (do_diffuse && !do_mol_AD){
      src_list.push_back(diff_src);
    }

    // optional external source
    if (add_ext_src == 1){
      src_list.push_back(ext_src);
    }

    // optional forcing source
    if (add_forcing_src == 1){
      src_list.push_back(forcing_src);
    }

#ifdef AMREX_PARTICLES
    if (do_spray_particles) {
      src_list.push_back(spray_src);
    }
#endif

    // optional LES source
    if (do_les){
      src_list.push_back(les_src);
    }

#ifdef USE_MASA
    // optional MMS source
    if(do_mms){
      src_list.push_back(mms_src);
    }
#endif
}
