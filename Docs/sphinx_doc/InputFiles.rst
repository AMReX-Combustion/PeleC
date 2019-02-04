Input Files and Controls
========================


The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to specified along with the execuatable as an `argv` option, for example:


::

	mpirun -np 64 ./Pele2d.xxx,yyy.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inpute file, e.g.:

::

	mpirun -np 64 ./Pele2d.gnu.DEBUG.MPI.ex inputs amr.restart = sod_x_chk0030 pelec.riemann_solver=3

The available options are divided into groups: those that control primarily AMReX are prefaced with `amr.` while those that are specific to Pele are prefaced with `pelec.`.

A typical input file looks something like the example below; a full list of Pele-specific input parameters are in `PeleC/Source/_cpp_parameters`. 
These parameters, once read, are available in the `PeleC` object for use from c++ and are also copied to the module `prob_params_module` for use in FORTRAN. 

::

    # ------------------  INPUTS TO MAIN PROGRAM  -------------------
    #absolute stop time for the simulation
    stop_time = 6 

    #maximum number of time steps at base AMR level
    max_step = 30 
    # ---------------------------------------------------------------
    
    #------------------------
    # PROBLEM SIZE & GEOMETRY
    # -----------------------

    #flag for periodicity (here x direction is periodic)
    geometry.is_periodic = 1 0 0  
    
    #0 => cart, 1 => RZ  2=>spherical
    geometry.coord_sys   = 0      

    #coordinates of domain's lower corner
    geometry.prob_lo     =   -0.3     0.0   0.0     

    #coordiates of domain's upper corner
    geometry.prob_hi     =    0.3     0.3   0.15  

    #number of cells along each direction at base level (note: dx=dy=dz)
    amr.n_cell           =    128     64    32   
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    PeleC specific inputs
    # ---------------------------------------------------------------

    # 0: Collela, Glaz and Ferguson (default)
    # 1: Collela and Glaz  
    # 2: HLLC
    pelec.riemann_solver    = 0     

    # >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
    # 0 = Interior           3 = Symmetry
    # 1 = Inflow             4 = SlipWall
    # 2 = Outflow            5 = NoSlipWall
    # >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<

    #boundary condition at the lower face of each coordinate direction
    pelec.lo_bc       =  0   1  0          
    
    #boundary condition at the upper face of each coordinate direction
    pelec.hi_bc       =  0   2  0          
    
    # TIME STEP CONTROL
    pelec.cfl            = 0.5     # cfl number for hyperbolic system
    pelec.init_shrink    = 0.3     # scale back initial timestep
    pelec.change_max     = 1.1     # scale back initial timestep
    pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt
    
    # DIAGNOSTICS & VERBOSITY
    pelec.sum_interval = 1        # coarse time steps between computing mass on domain
    pelec.v            = 1        # verbosity in PeleC cpp files
    amr.v              = 1        # verbosity in Amr.cpp
    #amr.grid_log       = grdlog  # name of grid logging file
    
    # REFINEMENT / REGRIDDING 
    amr.max_level       = 2       # maximum level number allowed
    amr.ref_ratio       = 2 2 2 2 # refinement ratio
    amr.regrid_int      = 2 2 2 2 # how often to regrid
    amr.blocking_factor = 8       # block factor in grid generation
    amr.max_grid_size   = 64
    amr.n_error_buf     = 2 8 2 2 # number of buffer cells in error est
    
    # CHECKPOINT FILES
    amr.checkpoint_files_output = 1
    amr.check_file              = chk    # root name of checkpoint file
    amr.check_int               = 500    # number of timesteps between     checkpoints
    
    # PLOTFILES
    amr.plot_files_output = 1
    amr.plot_file         = plt     # root name of plotfile
    amr.plot_int          = 100     # number of timesteps between plotfiles
    
    #PROBIN FILENAME
    amr.probin_file = probin-regt
    
    amr.derive_plot_vars = density xmom ymom zmom eden Temp pressure     x_velocity y_velocity
    pelec.plot_rhoy = 0
    pelec.plot_massfrac = 1
    pelec.do_react = 1
    pelec.diffuse_temp=1
    pelec.diffuse_enth=1
    pelec.diffuse_spec=1
    pelec.diffuse_vel=1
    pelec.sdc_iters = 2
    pelec.flame_trac_name = HO2
