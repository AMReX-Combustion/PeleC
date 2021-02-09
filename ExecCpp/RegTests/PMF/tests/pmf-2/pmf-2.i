# ------------------  INPUTS TO MAIN PROGRAM  -------------------
stop_time = 6
max_step = 10

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     =   0.0        0.0       1.0
geometry.prob_hi     =   0.3125     0.3125    6.0
amr.n_cell           =   8          8         128

#pelec.Riemann    = 0     # 0: HLL,  1: JBB,  2: HLLC
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       =  "Interior"  "Interior"  "Hard"
pelec.hi_bc       =  "Interior"  "Interior"  "Hard"

# TIME STEP CONTROL
pelec.cfl            = 0.1     # cfl number for hyperbolic system
pelec.init_shrink    = 0.1     # scale back initial timestep
pelec.change_max     = 1.1     # scale back initial timestep
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval = 1       # coarse time steps between computing mass on domain
pelec.v            = 1       # verbosity in PeleC cpp files
amr.v              = 1       # verbosity in Amr.cpp
#amr.grid_log       = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING 
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 8       # block factor in grid generation
amr.max_grid_size   = 32
amr.n_error_buf     = 2 8 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file              = chk    # root name of checkpoint file
amr.check_int               = 500    # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file         = plt     # root name of plotfile
amr.plot_int          = 10   # number of timesteps between plotfiles

# PROBLEM PARAMETERS
prob.pamb = 1013250.0  
prob.phi_in = -0.5
prob.pertmag = 0.005
prob.pmf_datafile = "LiDryer_H2_p1_phi0_4000tu0300.dat"
#prob.pmf_datafile = "PMF_CH4_1bar_300K_DRM_MixAvg.dat"

tagging.max_ftracerr_lev = 4
tagging.ftracerr = 150.e-6

extern.new_Jacobian_each_cell = 0

amr.derive_plot_vars = density xmom ymom zmom eden Temp pressure x_velocity y_velocity z_velocity
pelec.plot_rhoy = 0
pelec.plot_massfrac = 1
pelec.do_react = 1
pelec.diffuse_temp=1
pelec.diffuse_enth=1
pelec.diffuse_spec=1
pelec.diffuse_vel=1
pelec.sdc_iters = 2
pelec.flame_trac_name = HO2
pelec.do_mol=0

eb2.use_eb2 = 1
eb2.geom_type = "all_regular"
ebd.boundary_grad_stencil_type = 0

pelec.chem_integrator=2
