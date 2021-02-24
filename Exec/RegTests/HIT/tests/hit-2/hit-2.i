# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#stop_time = 0.00026398069024412264
max_step = 10

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 1 1 1
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     =   0.0  0.0  0.0
geometry.prob_hi     =   6.283185307179586232  6.283185307179586232  6.283185307179586232
amr.n_cell           =  32 32 32

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       =  "Interior"  "Interior"  "Interior" 
pelec.hi_bc       =  "Interior"  "Interior"  "Interior"

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.diffuse_vel = 1
pelec.diffuse_temp = 1
pelec.do_react = 0
pelec.do_grav = 0

# TIME STEP CONTROL
pelec.cfl            = 0.9     # cfl number for hyperbolic system
pelec.init_shrink    = 0.3     # scale back initial timestep
pelec.change_max     = 1.1     # max time step growth
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in Castro.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog
#amr.grid_log        = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING 
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 64
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt        # root name of plotfile
amr.plot_per        = 0.00006599517256103066       # time between plotfiles
amr.plot_vars  =  density Temp
amr.derive_plot_vars = x_velocity y_velocity z_velocity magvel magvort pressure

# PROBLEM PARAMETERS
prob.iname = "hit-2.ic"
prob.binfmt = true
prob.lambda0 = 0.2645751311064591
prob.reynolds_lambda0 = 133.6306209562122262
prob.mach_t0 = 0.1
prob.prandtl = 0.71
prob.inres = 32
prob.uin_norm = 1.4142135623730950

# EB
eb2.geom_type = "all_regular"
ebd.boundary_grad_stencil_type = 0
