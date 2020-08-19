# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10000
stop_time =  1.959e-6 #final time is 0.2*L*sqrt(rhoL/pL)

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     =  0     0     0
geometry.prob_hi     =  1     0.0625 0.0625
amr.n_cell           = 128     8     8

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       = "SlipWall"   "SlipWall"   "SlipWall"
pelec.hi_bc       = "SlipWall"   "SlipWall"   "SlipWall"

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.diffuse_vel = 0
pelec.diffuse_temp = 0
pelec.diffuse_spec = 0
pelec.do_react = 0

# TIME STEP CONTROL
pelec.cfl            = 0.9     # cfl number for hyperbolic system
pelec.init_shrink    = 0.1     # scale back initial timestep
pelec.change_max     = 1.05    # scale back initial timestep
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in Castro.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog

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
amr.plot_int        = 100        # number of timesteps between plotfiles
amr.plot_vars  =  density Temp
amr.derive_plot_vars = x_velocity y_velocity z_velocity magvel magvort pressure

# PROBLEM PARAMETERS
prob.p_l = 1e7
prob.u_l = 0.0
prob.rho_l = 9.6e-4
prob.p_r = 1e6
prob.u_r = 0.0
prob.rho_r = 1.2e-4
prob.idir = 1
prob.frac = 0.5
prob.left_gas = N2
prob.right_gas = HE
