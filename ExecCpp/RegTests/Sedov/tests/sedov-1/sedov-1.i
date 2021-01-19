# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 250

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical

geometry.prob_lo     =  0.   0.   0.
geometry.prob_hi     =  1.   1.   1.

#amr.n_cell           = 128  128  128
amr.n_cell           = 64 64 64

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
pelec.cfl            = 0.5     # cfl number for hyperbolic system
pelec.init_shrink    = 0.01    # scale back initial timestep
pelec.change_max     = 1.1     # scale back initial timestep
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 10      # timesteps between computing mass
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
amr.plot_files_output = 0
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100        # number of timesteps between plotfiles
amr.plot_vars  =  density Temp
amr.derive_plot_vars = x_velocity y_velocity z_velocity magvel magvort pressure

# PROBLEM PARAMETERS
prob.r_init = 0.01
prob.p_ambient = 1.e-5
prob.dens_ambient = 1.0
prob.exp_energy = 1.0
prob.nsub = 10

eb2.use_eb2 = 1
eb2.geom_type = "all_regular"
ebd.boundary_grad_stencil_type = 0
