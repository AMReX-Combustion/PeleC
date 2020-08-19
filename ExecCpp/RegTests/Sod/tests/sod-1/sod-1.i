# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#stop_time =  0.2
max_step = 10

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     =  0     0     0
geometry.prob_hi     =  1     0.25  0.25
amr.n_cell           = 32     8     8

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       =     "UserBC"   "SlipWall"     "SlipWall"
pelec.hi_bc       =     "UserBC"   "SlipWall"     "SlipWall"

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.do_react = 0
pelec.ppm_type = 1

# TIME STEP CONTROL
pelec.cfl            = 0.9     # cfl number for hyperbolic system
pelec.init_shrink    = 0.1     # scale back initial timestep
pelec.change_max     = 1.05    # scale back initial timestep
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in PeleC cpp files
amr.v                 = 1       # verbosity in Amr.cpp
#amr.grid_log        = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING 
amr.max_level       = 2       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 8       # block factor in grid generation
amr.max_grid_size   = 64
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 10         # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file         = plt      # root name of plotfile
amr.plot_int          = 10       # number of timesteps between plotfiles
amr.derive_plot_vars  = ALL # density xmom ymom zmom eden Temp pressure  # these variables appear in the plotfile

# PROBLEM PARAMETERS
prob.p_l = 1.0
prob.u_l = 0.0
prob.rho_l = 1.0
prob.p_r = 0.1
prob.u_r = 0.0
prob.rho_r = 0.125
prob.idir = 1
prob.frac = 0.5

# TAGGING
tagging.denerr = 3
tagging.dengrad = 0.01
tagging.max_denerr_lev = 3
tagging.max_dengrad_lev = 3
tagging.presserr = 3
tagging.pressgrad = 0.01
tagging.max_presserr_lev = 3
tagging.max_pressgrad_lev = 3

# EB
eb2.geom_type = "all_regular"
ebd.boundary_grad_stencil_type = 0
