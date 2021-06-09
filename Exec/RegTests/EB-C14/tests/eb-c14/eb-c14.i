# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 100000
stop_time = 1.6e-5
#stop_time = 0.0034 # 1 flowthroughs (1 flowthrough = L / umax)

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic =  0  0  0
geometry.coord_sys   =  0       # 0 => cart
geometry.prob_lo     =  0   0.0  0.0
#geometry.prob_hi     =  7.5 3.75 0.1171875
#amr.n_cell           =  512 256 8
geometry.prob_hi     =  7.5 3.75 1.875
amr.n_cell           =  32 16 8

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<

pelec.lo_bc       =  "FOExtrap" "NoSlipWall" "FOExtrap"
pelec.hi_bc       =  "FOExtrap" "FOExtrap"   "FOExtrap"

# Problem setup
pelec.eb_boundary_T = 300.
pelec.eb_isothermal = 0

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.do_mol = 1
pelec.do_react = 0
pelec.allow_negative_energy = 0
pelec.diffuse_temp = 1
pelec.diffuse_vel  = 1
pelec.diffuse_spec = 0
pelec.diffuse_enth = 0

# TIME STEP CONTROL
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt
pelec.cfl            = 0.3     # cfl number for hyperbolic system
pelec.init_shrink    = 0.8    # scale back initial timestep
pelec.change_max     = 1.05     # maximum increase in dt over successive steps

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in PeleC cpp files
amr.v                = 1       # verbosity in Amr.cpp
#amr.grid_log         = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2       # how often to regrid
amr.blocking_factor = 8      # block factor in grid generation
amr.max_grid_size   = 64

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk      # root name of checkpoint file
amr.check_int       = -1       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt
amr.plot_int        = 50
amr.derive_plot_vars=pressure x_velocity y_velocity z_velocity

eb2.geom_type = "plane"
eb2.plane_point = 5.2 0.0 0.0
eb2.plane_normal = 0.7547095802227719  -0.6560590289905074  0.0   # slope: 49 degree
ebd.boundary_grad_stencil_type = 0

tagging.max_dengrad_lev = 2
tagging.dengrad = 2.e-5
