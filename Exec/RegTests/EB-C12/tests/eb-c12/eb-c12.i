# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10000000
stop_time = 2.0

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic =  1   1  1
geometry.coord_sys   =  0       # 0 => cart
geometry.prob_lo     = -1.0 -0.2 -0.2
geometry.prob_hi     = 1.0 0.2 0.2
amr.n_cell           = 40 8 8


# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       =  "Interior"  "Interior"  "Interior"
pelec.hi_bc       =  "Interior"  "Interior"  "Interior"

# Problem setup
pelec.eb_boundary_T = 300.
pelec.eb_isothermal = 0
eb_verbosity = 1
eos_gamma = 1.4

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.do_mol = 1
pelec.do_react = 0
pelec.allow_negative_energy = 0
pelec.diffuse_temp = 0
pelec.diffuse_vel  = 0
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
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 16

# CHECKPOINT FILES
amr.checkpoint_files_output = 1
amr.check_file      = chk      # root name of checkpoint file
amr.check_int       = 1000       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt
amr.plot_int        = 1000
amr.derive_plot_vars=ALL

# EB
#eb2.geom_type = "all_regular"
eb2.geom_type = "cylinder"
eb2.cylinder_direction = 0
eb2.cylinder_center = 0.0 0.0 0.0
eb2.cylinder_radius = 0.1
eb2.cylinder_height = 1000.0
eb2.cylinder_has_fluid_inside = 1
ebd.boundary_grad_stencil_type = 0
