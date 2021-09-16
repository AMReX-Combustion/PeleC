# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time =  10e-3

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0  0  0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical

# full geometry:
# geometry.prob_lo     =  -0.15 -4.7025 -4.7025
# geometry.prob_hi     =  24.93 4.7025 4.7025
geometry.prob_lo     =  0.0 -4.7025 -4.7025
geometry.prob_hi     =  14.1075 4.7025 4.7025
amr.n_cell           =  144 96 96
amr.n_cell           =  72 48 48

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
pelec.lo_bc       =  "Hard"  "FOExtrap"  "FOExtrap"
pelec.hi_bc       =  "FOExtrap"  "FOExtrap"  "FOExtrap"

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.do_mol = 1
pelec.diffuse_vel = 0
pelec.diffuse_temp = 0
pelec.diffuse_spec = 0
pelec.do_react = 0
pelec.chem_integrator = 1
pelec.diffuse_enth = 0
pelec.add_ext_src = 0
pelec.clean_massfrac=0

# TIME STEP CONTROL
pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt
pelec.cfl            = 0.15     # cfl number for hyperbolic system
#pelec.fixed_dt       = 2e-7
pelec.init_shrink    = 1.0     # scale back initial timestep
pelec.change_max     = 1.05    # scale back initial timestep

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in Castro.cpp
amr.v                = 1       # verbosity in Amr.cpp
amr.data_log         = datlog

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 4 4 4 4 # how often to regrid
amr.blocking_factor = 8       # block factor in grid generation
amr.max_grid_size   = 64
# amr.loadbalance_with_workestimates = 1

# CHECKPOINT FILES
amr.checkpoint_files_output = 1
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 500        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt        # root name of plotfile
amr.plot_int        = 100       # number of timesteps between plotfiles
# amr.plot_int        = 50       # number of timesteps between plotfiles
amr.derive_plot_vars = ALL
pelec.plot_massfrac = 1

# PROBLEM PARAMETERS
prob.p0 = 1.0e6
prob.T0 = 300
prob.u_swirl_ax = 500.0
prob.T_swirl = 300.0
prob.l_exit = 30.0
prob.l_combustor = 10.0

# TAGGING
tagging.denerr = 1e20
tagging.dengrad = 8e-3
tagging.max_denerr_lev = 3
tagging.max_dengrad_lev = 3
tagging.max_vfracerr_lev = 1

eb2.geom_type = "converging-nozzle"
eb2.small_volfrac=1.e-6
ebd.boundary_grad_stencil_type = 0
pelec.eb_boundary_T = 300.
pelec.eb_isothermal = 0
pelec.eb_noslip = 0
