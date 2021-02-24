# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10000
#max_step = 2
stop_time = 3e-2

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic =  0   0  0
geometry.coord_sys   =  0       # 0 => cart
geometry.prob_lo     =  0    0    0 

geometry.prob_lo     =   0.0 -0.05 -0.05
geometry.prob_hi     =   0.1 0.05  0.05
amr.n_cell           =   16     16    16

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<

pelec.lo_bc       =  "FOExtrap" "FOExtrap" "FOExtrap"
pelec.hi_bc       =  "FOExtrap" "FOExtrap" "FOExtrap"

# Problem setup
pelec.eb_boundary_T = 1000.
pelec.eb_isothermal = 1
pelec.eb_small_vfrac = 0.0

# WHICH PHYSICS
pelec.do_hydro = 0
pelec.do_mol = 1
pelec.do_react = 0
pelec.allow_negative_energy = 0
pelec.diffuse_temp = 1
pelec.diffuse_vel  = 0
pelec.diffuse_spec = 0
pelec.diffuse_enth = 1
transport.const_conductivity = 2.7271624e+04

# TIME STEP CONTROL
pelec.dt_cutoff      = 5.e-20   # level 0 timestep below which we halt
pelec.cfl            = 0.2     # cfl number for hyperbolic system
# pelec.cfl            = 0.1     # cfl number for hyperbolic system
pelec.init_shrink    = 0.8      # scale back initial timestep
pelec.change_max     = 1.05     # maximum increase in dt over successive steps
pelec.max_dt         = 1e-5

# DIAGNOSTICS & VERBOSITY
pelec.sum_interval   = 1       # timesteps between computing mass
pelec.v              = 1       # verbosity in PeleC cpp files
amr.v                = 1       # verbosity in Amr.cpp
#amr.grid_log         = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2       # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 4

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk      # root name of checkpoint file
amr.check_int       = -1       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 1
amr.plot_file       = plt
amr.plot_int        = 100
#amr.plot_int        = 2
amr.derive_plot_vars=ALL

#extruded triangles lets the user create a maximum of 5 triangles 
#in 2D that will be extruded in the z direction
#make sure the coordinates are in anti-clockwise direction 
eb2.geom_type = "plane"
eb2.plane_point = 0.053125 0 0
eb2.plane_normal = 1.0 -1.0 0
# eb2.geom_type = "sphere"
# eb2.sphere_center = 0.05 0 0
# eb2.sphere_radius = 0.025
# eb2.sphere_has_fluid_inside = 1

ebd.boundary_grad_stencil_type = 0
#fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM PARAMETERS
prob.p_init    = 1013250.0
prob.T_init    = 900.0
