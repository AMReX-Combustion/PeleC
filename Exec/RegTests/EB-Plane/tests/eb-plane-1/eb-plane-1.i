# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time = 6.0

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic =  0   0  1
geometry.coord_sys   =  0       # 0 => cart
geometry.prob_lo     =  0    0    0 
geometry.prob_hi     =  7.5  0.9375  0.9375
amr.n_cell           =  64   8   8

# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
# Interior, UserBC, Symmetry, SlipWall, NoSlipWall
# >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<

pelec.lo_bc       =  "Hard"     "FOExtrap" "Interior"
pelec.hi_bc       =  "FOExtrap" "FOExtrap" "Interior"

# Problem setup
pelec.eb_boundary_T = 300.
pelec.eb_isothermal = 1

# WHICH PHYSICS
pelec.do_hydro = 1
pelec.do_mol = 1
pelec.do_react = 0
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
amr.max_level       = 2       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2       # how often to regrid
amr.blocking_factor = 8       # block factor in grid generation
amr.max_grid_size   = 8 

# CHECKPOINT FILES
amr.checkpoint_files_output = 0
amr.check_file      = chk      # root name of checkpoint file
amr.check_int       = -1       # number of timesteps between checkpoints

# PLOTFILES
amr.plot_files_output = 0
amr.plot_file       = plt
amr.plot_int        = 1
amr.derive_plot_vars=ALL

#extruded triangles lets the user create a maximum of 5 triangles 
#in 2D that will be extruded in the z direction
#make sure the coordinates are in anti-clockwise direction 
#eb2.geom_type = "all_regular"
eb2.geom_type = "moving_plane"
extruded_triangles.num_tri = 2
extruded_triangles.tri_0_point_0 = 4.25    3.01  0.0
extruded_triangles.tri_0_point_1 = 5.25    2.01  0.0
extruded_triangles.tri_0_point_2 = 5.25    4.01  0.0

extruded_triangles.tri_1_point_0 = 5.25   2.01  0.0
extruded_triangles.tri_1_point_1 = 6.25   3.01  0.0
extruded_triangles.tri_1_point_2 = 5.25   4.01   0.0
ebd.boundary_grad_stencil_type = 0

prob.p = 1013250.0
prob.rho = 0.00116
prob.vx_in =  0.0
prob.vy_in =  0.0
prob.Re_L = 625.0
prob.Pr = 0.7

#fabarray.mfiter_tile_size = 1024 1024 1024
