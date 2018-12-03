module probdata_module

  ! HIT parameters
  logical, save            :: restart
  double precision, save   :: lambda0, reynolds_lambda0, mach_t0, prandtl
  double precision, save   :: Lx, Ly, Lz
  double precision, save   :: Linput
  double precision, save   :: k0, rho0, urms0, tau, p0, T0, eint0

  integer, save :: spectrum_type, mode_start, nmodes
  double precision, save   :: turb_scale, force_scale, forcing_time_scale_min, forcing_time_scale_max
  double precision, save :: time_offset

  integer, save :: div_free_force

  integer, parameter :: blrandseed = 0
  integer, parameter :: moderate_zero_modes = 0
  double precision, parameter :: forcing_epsilon = 1.0d-4

  double precision, save, dimension(:,:,:), allocatable :: FTX, FTY, FTZ, TAT, TAP, &
                                                FPX, FPY, FPZ, FAX, FAY, FAZ, &
                                                FPXX, FPYX, FPZX, FPXY, FPYY, FPZY, FPXZ, FPYZ, FPZZ

end module probdata_module
