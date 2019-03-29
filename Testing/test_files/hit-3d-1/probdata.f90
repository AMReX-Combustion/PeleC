module probdata_module

  ! HIT parameters
  character(len=255), save :: iname
  logical, save            :: binfmt
  logical, save            :: restart
  double precision, save   :: lambda0, reynolds_lambda0, mach_t0, prandtl
  integer, save            :: inres
  double precision, save   :: uin_norm
  double precision, save   :: L_x, L_y, L_z
  double precision, save   :: Linput
  double precision, save, dimension(:,:,:), allocatable :: xinput, yinput, zinput, &
       uinput, vinput, winput
  double precision, save, dimension(:), allocatable :: xarray, xdiff
  double precision, save   :: k0, rho0, urms0, tau, p0, T0, eint0

end module probdata_module
