module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  real(rt), save :: rpulse = 0.5d0
  real(rt), save :: rho0   = 1.2d-3
  real(rt), save :: drho0  = 1.2d-4
  real(rt), save :: p0     = 1.01325d6
end module probdata_module
