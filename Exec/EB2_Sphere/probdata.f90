module probdata_module

!     Combustor variables
      double precision, save :: pamb, Tinit, uinit
      double precision, save :: Tin
      double precision, save :: Rfrac
      double precision, save :: swrlang
      double precision, save :: pertmag


      double precision, save :: uin


      double precision :: rho_l, rhoe_l

!     These help specify which specific problem
      integer        , save ::  probtype,idir

      double precision, save :: split(3)
      
end module probdata_module
