module probdata_module

!     Sod variables
      double precision, save ::  p_l, rho_l, p_r, rho_r,T_l,T_r,rhoe_l,rhoe_r

      integer        , save ::  idir

      double precision, save :: split(3), frac
      
end module probdata_module
