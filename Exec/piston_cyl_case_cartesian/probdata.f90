module probdata_module

      double precision, save ::  Pres_domain, Temp_domain,Yfuel_domain,Yox_domain,YN2_domain
      double precision, save ::  dens_jet,vel_jet
      double precision, save ::  Yfuel_jet,Yox_jet, YN2_jet
      integer, save :: fuel_id, ox_id, N2_id
      
      !geometric parameters
      double precision, save :: centx,centz,r_circ,r_hole,cone_angle
      integer, save :: nholes

end module probdata_module
