module probdata_module

    !     Shu-Osher variables
    double precision, save :: p_l, u_l, rho_l, p_r, u_r, rho_r_base, rho_r_amp, &
        rho_r_osc, frac, use_Tinit

    !     These help specify which specific problem
    integer        , save ::  probtype,idir

    double precision, save :: split

end module probdata_module
