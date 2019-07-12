module rk_params_module

    !RK4 with 6 stages and error correction
    integer, parameter :: rk64_stages=6 
    double precision,parameter,dimension(rk64_stages) ::      &
        alpha_rk64 = (/+  3296351145737.d0/15110423921029.d0, &
                       +  1879360555526.d0/ 7321162733569.d0, &
                       + 10797097731880.d0/20472212111779.d0, &
                       +   754636544611.d0/15563872110659.d0, &
                       +  3260218886217.d0/ 2618290685819.d0, &
                       +  5069185909380.d0/12292927838509.d0/)

    double precision, parameter, dimension(rk64_stages) ::    &
        beta_rk64 = (/-  1204558336989.d0/10607789004752.d0,  &
                      -  3028468927040.d0/14078136890693.d0,  &
                      -   455570672869.d0/ 8930094212428.d0,  &
                      - 17275898420483.d0/15997285579755.d0,  &
                      -  2453906524165.d0/ 9868353053862.d0,  &
                      0.d0/)
      double precision, parameter, dimension(rk64_stages) ::  &
      err_rk64 = (/-   530312978447.d0/ 9560368366154.d0,     &
                   +    473021958881.d0/ 2984707536468.d0,    &
                   -    947229622805.d0/10456009803779.d0,    &
                   -   2921473878215.d0/13334914072261.d0,    &
                   +   1519535112975.d0/ 9264196100452.d0,    &
                   +    167623581683.d0/ 3930932046784.d0/) 


  !$acc declare create(alpha_rk64,beta_rk64,err_rk64)
end module rk_params_module
