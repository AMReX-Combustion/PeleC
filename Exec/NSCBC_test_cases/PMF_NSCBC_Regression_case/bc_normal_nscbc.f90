  subroutine bcnormal_nscbc(x,u_int,u_ext,dir,sgn,bc_type,bc_params,rho_only)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, NVAR
    use meth_params_module, only: nb_nscbc_params
    use network, only: nspec, naux
    use prob_params_module, only : Interior, Inflow, Outflow, SlipWall, NoSlipWall, &
                                   problo, probhi
 
    
    
    use bl_constants_module, only: M_PI
    
    implicit none

    double precision :: x(3)
    double precision :: u_int(NVAR),u_ext(NVAR)
    logical rho_only
    integer :: dir,sgn, nPMF
    integer, intent(out) :: bc_type
    double precision, intent(out) :: bc_params(nb_nscbc_params)

    type (eos_t) :: eos_state
    double precision :: u(3), rho_inv
    double precision :: relax_U, relax_V, relax_T, beta, sigma_out
    double precision, allocatable :: pmf_vals(:)

    call build(eos_state)
    allocate(pmf_vals(nspec+4))
!    if (.not. pmf_initialized) then
!       call init_pmf()
!    end if



    if (rho_only .EQV. .TRUE. ) then
       if (dir.eq.prob_axis .and. sgn.eq.+1) then
          u_ext(1) = fuel_state(URHO)
       endif
    else
       if (sgn.eq.-1) then


          ! Set outflow pressure
           bc_type = Outflow
           sigma_out = 0.25d0
           bc_params(1) = sigma_out

!           call pmf(probhi(1),probhi(1),pmf_vals,nPMF)

!           eos_state % molefrac(1:nspec) = pmf_vals(4:3+nspec)
!           eos_state % T = pmf_vals(1)

eos_state % molefrac(1) = 1.75291091784e-09 
eos_state % molefrac(2) = 3.20874455466e-33 
eos_state % molefrac(3) = 1.84876115857e-20 
eos_state % molefrac(4) = 0.195621797515 
eos_state % molefrac(5) = 1.38994428424e-23 
eos_state % molefrac(6) = 2.27816535353e-15 
eos_state % molefrac(7) = 3.94814815092e-19 
eos_state % molefrac(8) = 5.24062291401e-19 
eos_state % molefrac(9) = 3.99849375957e-55 
eos_state % molefrac(10) = 0 
eos_state % molefrac(11) = 1.85356538075e-41 
eos_state % molefrac(12) = 6.30833028933e-43 
eos_state % molefrac(13) = 2.34400713209e-20 
eos_state % molefrac(14) = 0.0684676291296 
eos_state % molefrac(15) = 6.0981134532e-17 
eos_state % molefrac(16) = 1.1094324385e-19 
eos_state % molefrac(17) = 3.17519684053e-41 
eos_state % molefrac(18) = 8.05466547043e-21 
eos_state % molefrac(19) = 9.248240059289999e-44 
eos_state % molefrac(20) = 3.96678395197e-32 
eos_state % molefrac(21) = 1.07250785391e-21 
eos_state % molefrac(22) = 2.20913371338e-55 
eos_state % molefrac(23) = 6.87546428076e-22 
eos_state % molefrac(24) = 2.07174123823e-46 
eos_state % molefrac(25) = 2.83642180682e-21 
eos_state % molefrac(26) = 0 
eos_state % molefrac(27) = 1.19367255742e-21 
eos_state % molefrac(28) = 4.25157392917e-42 
eos_state % molefrac(29) = 2.72511545355e-22 
eos_state % molefrac(30) = 2.70531797932e-23 
eos_state % molefrac(31) = 5.42811126531e-28 
eos_state % molefrac(32) = 4.68630667768e-36 
eos_state % molefrac(33) = 1.6705613424e-20 
eos_state % molefrac(34) = 5.88665298123e-22 
eos_state % molefrac(35) = 3.39277716248e-43 
eos_state % molefrac(36) = 8.4121111362e-22 
eos_state % molefrac(37) = 3.99083193418e-22 
eos_state % molefrac(38) = 2.1251891218e-22 
eos_state % molefrac(39) = 2.63395838447e-22 
eos_state % molefrac(40) = 3.82001624801e-52 
eos_state % molefrac(41) = 3.44997311808e-22 
eos_state % molefrac(42) = 1.36116756229e-24 
eos_state % molefrac(43) = 0 
eos_state % molefrac(44) = 1.98984120784e-22 
eos_state % molefrac(45) = 3.65047146931e-24 
eos_state % molefrac(46) = 2.12432804424e-22 
eos_state % molefrac(47) = 9.33069927205e-24 
eos_state % molefrac(48) = 0.735910571602 
eos_state % molefrac(49) = 0 
eos_state % molefrac(50) = 3.02228185193e-26 
eos_state % molefrac(51) = 2.32179619751e-22 
eos_state % molefrac(52) = 1.28451850281e-39 
eos_state % molefrac(53) = 2.45409055206e-22 


          eos_state % T = 298 !pmf_vals(1)

          eos_state % p = pamb
          u = 0
          u(prob_axis) = 18.7936517793 !pmf_vals(2)

           call eos_xty(eos_state) ! get mass fractions from mole fractions
           call eos_tp(eos_state)

           u_ext(URHO ) = eos_state % rho
           u_ext(UMX  ) = eos_state % rho  *  u(1)
           u_ext(UMY  ) = eos_state % rho  *  u(2)
           u_ext(UMZ  ) = eos_state % rho  *  u(3)
           u_ext(UEINT) = eos_state % rho  *  eos_state % e
           u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
           u_ext(UTEMP) = eos_state % T
           u_ext(UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)

       else

          relax_U = 1.0d0
          relax_T = - 1.0d0

          bc_type = Inflow
          bc_params(1) = relax_U
          bc_params(2) = relax_V
          bc_params(3) = relax_T

!          call pmf(problo(1),problo(1),pmf_vals,nPMF)

!          eos_state % molefrac(1:nspec) = pmf_vals(4:3+nspec)

eos_state % molefrac(1) = 1.75291091784e-09
eos_state % molefrac(2) = 3.20874455466e-33
eos_state % molefrac(3) = 1.84876115857e-20
eos_state % molefrac(4) = 0.195621797515
eos_state % molefrac(5) = 1.38994428424e-23
eos_state % molefrac(6) = 2.27816535353e-15
eos_state % molefrac(7) = 3.94814815092e-19
eos_state % molefrac(8) = 5.24062291401e-19
eos_state % molefrac(9) = 3.99849375957e-55
eos_state % molefrac(10) = 0
eos_state % molefrac(11) = 1.85356538075e-41
eos_state % molefrac(12) = 6.30833028933e-43
eos_state % molefrac(13) = 2.34400713209e-20
eos_state % molefrac(14) = 0.0684676291296
eos_state % molefrac(15) = 6.0981134532e-17
eos_state % molefrac(16) = 1.1094324385e-19
eos_state % molefrac(17) = 3.17519684053e-41
eos_state % molefrac(18) = 8.05466547043e-21
eos_state % molefrac(19) = 9.248240059289999e-44
eos_state % molefrac(20) = 3.96678395197e-32
eos_state % molefrac(21) = 1.07250785391e-21
eos_state % molefrac(22) = 2.20913371338e-55
eos_state % molefrac(23) = 6.87546428076e-22
eos_state % molefrac(24) = 2.07174123823e-46
eos_state % molefrac(25) = 2.83642180682e-21
eos_state % molefrac(26) = 0
eos_state % molefrac(27) = 1.19367255742e-21
eos_state % molefrac(28) = 4.25157392917e-42
eos_state % molefrac(29) = 2.72511545355e-22
eos_state % molefrac(30) = 2.70531797932e-23
eos_state % molefrac(31) = 5.42811126531e-28
eos_state % molefrac(32) = 4.68630667768e-36
eos_state % molefrac(33) = 1.6705613424e-20
eos_state % molefrac(34) = 5.88665298123e-22
eos_state % molefrac(35) = 3.39277716248e-43
eos_state % molefrac(36) = 8.4121111362e-22
eos_state % molefrac(37) = 3.99083193418e-22
eos_state % molefrac(38) = 2.1251891218e-22
eos_state % molefrac(39) = 2.63395838447e-22
eos_state % molefrac(40) = 3.82001624801e-52
eos_state % molefrac(41) = 3.44997311808e-22
eos_state % molefrac(42) = 1.36116756229e-24
eos_state % molefrac(43) = 0
eos_state % molefrac(44) = 1.98984120784e-22
eos_state % molefrac(45) = 3.65047146931e-24
eos_state % molefrac(46) = 2.12432804424e-22
eos_state % molefrac(47) = 9.33069927205e-24
eos_state % molefrac(48) = 0.735910571602
eos_state % molefrac(49) = 0
eos_state % molefrac(50) = 3.02228185193e-26
eos_state % molefrac(51) = 2.32179619751e-22
eos_state % molefrac(52) = 1.28451850281e-39
eos_state % molefrac(53) = 2.45409055206e-22


          eos_state % T = 298 !pmf_vals(1)

          eos_state % p = pamb
          u = 0
          u(prob_axis) = 18.7936517793 !pmf_vals(2)

          call eos_xty(eos_state) ! get mass fractions from mole fractions
          call eos_tp(eos_state)

          u_ext(URHO ) = eos_state % rho
          u_ext(UMX  ) = eos_state % rho  *  u(1)
          u_ext(UMY  ) = eos_state % rho  *  u(2)
          u_ext(UMZ  ) = eos_state % rho  *  u(3)
          u_ext(UEINT) = eos_state % rho  *  eos_state % e
          u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
          u_ext(UTEMP) = eos_state % T
          u_ext(UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)

       endif
    endif

    deallocate(pmf_vals)
    call destroy(eos_state)

  end subroutine bcnormal_nscbc
