module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use probdata_module
    use amrex_fort_module
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ reynolds, mach, prandtl, rho_x_fact, rho_y_fact, rho_z_fact, u_0_fact, v_0_fact, w_0_fact, u_r_fact, v_r_fact, w_r_fact, p_r_fact, a_rhox, a_rhoy, a_rhoz, a_ur, a_vr, a_wr, a_pr

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    reynolds = 1.d0
    mach = 1.d0
    prandtl = 1.d0
    rho_x_fact = 0.1d0
    rho_y_fact = 1.d0
    rho_z_fact = 1.d0
    u_0_fact = 0.0d0
    v_0_fact = 0.0d0
    w_0_fact = 0.0d0
    u_r_fact = 0.1d0
    v_r_fact = 0.1d0
    w_r_fact = 0.1d0
    p_r_fact = 0.0d0
    a_rhox = 2.0d0
    a_rhoy = 2.0d0
    a_rhoz = 2.0d0
    a_ur = 2.0d0
    a_vr = 2.0d0
    a_wr = 2.0d0
    a_pr = 0.0d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    L_x = probhi(1) - problo(1)
    L_y = probhi(2) - problo(2)
    L_z = probhi(3) - problo(3)

  end subroutine amrex_probinit


  ! ::: -----------------------------------------------------------
  ! ::: This routine is called at problem setup time and is used
  ! ::: to initialize data on each grid.
  ! :::
  ! ::: NOTE:  all arrays have one cell of ghost zones surrounding
  ! :::        the grid interior.  Values in these cells need not
  ! :::        be set here.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: level     => amr level of grid
  ! ::: time      => time at which to init data
  ! ::: lo,hi     => index limits of grid interior (cell centered)
  ! ::: nstate    => number of state components.  You should know
  ! :::		   this already!
  ! ::: state     <=  Scalar array
  ! ::: delta     => cell size
  ! ::: xlo,xhi   => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name="pc_initdata")

    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
    use probdata_module
    use network, only: nspec, naux, molec_wt
    use fundamental_constants_module, only: k_B, n_A
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP, small_temp
    use prob_params_module, only: dim
    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use eos_module

#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    integer :: i,j,k
    double precision :: x,y,z,rho,u,v,w,p,eint
    double precision :: rho0, u0, p0, T0
    type(eos_t) :: eos_state

#ifdef USE_MASA

    ! Define the molecular weight for air
    molec_wt = 28.97d0

    ! Initial pressure and temperature
    p0 = 1.013d6 ! [erg cm^-3]
    T0 = 300.d0

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0  = eos_state % rho
    u0 = mach * eos_state % cs
    const_bulk_viscosity = rho0 * u0 * L_x / reynolds
    const_diffusivity = ZERO
    const_viscosity = rho0 * u0 * L_x / reynolds
    const_conductivity = const_viscosity * eos_state % cp / prandtl

    ! MASA parameters for the following functions
    ! rho = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L);
    ! u = u_0 + u_r * cos(a_ux * PI * r / L);
    ! v = v_0 + v_r * cos(a_vx * PI * r / L);
    ! w = w_0 + w_r * cos(a_wx * PI * r / L);
    ! p = p_0 + p_r * cos(a_px * PI * r / L);
    call masa_set_param("L", L_x)
    call masa_set_param("R", (k_B*n_A)/eos_state % wbar)
    call masa_set_param("k", const_conductivity)
    call masa_set_param("Gamma", gamma_const)
    call masa_set_param("mu", const_viscosity)
    call masa_set_param("mu_bulk", const_bulk_viscosity)
    call masa_set_param("rho_0",rho0)
    call masa_set_param("rho_x",rho_x_fact*rho0)
    call masa_set_param("rho_y",rho_y_fact)
    call masa_set_param("rho_z",rho_z_fact)
    call masa_set_param("u_0",u_0_fact*u0)
    call masa_set_param("v_0",v_0_fact*u0)
    call masa_set_param("w_0",w_0_fact*u0)
    call masa_set_param("u_r",u_r_fact*u0)
    call masa_set_param("v_r",v_r_fact*u0)
    call masa_set_param("w_r",w_r_fact*u0)
    call masa_set_param("p_0",p0)
    call masa_set_param("p_r",p_r_fact*p0)
    call masa_set_param("a_rhox",a_rhox)
    call masa_set_param("a_rhoy",a_rhoy)
    call masa_set_param("a_rhoz",a_rhoz)
    call masa_set_param("a_ur",a_ur)
    call masa_set_param("a_vr",a_vr)
    call masa_set_param("a_wr",a_wr)
    call masa_set_param("a_pr",a_pr)
    call masa_set_param("Cs",ZERO)
    call masa_set_param("CI",ZERO)
    call masa_set_param("PrT",ONE)
    call masa_set_param("deltabar",delta(1))

    ! Display and check
    if ( amrex_pd_ioprocessor() ) then
       call masa_display_param()
    endif
    call masa_sanity_check()

    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             rho = masa_eval_3d_exact_rho(x,y,z)
             u = masa_eval_3d_exact_u(x,y,z)
             v = masa_eval_3d_exact_v(x,y,z)
             w = masa_eval_3d_exact_w(x,y,z)
             p = masa_eval_3d_exact_p(x,y,z)

             ! Call EOS by specifying the density and pressure
             eos_state % rho = rho
             eos_state % p = p
             eos_state % massfrac    = 0.d0
             eos_state % massfrac(1) = 1.d0
             call eos_rp(eos_state)
             eint = eos_state % e

             ! Fill the states
             state(i,j,k,URHO)            = rho
             state(i,j,k,UFS:UFS+nspec-1) = rho * eos_state % massfrac(1:nspec)
             state(i,j,k,UMX)             = rho * u
             if (dim .ge. 2) then
                state(i,j,k,UMY)             = rho * v
                if (dim .ge. 3) then
                   state(i,j,k,UMZ)             = rho * w
                endif
             endif
             state(i,j,k,UTEMP)           = eos_state % T
             state(i,j,k,UEINT)           = rho * eint
             state(i,j,k,UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))

          enddo
       enddo
    enddo

#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif
  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
