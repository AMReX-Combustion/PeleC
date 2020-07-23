module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use amrex_fort_module
    use fuego_chemistry, only: molecular_weight
    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ reynolds, mach, prandtl, rho_x_fact, rho_y_fact, rho_z_fact, u_0_fact, u_x_fact, u_y_fact, u_z_fact, v_0_fact, v_x_fact, v_y_fact, v_z_fact, w_0_fact, w_x_fact, w_y_fact, w_z_fact, p_x_fact, p_y_fact, p_z_fact, a_rhox, a_rhoy, a_rhoz, a_ux, a_uy, a_uz, a_vx, a_vy, a_vz, a_wx, a_wy, a_wz, a_px, a_py, a_pz

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    ! Define the molecular weight for air
    molecular_weight = 28.97d0

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
    u_0_fact = 1.d0
    u_x_fact = 1.d0
    u_y_fact = 1.d0
    u_z_fact = 1.d0
    v_0_fact = 1.d0
    v_x_fact = 1.d0
    v_y_fact = 1.d0
    v_z_fact = 1.d0
    w_0_fact = 1.d0
    w_x_fact = 1.d0
    w_y_fact = 1.d0
    w_z_fact = 1.d0
    p_x_fact = 0.2d0
    p_y_fact = 1.d0
    p_z_fact = 1.d0
    a_rhox = 2.0d0
    a_rhoy = 4.0d0
    a_rhoz = 6.0d0
    a_ux = 2.0d0
    a_uy = 4.0d0
    a_uz = 6.0d0
    a_vx = 4.0d0
    a_vy = 2.0d0
    a_vz = 6.0d0
    a_wx = 6.0d0
    a_wy = 4.0d0
    a_wz = 2.0d0
    a_px = 6.0d0
    a_py = 2.0d0
    a_pz = 4.0d0

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
    use fuego_chemistry, only: Ru, nspecies, naux, molecular_weight
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP, small_temp, Cs, CI, PrT
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
    double precision :: rho0, u0, v0, w0, p0, T0
    type(eos_t) :: eos_state

#ifdef USE_MASA

    ! Define the molecular weight for air
    molecular_weight = 28.97d0

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
    ! u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * u_z * cos(a_uz * PI * z / L);
    ! v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * v_z * cos(a_vz * PI * z / L);
    ! w = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * w_z * cos(a_wz * PI * z / L);
    ! p = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * p_z * cos(a_pz * PI * z / L);
    call masa_set_param("L", L_x)
    call masa_set_param("R", Ru/eos_state % wbar)
    call masa_set_param("k", const_conductivity)
    call masa_set_param("Gamma", gamma_const)
    call masa_set_param("mu", const_viscosity)
    call masa_set_param("mu_bulk", const_bulk_viscosity)
    call masa_set_param("rho_0",rho0)
    call masa_set_param("rho_x",rho_x_fact*rho0)
    call masa_set_param("rho_y",rho_y_fact)
    call masa_set_param("rho_z",rho_z_fact)
    call masa_set_param("u_0",u_0_fact*u0)
    call masa_set_param("u_x",u_x_fact*u0)
    call masa_set_param("u_y",u_y_fact)
    call masa_set_param("u_z",u_z_fact)
    call masa_set_param("v_0",v_0_fact*u0)
    call masa_set_param("v_x",v_x_fact)
    call masa_set_param("v_y",v_y_fact*u0)
    call masa_set_param("v_z",v_z_fact)
    call masa_set_param("w_0",w_0_fact*u0)
    call masa_set_param("w_x",w_x_fact)
    call masa_set_param("w_y",w_y_fact)
    call masa_set_param("w_z",w_z_fact*u0)
    call masa_set_param("p_0",p0)
    call masa_set_param("p_x",p_x_fact*p0)
    call masa_set_param("p_y",p_y_fact)
    call masa_set_param("p_z",p_z_fact)
    call masa_set_param("a_rhox",a_rhox)
    call masa_set_param("a_rhoy",a_rhoy)
    call masa_set_param("a_rhoz",a_rhoz)
    call masa_set_param("a_ux",a_ux)
    call masa_set_param("a_uy",a_uy)
    call masa_set_param("a_uz",a_uz)
    call masa_set_param("a_vx",a_vx)
    call masa_set_param("a_vy",a_vy)
    call masa_set_param("a_vz",a_vz)
    call masa_set_param("a_wx",a_wx)
    call masa_set_param("a_wy",a_wy)
    call masa_set_param("a_wz",a_wz)
    call masa_set_param("a_px",a_px)
    call masa_set_param("a_py",a_py)
    call masa_set_param("a_pz",a_pz)
    call masa_set_param("Cs",Cs)
    call masa_set_param("CI",CI)
    call masa_set_param("PrT",PrT)
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
             state(i,j,k,UFS:UFS+nspecies-1) = rho * eos_state % massfrac(1:nspecies)
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
