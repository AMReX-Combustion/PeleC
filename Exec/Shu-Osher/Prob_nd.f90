module pc_prob_module

    implicit none

    private

    public :: amrex_probinit, pc_initdata, pc_prob_close

contains

    subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

        use amrex_fort_module, only : amrex_real
        use bl_error_module 
        use network
        use probdata_module

        implicit none

        integer :: init, namlen
        integer :: name(namlen)
        real(amrex_real) :: problo(3), probhi(3)

        integer :: untin, i

        namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r_base, rho_r_amp, rho_r_osc

        !
        !     Build "probin" filename -- the name of file containing fortin namelist.
        !     
        integer maxlen
        parameter (maxlen=256)
        character probin*(maxlen)

        if (namlen .gt. maxlen) then
            call bl_error("probin file name too long")
        end if

        do i = 1, namlen
        probin(i:i) = char(name(i))
        end do

        ! set namelist defaults
        rho_l       = 3.857143d0      ! left density  (g / cc)
        u_l         = 2.629369d0      ! left velocity (cm / s)
        p_l         = 10.33333d0      ! left pressure (erg / cc)

        rho_r_base   = 1.d0
        rho_r_amp    = 0.2d0
        rho_r_osc    = 5.d0
        u_r          = 0.d0
        p_r          = 1.d0

        !     Read namelists
        untin = 9
        open(untin,file=probin(1:namlen),form='formatted',status='old')
        read(untin,fortin)
        close(unit=untin)
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
    ! ::: nvar      => number of state components.
    ! ::: state     <= scalar array
    ! ::: delta     => cell size
    ! ::: xlo, xhi  => physical locations of lower left and upper
    ! :::              right hand corner of grid.  (does not include
    ! :::		   ghost region).
    ! ::: -----------------------------------------------------------

    subroutine pc_initdata(level,time,lo,hi,nvar, &
            state,state_lo,state_hi, &
            delta,xlo,xhi) bind(C, name = "pc_initdata")

        use amrex_fort_module, only : amrex_real
        use network, only : nspec
        use eos_module, only : eos_rp
        use eos_type_module, only : eos_t, build
        use probdata_module
        use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS

        implicit none

        integer :: level, nvar
        integer :: lo(3), hi(3)
        integer :: state_lo(3), state_hi(3)
        real(amrex_real) :: xlo(3), xhi(3), time, delta(3)
        real(amrex_real) :: state(state_lo(1):state_hi(1), &
            state_lo(2):state_hi(2), &
            state_lo(3):state_hi(3),nvar)
        real(amrex_real) :: cen(3), cdist
        real(amrex_real) :: u
        type(eos_t)      :: eos_state
        real(amrex_real) :: massfrac(nspec)
        integer          :: i,j,k

        call build(eos_state)

        massfrac(:) = 0.d0
        massfrac(1) = 1.d0

        do k = lo(3), hi(3)

            cen(3) = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)

            do j = lo(2), hi(2)

                cen(2) = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

                do i = lo(1), hi(1)

                    cen(1) = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

                    if (cen(1) < 1.d0) then      ! left side
                        u               = u_l
                        eos_state%rho   = rho_l
                        eos_state%p     = p_l
                        eos_state%T     = 100000.d0  ! initial guess
                    else                      ! right side
                        u   = u_r
                        eos_state%rho   = rho_r_base + &
                            rho_r_amp * sin(rho_r_osc * cen(1))
                        eos_state%p     = p_r
                        eos_state%T     = 100000.d0  ! initial guess
                    end if

                    eos_state%massfrac(:) = massfrac(:)

                    call eos_rp(eos_state)

                    state(i,j,k,URHO)   = eos_state%rho
                    state(i,j,k,UMX)    = eos_state%rho * u
                    state(i,j,k,UMY)    = 0.d0
                    state(i,j,k,UMZ)    = 0.d0
                    state(i,j,k,UEDEN)  = eos_state%rho * eos_state%e + &
                        0.5 * eos_state%rho * u * u
                    state(i,j,k,UEINT)  = eos_state%rho * eos_state%e
                    state(i,j,k,UTEMP)  = eos_state%T

                    state(i,j,k,UFS:UFS-1+nspec)    = 0.0d0
                    state(i,j,k,UFS  )              = state(i,j,k,URHO)

                enddo

            enddo

        enddo

    end subroutine pc_initdata

    subroutine pc_prob_close() &
            bind(C, name="pc_prob_close")
    end subroutine pc_prob_close

end module pc_prob_module
