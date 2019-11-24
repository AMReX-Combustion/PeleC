module reactions_module

  implicit none

  public

contains

  subroutine pc_react_state(lo,hi, &
                            uold,uo_lo,uo_hi, &
                            unew,un_lo,un_hi, &
                            asrc,as_lo,as_hi, &
                            mask,m_lo,m_hi, &
                            cost,c_lo,c_hi, &
                            IR,IR_lo,IR_hi, &
#ifdef PELEC_USE_EB
                            flag, fglo, fghi, &
#endif
                            time,dt_react,do_update) bind(C, name="pc_react_state")

#ifdef PELEC_USE_EB
use amrex_ebcellflag_module, only : is_covered_cell
#endif

    use network, only : nspecies
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS
#ifdef USE_SUNDIALS_PP
    use cvode_module, only : react
#else
    use reactor_module, only : react
#endif
    use amrex_fort_module, only : amrex_real
    use amrex_constants_module, only : HALF

    implicit none

    integer          ::    lo(3),    hi(3)
    integer          :: uo_lo(3), uo_hi(3)
    integer          :: un_lo(3), un_hi(3)
    integer          :: as_lo(3), as_hi(3)
    integer          ::  m_lo(3),  m_hi(3)
    integer          ::  c_lo(3),  c_hi(3)
    integer          :: IR_lo(3), IR_hi(3)
    double precision :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    double precision :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)
    double precision :: asrc(as_lo(1):as_hi(1),as_lo(2):as_hi(2),as_lo(3):as_hi(3),NVAR)
    integer          :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    double precision :: cost(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision :: IR(IR_lo(1):IR_hi(1),IR_lo(2):IR_hi(2),IR_lo(3):IR_hi(3),nspecies+1)
    double precision :: time, dt_react
    integer          :: do_update
#ifdef PELEC_USE_EB
    integer, intent(in) :: fglo(3),fghi(3)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
#endif
    integer          :: i, j, k
    double precision :: rho_e_K_old,rho_e_K_new, rhoE_old, rhoE_new, rho_new, mom_new(3)

    real(amrex_real) ::    rY(nspecies+1), rY_src(nspecies)
    real(amrex_real) ::    energy, energy_src, pressure, rho
#ifdef USE_SUNDIALS_PP
    real(amrex_real) ::    nrg(1), nrg_src(1)
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#ifdef PELEC_USE_EB
             if (mask(i,j,k) .eq. 1 .and. .not. is_covered_cell(flag(i,j,k))) then
#else
             if (mask(i,j,k) .eq. 1) then
#endif

                rhoE_old                        = uold(i,j,k,UEDEN)
                rho_e_K_old                     = HALF * sum(uold(i,j,k,UMX:UMZ)**2) / uold(i,j,k,URHO)
                rho                             = sum(uold(i,j,k,UFS:UFS+nspecies-1))

                energy           = (rhoE_old - rho_e_K_old) / uold(i,j,k,URHO)
                rY(nspecies+1)      = uold(i,j,k,UTEMP)
                rY(1:nspecies)      = uold(i,j,k,UFS:UFS+nspecies-1)

                ! rho.e source term computed using (rho.E,rho.u,rho)_new rather than pulling from UEINT comp of asrc
                rho_e_K_new = HALF * sum(unew(i,j,k,UMX:UMZ)**2) / unew(i,j,k,URHO)
                energy_src       = ( (unew(i,j,k,UEDEN) - rho_e_K_new) &
                       -                (rho  *  energy) ) / dt_react

#ifdef USE_SUNDIALS_PP
                nrg(1)           = rho  *  energy
                nrg_src(1)       = energy_src
#endif

                rY_src(1:nspecies)  = asrc(i,j,k,UFS:UFS+nspecies-1)
                !react_state_in % i = i
                !react_state_in % j = j
                !react_state_in % k = k

                pressure         = 1013250.d0
                
#ifdef USE_SUNDIALS_PP
                cost(i,j,k) = react(rY, rY_src, nrg, nrg_src,&
#else
                cost(i,j,k) = react(rY, rY_src, energy, energy_src,&
                                    pressure,&
#endif
                                    dt_react,time)


                rho_new = sum(rY(1:nspecies))
                mom_new = uold(i,j,k,UMX:UMZ) + dt_react*asrc(i,j,k,UMX:UMZ)
#ifdef USE_SUNDIALS_PP
                rhoE_new = nrg(1)
#else
                rhoE_new = rho_new  *  energy
#endif
                rho_e_K_new = HALF * sum(mom_new**2) / rho_new
                rhoE_new    = rhoE_new + rho_e_K_new
                
                if (do_update .eq. 1) then

                   unew(i,j,k,URHO)            = rho_new
                   unew(i,j,k,UMX:UMZ)         = mom_new
#ifdef USE_SUNDIALS_PP
                   unew(i,j,k,UEINT)           = nrg(1)
#else
                   unew(i,j,k,UEINT)           = rho_new*energy
#endif
                   unew(i,j,k,UEDEN)           = rhoE_new
                   unew(i,j,k,UTEMP)           = rY(nspecies+1)
                   unew(i,j,k,UFS:UFS+nspecies-1) = rY(1:nspecies)

                endif

                ! Add drhoY/dt to reactions MultiFab, but be
                ! careful because the reactions and state MFs may
                ! not have the same number of ghost cells.
                ! Also, be careful because the container for uold might be the same as that for unew
                if ( i .ge. IR_lo(1) .and. i .le. IR_hi(1) .and. &
                     j .ge. IR_lo(2) .and. j .le. IR_hi(2) .and. &
                     k .ge. IR_lo(3) .and. k .le. IR_hi(3) ) then

                   IR(i,j,k,1:nspecies) = (rY(1:nspecies) - uold(i,j,k,UFS:UFS+nspecies-1)) / dt_react - asrc(i,j,k,UFS:UFS+nspecies-1)
                   IR(i,j,k,nspecies+1) = (rhoE_new          -         rhoE_old     ) / dt_react - asrc(i,j,k,UEDEN          )

                endif

             end if
          end do
       enddo
    enddo


  end subroutine pc_react_state

#ifdef REACTIONS
  subroutine pc_react_state_expl(lo,hi, &
                            uold,uo_lo,uo_hi, &
                            unew,un_lo,un_hi, &
                            asrc,as_lo,as_hi, &
                            mask,m_lo,m_hi, &
                            cost,c_lo,c_hi, &
                            IR,IR_lo,IR_hi, &
                            time,dt_react,do_update,&
                            nsubsteps_min,nsubsteps_max,nsubsteps_guess,errtol) bind(C, name="pc_react_state_expl")

    use eos_type_module
    use eos_module, only : eos_t,eos_rt
    use network, only : nspecies
    use chemistry_module  , only : molecular_weight
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS
    use amrex_fort_module, only : amrex_real
    use amrex_constants_module, only : HALF
    use rk_params_module
    use fuego_chemistry, only : VCKWYR

    implicit none

    integer           ::    lo(3),    hi(3)
    integer           :: uo_lo(3), uo_hi(3)
    integer           :: un_lo(3), un_hi(3)
    integer           :: as_lo(3), as_hi(3)
    integer           ::  m_lo(3),  m_hi(3)
    integer           ::  c_lo(3),  c_hi(3)
    integer           :: IR_lo(3), IR_hi(3)
    double precision  :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    double precision  :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)
    double precision  :: asrc(as_lo(1):as_hi(1),as_lo(2):as_hi(2),as_lo(3):as_hi(3),NVAR)
    integer           :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    double precision  :: cost(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision  :: IR(IR_lo(1):IR_hi(1),IR_lo(2):IR_hi(2),IR_lo(3):IR_hi(3),nspecies+1)
    double precision  :: time, dt_react
    integer           :: do_update
    integer           :: nsubsteps_min,nsubsteps_max,nsubsteps_guess
    double precision  :: errtol

    integer           :: i, j, k

    !Euler explicit for testing
    !integer,parameter :: nrkstages=1
    !double precision,parameter,dimension(nrkstages)  :: rkcoeffs=(/1.d0/)

    integer :: npts,steps,stage,ns

    double precision :: saneval(NVAR)
    double precision :: urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    double precision :: yrk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspecies)
    double precision :: urk_carryover(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    double precision :: urk_err(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    double precision :: wdot(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspecies)
    double precision :: eint(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspecies)
    double precision :: cv(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: rhoedot_ext(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: rhoydot_ext(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspecies)
    double precision :: re_old(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: tempsrc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    type (eos_t) :: eos_state

    double precision :: rhoE_old,rho_e_K_old,rho,energy,rho_e_K_new,rho_new
    double precision :: rhoe_new,rhoet_new,rhoInv,mom_new(3)
    double precision :: dt_rk,dt_rk_max,dt_rk_min,updt_time
    
    call build(eos_state)

    !initialize all local arrays
    yrk         = 0.d0
    urk         = uold(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)
    urk_err     = 0.d0
    wdot        = 0.d0
    eint        = 0.d0
    cv          = 0.d0
    rhoedot_ext = 0.d0
    rhoydot_ext = 0.d0
    re_old      = 0.d0
    tempsrc     = 0.d0
    mom_new     = 0.d0

    !==============================================================
    !sanitize urk so that masked values are sane and not NaN
    !May be removed if we are confident this is not the case
    !with uninitialized fortran arrays I am not sure!
    !==============================================================
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                if(mask(i,j,k) .eq. 1) then
                    saneval(:)=urk(i,j,k,:)
                    exit
                endif
            enddo
        enddo
    enddo
    
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                if(mask(i,j,k) .ne. 1) then
                    urk(i,j,k,:)=saneval(:)
                endif
            enddo
        enddo
    enddo
    !===============================================================

    !setting urk_carryover
    urk_carryover = urk
    
    !compute rhoe_ext/rhoy_ext - one time operation, keeping original 
    !pc_react_state code ===========================================
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                if(mask(i,j,k) .eq. 1) then

                    rhoE_old      = uold(i,j,k,UEDEN)
                    rho_e_K_old   = HALF * sum(uold(i,j,k,UMX:UMZ)**2) / uold(i,j,k,URHO)
                    rho           = sum(uold(i,j,k,UFS:UFS+nspecies-1))

                    energy        = (rhoE_old - rho_e_K_old) / uold(i,j,k,URHO)
                    re_old(i,j,k) = energy*uold(i,j,k,URHO)

                    ! rho.e source term computed using 
                    !(rho.E,rho.u,rho)_new rather than pulling from UEINT comp of asrc
                    rho_e_K_new   = HALF * sum(unew(i,j,k,UMX:UMZ)**2)/unew(i,j,k,URHO)
                    rhoedot_ext(i,j,k)  = ( (unew(i,j,k,UEDEN) - rho_e_K_new) &
                                        -   (rho  *  energy) ) / dt_react
                    rhoydot_ext(i,j,k,1:nspecies)  = asrc(i,j,k,UFS:UFS+nspecies-1)
                endif
            enddo
        enddo
    enddo
    !===============================================================
    dt_rk     = dt_react/nsubsteps_guess
    dt_rk_min = dt_react/nsubsteps_max
    dt_rk_max = dt_react/nsubsteps_min
    !===============================================================

    !write(6,*) "dt_rk_guess:",dt_rk
    !write(6,*) "============================"
    !write(6,*) "============================"
    !flush(6)

    npts=(hi(3)-lo(3)+1)*(hi(2)-lo(2)+1)*(hi(1)-lo(1)+1)
    updt_time=0.d0
    steps=0

    do while(updt_time .lt. dt_react)
    !do steps=1,nsubsteps

        urk_carryover   = urk
        urk_err     = 0.d0
        wdot        = 0.d0
        tempsrc     = 0.d0
        updt_time = updt_time+dt_rk
        steps     = steps+1

        !write(6,*) "updt_time,dt_rk:",updt_time,dt_rk,dt_rk_min,dt_rk_max,dt_react
        !flush(6)

        do stage=1,rk64_stages

           !computing mass fractions
           !hope the compiler can optimize this better
           do k=lo(3),hi(3)
              do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    yrk(i,j,k,1:nspecies)=urk(i,j,k,UFS:UFS+nspecies-1)/urk(i,j,k,URHO)
                enddo
              enddo
           enddo

           !Returns the molar production rate of species
           !Given rho, T, and mass fractions
           call VCKWYR(npts, urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO), &
               urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UTEMP), &
               yrk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nspecies), &
               wdot)

            !Ideally this has to come from a vector call to fuego
            !will replace after pelephysics has this facility
            do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                    do i=lo(1),hi(1)

                            eos_state % rho               = urk(i,j,k,URHO)
                            rhoInv                        = 1.d0 / eos_state % rho
                            eos_state % massfrac(1:nspecies) = urk(i,j,k,UFS:UFS+nspecies-1) * rhoInv
                            eos_state % T                 = urk(i,j,k,UTEMP)
                            call eos_rt(eos_state)
                            do ns=1,nspecies
                                eint(i,j,k,ns)=eos_state%ei(ns)
                            enddo
                            cv(i,j,k)=eos_state%cv

                    enddo
                enddo
            enddo

           !convert wdot to gms/s and add external source
           !hope the compiler can optimize this better
           do k=lo(3),hi(3)
            do j=lo(2),hi(2)
              do i=lo(1),hi(1)
                wdot(i,j,k,1:nspecies) = wdot(i,j,k,1:nspecies)*molecular_weight(1:nspecies) + rhoydot_ext(i,j,k,1:nspecies)
              enddo
            enddo
           enddo
           
           !update temperature src
           !hope the compiler can optimize this better
           do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                    do i=lo(1),hi(1)
                        tempsrc(i,j,k)=rhoedot_ext(i,j,k)
                        do ns=1,nspecies
                            tempsrc(i,j,k)=tempsrc(i,j,k)-wdot(i,j,k,ns)*eint(i,j,k,ns)
                        enddo
                        tempsrc(i,j,k)=tempsrc(i,j,k)/urk(i,j,k,URHO)/cv(i,j,k)
                    enddo
                enddo
            enddo
           
           !update species and temperature
           !hope the compiler can optimize this better
           do k=lo(3),hi(3)
            do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                !=====================
                !update urk_err
                !=====================

                !update species
                urk_err(i,j,k,UFS:UFS+nspecies-1)  = urk_err(i,j,k,UFS:UFS+nspecies-1) + &
                err_rk64(stage)*dt_rk*wdot(i,j,k,1:nspecies)*mask(i,j,k)

                !update temperature 
                urk_err(i,j,k,UTEMP) = urk_err(i,j,k,UTEMP) + err_rk64(stage)*dt_rk*tempsrc(i,j,k)*mask(i,j,k)
                !===============================================================================================

                !=====================
                !update stage solution
                !=====================

                !update species
                urk(i,j,k,UFS:UFS+nspecies-1)  = urk_carryover(i,j,k,UFS:UFS+nspecies-1) + &
                alpha_rk64(stage)*dt_rk*wdot(i,j,k,1:nspecies)*mask(i,j,k)

                !update temperature 
                urk(i,j,k,UTEMP) = urk_carryover(i,j,k,UTEMP) + alpha_rk64(stage)*dt_rk*tempsrc(i,j,k)*mask(i,j,k)
                !===============================================================================================

                !=====================
                !update urk_carryover
                !=====================
                
                !update species
                urk_carryover(i,j,k,UFS:UFS+nspecies-1)  = urk(i,j,k,UFS:UFS+nspecies-1) + &
                beta_rk64(stage)*dt_rk*wdot(i,j,k,1:nspecies)*mask(i,j,k)

                !update temperature 
                urk_carryover(i,j,k,UTEMP) = urk(i,j,k,UTEMP) + beta_rk64(stage)*dt_rk*tempsrc(i,j,k)*mask(i,j,k)
                !===============================================================================================

            enddo
           enddo
          enddo
        
          !update density
            do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                    do i=lo(1),hi(1)
                    !update density
                    urk(i,j,k,URHO)  = sum(urk(i,j,k,UFS:UFS+nspecies-1))
                    enddo
                enddo
            enddo

        enddo !stage loop

        call adapt_timestep(lo,hi,urk_err,dt_rk,dt_rk_max,dt_rk_min,errtol)
    
        !write(6,*)"dt_rk:",dt_rk,dt_rk_min,dt_rk_max
        !write(6,*)"================================"
        !flush(6)
 
    enddo !substep loop
    cost(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=steps

    !write(6,*)"No: of chemistry substeps:",steps
    !write(6,*)"================================"
    !flush(6)
    
    !this is a one time operation, so keeping original pc_react_state code
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)

                if(mask(i,j,k) .eq. 1) then

                    rho_new  = urk(i,j,k,URHO)
                    mom_new  = uold(i,j,k,UMX:UMZ) + dt_react*asrc(i,j,k,UMX:UMZ)
                    rhoe_new = re_old(i,j,k) + dt_react*rhoedot_ext(i,j,k)
                    rho_e_K_new = HALF * sum(mom_new**2) / rho_new
                    rhoet_new    = rhoe_new + rho_e_K_new

    
                    if (do_update .eq. 1) then

                        unew(i,j,k,URHO)            = rho_new
                        unew(i,j,k,UMX:UMZ)         = mom_new
                        unew(i,j,k,UEINT)           = rhoe_new
                        unew(i,j,k,UEDEN)           = rhoet_new
                        unew(i,j,k,UTEMP)           = urk(i,j,k,UTEMP)
                        unew(i,j,k,UFS:UFS+nspecies-1) = urk(i,j,k,UFS:UFS+nspecies-1)

                    endif

                    ! Add drhoY/dt to reactions MultiFab, but be
                    ! careful because the reactions and state MFs may
                    ! not have the same number of ghost cells.
                    ! Also, be careful because the container for uold might be the same as that for unew
                    if ( i .ge. IR_lo(1) .and. i .le. IR_hi(1) .and. &
                        j .ge. IR_lo(2) .and. j .le. IR_hi(2) .and. &
                        k .ge. IR_lo(3) .and. k .le. IR_hi(3) ) then

                        IR(i,j,k,1:nspecies) = (urk(i,j,k,UFS:UFS+nspecies-1) - uold(i,j,k,UFS:UFS+nspecies-1)) / dt_react -&
                            asrc(i,j,k,UFS:UFS+nspecies-1)
                        IR(i,j,k,nspecies+1) =(rhoet_new - uold(i,j,k,UEDEN)) / dt_react - asrc(i,j,k,UEDEN)

                    endif

                end if
            enddo
        enddo
    enddo
    

  end subroutine pc_react_state_expl

  subroutine adapt_timestep(lo,hi,urk_err,dt_rk4,dt_rk4_max,dt_rk4_min,tol)

      use meth_params_module, only : NVAR

      implicit none
    
      integer           ::    lo(3), hi(3)
      double precision  :: urk_err(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
      double precision  :: dt_rk4, dt_rk4_max, dt_rk4_min, tol

      integer :: i,j,k
      double precision :: max_err,change_factor
      double precision,parameter :: safety_fac=1e4
      double precision,parameter :: exp1=0.25
      double precision,parameter :: exp2=0.2
      double precision,parameter :: beta=1.d0


      max_err=tiny(max_err)

      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                
                if(maxval(abs(urk_err(i,j,k,:))) .gt. max_err) then
                    max_err=maxval(abs(urk_err(i,j,k,:)))
                endif
            enddo
        enddo
    enddo
    
    !write(6,*)"max_err:",max_err
    !flush(6)

    !chance to increase time step
    if(max_err .lt. tol) then

        !limit max_err,can't be 0
        max_err=max(max_err,tol/safety_fac)
        change_factor=beta*(tol/max_err)**(exp1)
        dt_rk4=min(dt_rk4_max,dt_rk4*change_factor)

    !reduce time step (error is high!)
    else
        change_factor=beta*(tol/max_err)**(exp2)
        dt_rk4=max(dt_rk4_min,dt_rk4*change_factor)
    endif

  end subroutine adapt_timestep
#endif

end module reactions_module
