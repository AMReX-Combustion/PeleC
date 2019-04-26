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
                            time,dt_react,do_update) bind(C, name="pc_react_state")

    use network           , only : nspec
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS
    use reactor_module, only : react
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
    double precision :: IR(IR_lo(1):IR_hi(1),IR_lo(2):IR_hi(2),IR_lo(3):IR_hi(3),nspec+1)
    double precision :: time, dt_react
    integer          :: do_update

    integer          :: i, j, k
    double precision :: rho_e_K_old,rho_e_K_new, rhoE_old, rhoE_new, rho_new, mom_new(3)

    real(amrex_real) ::    rY(nspec+1), rY_src(nspec)
    real(amrex_real) ::    energy, energy_src, pressure, rho


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (mask(i,j,k) .eq. 1) then

                rhoE_old                        = uold(i,j,k,UEDEN)
                rho_e_K_old                     = HALF * sum(uold(i,j,k,UMX:UMZ)**2) / uold(i,j,k,URHO)
                rho                             = sum(uold(i,j,k,UFS:UFS+nspec-1))

                energy           = (rhoE_old - rho_e_K_old) / uold(i,j,k,URHO)
                rY(nspec+1)      = uold(i,j,k,UTEMP)
                rY(1:nspec)      = uold(i,j,k,UFS:UFS+nspec-1)

                ! rho.e source term computed using (rho.E,rho.u,rho)_new rather than pulling from UEINT comp of asrc
                rho_e_K_new = HALF * sum(unew(i,j,k,UMX:UMZ)**2) / unew(i,j,k,URHO)
                energy_src       = ( (unew(i,j,k,UEDEN) - rho_e_K_new) &
                       -                (rho  *  energy) ) / dt_react

                rY_src(1:nspec)  = asrc(i,j,k,UFS:UFS+nspec-1)
                !react_state_in % i = i
                !react_state_in % j = j
                !react_state_in % k = k

                pressure         = 1013250.d0

                cost(i,j,k) = react(rY, rY_src,&
                                    energy, energy_src,&
                                    pressure,dt_react,time,0)

                rho_new = sum(rY(1:nspec))
                mom_new = uold(i,j,k,UMX:UMZ) + dt_react*asrc(i,j,k,UMX:UMZ)
                rhoe_new = rho_new  *  energy
                rho_e_K_new = HALF * sum(mom_new**2) / rho_new
                rhoE_new    = rhoe_new + rho_e_K_new
                
                if (do_update .eq. 1) then

                   unew(i,j,k,URHO)            = rho_new
                   unew(i,j,k,UMX:UMZ)         = mom_new
                   unew(i,j,k,UEINT)           = rhoe_new
                   unew(i,j,k,UEDEN)           = rhoE_new
                   unew(i,j,k,UTEMP)           = rY(nspec+1)
                   unew(i,j,k,UFS:UFS+nspec-1) = rY(1:nspec)

                endif

                ! Add drhoY/dt to reactions MultiFab, but be
                ! careful because the reactions and state MFs may
                ! not have the same number of ghost cells.
                ! Also, be careful because the container for uold might be the same as that for unew
                if ( i .ge. IR_lo(1) .and. i .le. IR_hi(1) .and. &
                     j .ge. IR_lo(2) .and. j .le. IR_hi(2) .and. &
                     k .ge. IR_lo(3) .and. k .le. IR_hi(3) ) then

                   IR(i,j,k,1:nspec) = (rY(1:nspec) - uold(i,j,k,UFS:UFS+nspec-1)) / dt_react - asrc(i,j,k,UFS:UFS+nspec-1)
                   IR(i,j,k,nspec+1) = (        rhoE_new          -         rhoE_old        ) / dt_react - asrc(i,j,k,UEDEN          )

                endif

             end if
          end do
       enddo
    enddo

  end subroutine pc_react_state

  subroutine pc_react_state_expl(lo,hi, &
                            uold,uo_lo,uo_hi, &
                            unew,un_lo,un_hi, &
                            asrc,as_lo,as_hi, &
                            mask,m_lo,m_hi, &
                            cost,c_lo,c_hi, &
                            IR,IR_lo,IR_hi, &
                            time,dt_react,do_update,nsubsteps) bind(C, name="pc_react_state_expl")

    use network           , only : nspec
    use chemistry_module  , only : molecular_weight
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS
    use reactor_module, only : react
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
    double precision :: IR(IR_lo(1):IR_hi(1),IR_lo(2):IR_hi(2),IR_lo(3):IR_hi(3),nspec+1)
    double precision :: time, dt_react
    integer          :: do_update,nsubsteps

    integer          :: i, j, k
    integer,parameter :: nrkstages=2
    double precision,parameter,dimension(nrkstages)  :: rkcoeffs=(/0.5d0,1.d0/)
    integer :: npts,nsubsteps,dt_rk,steps,stage,ierr,updt_time

    double precision :: urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    double precision :: urk_old(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    double precision :: wdot(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec)
    double precision :: eint(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec)
    double precision :: cv(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: mom_new(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3)
    double precision :: rhoedot_ext(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: re_old(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    type (eos_t) :: eos_state

    double precision :: rhoE_old,rho_e_K_old,rho,energy,rho_e_K_new

    !initialize all local arrays
    urk    =unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    urk_old=unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)
    wdot=0.d0
    eint=0.d0
    cv=0.d0
    mom_new=0.d0
    rhoedot_ext=0.d0
    re_old=0.d0

    call build(eos_state)
    
    
    
    !compute rhoe_ext
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                if(mask(i,j,k) .eq. 1) then

                    rhoE_old     = uold(i,j,k,UEDEN)
                    rho_e_K_old  = HALF * sum(uold(i,j,k,UMX:UMZ)**2) / uold(i,j,k,URHO)
                    rho          = sum(uold(i,j,k,UFS:UFS+nspec-1))

                    energy       = (rhoE_old - rho_e_K_old) / uold(i,j,k,URHO)
                    re_old(i,j,k) = energy*uold(i,j,k,URHO)

                    ! rho.e source term computed using (rho.E,rho.u,rho)_new rather than pulling from UEINT comp of asrc
                    rho_e_K_new  = HALF * sum(unew(i,j,k,UMX:UMZ)**2)/unew(i,j,k,URHO)
                    rhoedot_ext(i,j,k)  = ( (unew(i,j,k,UEDEN) - rho_e_K_new) &
                                        -   (rho  *  energy) ) / dt_react
                endif
            enddo
        enddo
    enddo

    
    nsubsteps=100
    dt_rk=dt_react/nsubsteps
    npts=(lo(3)-hi(3)+1)*(lo(2)-hi(2)+1)*(lo(1)-hi(1)+1)
    updt_time=0.0

    do steps=1,nsubsteps

        updt_time=updt_time+(steps-1)*dt_rk
        urk_old=urk
        do stage=1,nrkstages

           !Returns the molar production rate of species
           !Given rho, T, and mass fractions
           call VCKWYR(npts, urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO), &
               urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UTEMP), &
               urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UFS:UFS+nspec-1)/urk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO), &
               wdot)

            !Ideally this has to come from a vector call to fuego
            do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                    do i=lo(1),hi(1)

                        if(mask(i,j,k) .eq. 1) then
                            eos_state % rho               = sum(urk(i,j,k,UFS:UFS+nspec-1))
                            rhoInv                        = 1.d0 / eos_state % rho
                            eos_state % massfrac(1:nspec) = urk(i,j,k,UFS:UFS+nspec-1) * rhoInv
                            eos_state % T                 = urk(i,j,k,UTEMP) !guess
                            eos_state % e = (rhoe_init + updt_time * rhoedot_ext(i,j,k)) * rhoInv
                            call eos_re(eos_state)
                            do ns=1,nspec
                                eint(i,j,k,ns)=eos_state%ei(ns)
                            enddo
                            cv(i,j,k)=eos_state%cv
                        endif

                    enddo
                enddo
            enddo

           !convert wdot to gms/s
           wdot(:,:,:,1:nspec) = wdot(:,:,:,1:nspec)*molecular_weight(1:nspec) 

           !update species
           urk(:,:,:,UFS:UFS+nspec-1)  = urk_old(:,:,:,UFS:UFS+nspec-1) + rkcoeffs(stage)*dt_rk*wdot*mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

            !update temperature
            do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                    do i=lo(1),hi(1)
                        if(mask(i,j,k) .eq. 1) then
                            tempsrc=0.0
                            do ns=1,nspec
                                tempsrc=tempsrc+wdot(i,j,k,ns)*eint(i,j,k,ns)
                            enddo
                            tempsrc=-tempsrc/sum(urk(i,j,k,UFS:UFS+nspec-1))/cv(i,j,k)
                            urk(i,j,k,UTEMP) = urk(i,j,k,UTEMP) - rkcoeffs(stage)*tempsrc*mask(i,j,k)
                        endif
                    enddo
                enddo
            enddo

        enddo
    enddo
   
    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                if(mask(i,j,k) .eq. 1) then
                    mom_new = uold(i,j,k,UMX:UMZ) + dt_react*asrc(i,j,k,UMX:UMZ)
                    rho_new = sum(urk(i,j,k,UFS:UFS+nspec-1))
                    rhoe_new = re_old(i,j,k)+dt_react*rhoedot_ext(i,j,k)
                    rho_e_K_new = HALF * sum(mom_new**2) / rho_new
                    rhoE_new    = rhoe_new + rho_e_K_new

    
                    if (do_update .eq. 1) then
                        unew(i,j,k,URHO)            = rho_new
                        unew(i,j,k,UMX:UMZ)         = mom_new
                        unew(i,j,k,UEINT)           = rhoe_new
                        unew(i,j,k,UEDEN)           = rhoE_new
                        unew(i,j,k,UTEMP)           = urk(i,j,k,UTEMP)
                        unew(i,j,k,UFS:UFS+nspec-1) = urk(i,j,k,UFS:UFS+nspec-1)
                    endif
                    ! Add drhoY/dt to reactions MultiFab, but be
                    ! careful because the reactions and state MFs may
                    ! not have the same number of ghost cells.
                    ! Also, be careful because the container for uold might be the same as that for unew
                    if ( i .ge. IR_lo(1) .and. i .le. IR_hi(1) .and. &
                        j .ge. IR_lo(2) .and. j .le. IR_hi(2) .and. &
                        k .ge. IR_lo(3) .and. k .le. IR_hi(3) ) then

                        IR(i,j,k,1:nspec) = (urk(i,j,k,UFS:UFS+nspec-1) - uold(i,j,k,UFS:UFS+nspec-1)) / dt_react -&
                            asrc(i,j,k,UFS:UFS+nspec-1)
                        IR(i,j,k,nspec+1) =(rhoE_new-rhoE_old) / dt_react - asrc(i,j,k,UEDEN)

                    endif

                end if
            enddo
        enddo
    enddo
    

  end subroutine pc_react_state_expl

end module reactions_module
