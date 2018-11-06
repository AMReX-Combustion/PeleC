module bc_fill_module

  implicit none

  public

contains

  subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_hypfill")

    use meth_params_module, only: NVAR
    use prob_params_module, only: dim
    use eos_type_module

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    double precision :: x(3)
    integer :: i, j, k, n
    logical rho_only

    type (eos_t) :: eos_state

    call build(eos_state)

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the proper
    ! upwinding.  However, this decision makes this routine look
    ! somewhat non-orthodox, in that we need to set external values in
    ! either case....how do we know it's Outflow?  We have to assume
    ! that the setup routines converted Outflow to FOEXTRAP.

    ! Set flag for bc function
    rho_only = .FALSE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and. adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,bc(1,1,1),eos_state,rho_only,time)
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and. adv_hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, adv_hi(1)
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,bc(1,2,1),eos_state,rho_only,time)
             end do
          end do
       end do
    end if

    if (dim .gt. 1) then
       !     YLO
       if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and. adv_lo(2).lt.domlo(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = adv_lo(2), domlo(2)-1
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,bc(2,1,1),eos_state,rho_only,time)
                end do
             end do
          end do
       end if

       !     YHI
       if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and. adv_hi(2).gt.domhi(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = domhi(2)+1, adv_hi(2)
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,bc(2,2,1),eos_state,rho_only,time)
                end do
             end do
          end do
       end if

       if (dim .gt. 2) then
          !     ZLO
          if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and. adv_lo(3).lt.domlo(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), domlo(3)-1
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,bc(3,1,1),eos_state,rho_only,time)
                   end do
                end do
             end do
          end if

          !     ZHI
          if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and. adv_hi(3).gt.domhi(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = domhi(3)+1, adv_hi(3)
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,bc(3,2,1),eos_state,rho_only,time)
                   end do
                end do
             end do
          end if
       end if
    end if

    call destroy(eos_state)

  end subroutine pc_hypfill

  subroutine pc_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_denfill")

    use prob_params_module, only: dim
    use eos_type_module

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    double precision :: x(3)
    logical rho_only
    integer :: i,j,k

    type (eos_t) :: eos_state

    call build(eos_state)

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    rho_only = .TRUE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k),adv(i,j,k),1,+1,bc(1,1,1),eos_state,rho_only,time)
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, adv_hi(1)
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domhi(1),j,k),adv(i,j,k),1,-1,bc(1,2,1),eos_state,rho_only,time)
             end do
          end do
       end do
    end if

    if (dim > 1) then
       !     YLO
       if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_lo(2).lt.domlo(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = adv_lo(2), domlo(2)-1
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domlo(2),k),adv(i,j,k),2,+1,bc(2,1,1),eos_state,rho_only,time)
                end do
             end do
          end do
       end if

       !     YHI
       if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_hi(2).gt.domhi(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = domhi(2)+1, adv_hi(2)
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domhi(2),k),adv(i,j,k),2,-1,bc(2,2,1),eos_state,rho_only,time)
                end do
             end do
          end do
       end if

       if (dim > 2) then
          !     ZLO
          if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and.adv_lo(3).lt.domlo(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), domlo(3)-1
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domlo(3)),adv(i,j,k),3,+1,bc(3,1,1),eos_state,rho_only,time)
                   end do
                end do
             end do
          end if

          !     ZHI
          if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and.adv_hi(3).gt.domhi(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = domhi(3)+1, adv_hi(3)
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domhi(3)),adv(i,j,k),3,-1,bc(3,2,1),eos_state,rho_only,time)
                   end do
                end do
             end do
          end if
       end if
    end if

    call destroy(eos_state)

  end subroutine pc_denfill

  subroutine bcnormal(x,u_int,u_ext,dir,sgn,bc,eos_state,rho_only,time)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, NVAR
    use network, only: nspec, naux
    use prob_params_module, only : problo, probhi, dim

    use bl_constants_module, only: M_PI, HALF

#ifdef USE_MASA
    use masa
#endif

    implicit none

    include 'AMReX_bc_types.fi'

    double precision :: x(3),vint(3)
    double precision :: u_int(NVAR),u_ext(NVAR)
    double precision :: time
    logical rho_only
    integer :: dir,sgn,bc
    type (eos_t) :: eos_state
    double precision :: rho,u,v,w,p,eint

#ifdef USE_MASA

    call build(eos_state)


    ! inflow
    if(bc .eq. EXT_DIR .and. sgn .eq. 1) then

       rho = masa_eval_3d_exact_rho(x(1),x(2),x(3))

       if(rho_only .eqv. .true.) then
          u_ext(URHO) = rho
       else
          u = masa_eval_3d_exact_u(x(1),x(2),x(3))
          v = masa_eval_3d_exact_v(x(1),x(2),x(3))
          w = masa_eval_3d_exact_w(x(1),x(2),x(3))
          p = masa_eval_3d_exact_p(x(1),x(2),x(3))

          eos_state % rho = rho
          eos_state % p = p
          eos_state % massfrac    = 0.d0
          eos_state % massfrac(1) = 1.d0
          call eos_rp(eos_state)
          eint = eos_state % e

          u_ext(URHO)            = rho
          u_ext(UFS:UFS+nspec-1) = rho * eos_state % massfrac(1:nspec)
          u_ext(UMX)             = rho * u
          u_ext(UMY)             = rho * v
          u_ext(UMZ)             = rho * w
          u_ext(UTEMP)           = eos_state % T
          u_ext(UEINT)           = rho * eint
          u_ext(UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))
       endif

    ! outflow
    else if(bc .eq. EXT_DIR .and. sgn .eq. -1) then

       rho = masa_eval_3d_exact_rho(x(1),x(2),x(3))

       if(rho_only .eqv. .true.) then
          u_ext(URHO) = rho
       else
          u = masa_eval_3d_exact_u(x(1),x(2),x(3))
          v = masa_eval_3d_exact_v(x(1),x(2),x(3))
          w = masa_eval_3d_exact_w(x(1),x(2),x(3))
          p = masa_eval_3d_exact_p(x(1),x(2),x(3))

          eos_state % rho = rho
          eos_state % p = p
          eos_state % massfrac    = 0.d0
          eos_state % massfrac(1) = 1.d0
          call eos_rp(eos_state)
          eint = eos_state % e

          u_ext(URHO)            = rho
          u_ext(UFS:UFS+nspec-1) = rho * eos_state % massfrac(1:nspec)
          u_ext(UMX)             = rho * u
          u_ext(UMY)             = rho * v
          u_ext(UMZ)             = rho * w
          u_ext(UTEMP)           = eos_state % T
          u_ext(UEINT)           = rho * eint
          u_ext(UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))
       endif

    endif

#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif
  end subroutine bcnormal

  subroutine pc_reactfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_reactfill")

    use meth_params_module, only: NVAR
    use prob_params_module, only: dim

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer :: n

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
  end subroutine pc_reactfill

end module bc_fill_module
