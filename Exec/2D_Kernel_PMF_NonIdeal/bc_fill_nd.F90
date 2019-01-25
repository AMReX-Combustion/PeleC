module bc_fill_module

  implicit none

  public

contains

  subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_hypfill")

    use meth_params_module, only: NVAR,URHO
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

    type(eos_t) :: eos_state

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
                call bcnormal(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,rho_only,eos_state)
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
                call bcnormal(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,rho_only,eos_state)
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
                   call bcnormal(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,rho_only,eos_state)
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
                   call bcnormal(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,rho_only,eos_state)
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
                      call bcnormal(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,rho_only,eos_state)
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
                      call bcnormal(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,rho_only,eos_state)
                   end do
                end do
             end do
          end if
       end if
    end if

    call destroy(eos_state)

  end subroutine pc_hypfill


  subroutine bcnormal(x,u_int,u_ext,dir,sgn,rho_only,eos_state)

    use probdata_module
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS, NVAR
    use prob_params_module, only: dim
    use network, only : nspec
    
    implicit none

    double precision :: x(3)
    double precision :: u_int(NVAR),u_ext(NVAR)
    logical :: rho_only
    integer :: dir,sgn
    type(eos_t) :: eos_state
    double precision :: rho_inv, u(3)

    if (.not. pmf_initialized) then
       call init_pmf()
    end if

    if (rho_only .EQV. .TRUE. ) then
       if (dir.eq.dim .and. sgn.eq.+1) then
          u_ext(1) = fuel_state(URHO)
       endif
    else
       if (sgn.eq.-1) then
          !  Use extrapolated rho, rhoY/rho and ambient p (plus guess for T)
          rho_inv = 1.d0 / u_ext(URHO)
          eos_state % rho = u_ext(URHO)
          eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) * rho_inv
          eos_state % T = u_ext(UTEMP)
          eos_state % p = pamb
          u = u_ext(UMX:UMZ) * rho_inv

          call eos_rp(eos_state)

          u_ext(URHO) = eos_state % rho
          u_ext(UMX:UMZ) = eos_state % rho  *  u
          u_ext(UEINT) = eos_state % rho  *  eos_state % e
          u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
          u_ext(UTEMP) = eos_state % T
          u_ext(UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)          
       else
          u_ext(1:NVAR) = fuel_state(1:NVAR)

          rho_inv = 1.d0 / u_ext(URHO)
          eos_state % rho = u_ext(URHO)
          eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) * rho_inv
          eos_state % T = u_ext(UTEMP)
          eos_state % p = pamb

          call eos_rp(eos_state)

          u = 0.d0
          u(dim) = 50.d0
          u_ext(UMX:UMZ) = eos_state % rho  *  u
          u_ext(UEINT) = eos_state % rho  *  eos_state % e
          u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
          u_ext(UTEMP) = eos_state % T
          u_ext(UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec) 
       endif
    endif

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

    double precision :: x(3)
    integer :: i, j, k, n
    logical rho_only

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
  end subroutine pc_reactfill

end module bc_fill_module
