
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
    use meth_params_module, only : URHO, UMX, UMZ, UTEMP, UEDEN, UEINT, UFS, NVAR
    use network, only: nspec, naux
    use prob_params_module, only : problo, probhi, dim
        
    use bl_constants_module, only: M_PI
    
    implicit none

    include 'AMReX_bc_types.fi'

    double precision :: x(3)
    double precision :: u_int(NVAR),u_ext(NVAR)
    double precision :: time
    logical rho_only
    integer :: dir,sgn,bc
    type (eos_t) :: eos_state

    double precision :: u(3), rho_inv
    double precision :: r, rmax, xcen(dim)
    integer :: i

    double precision :: holespace(dim), dist(dim), radius, holeRad
    integer :: nHoles, idholes(dim), d(dim)
    double precision :: pert

    nHoles = 40
    holespace = (probhi - problo) / DBLE(nHoles)
    idholes = INT(x(1:dim)/holespace)

    dist = x - (idholes + 0.5d0)*holespace
    holeRad = .04d0

    d = 0
    d(3) = 1
    radius = SQRT(d(1)*(dist(2)**2 + dist(3)**2) &
         +        d(2)*(dist(1)**2 + dist(3)**2) &
         +        d(3)*(dist(1)**2 + dist(2)**2))


    if (bc.eq.EXT_DIR) then ! inflow: extrapolate none, set (Y,u,p,T), compute (rho,e)
       r = 0.0
       rmax = 0.0
       do i = 1, dim
          xcen(i) = 0.5*(problo(i)+probhi(i))
          if (i .ne. dir) then
             r = r + (x(i)-xcen(i))*(x(i)-xcen(i))
             rmax = rmax + (problo(i)-xcen(i))*(problo(i)-xcen(i))
          endif
       end do
       r = sqrt(r)
       rmax = sqrt(rmax)
       
       pert = sin(7*M_pi*x(1)/rmax + .5d0)*cos(3*M_pi*x(2)/rmax+7.d0)*sin(time/1.d-6) &
            + sin(3*M_pi*x(1)/rmax + .5d0)*cos(4*M_pi*x(2)/rmax-3d0)*sin(2.d0*time/1.d-6-.1d0) &
            + sin(2*M_pi*x(1)/rmax + 3*M_pi*x(2)/rmax+4.8d0)*cos(.33d0*time/1.d-6) &
            + sin(5*M_pi*x(1)/rmax - 3*M_pi*x(2)/rmax-.17d0)*cos(time/1.d-6+2.d0)

       eos_state % massfrac(1:nspec) = 0.0
       eos_state % massfrac(1) = 1.0
       eos_state % p = pamb

       if (r/Rmax .gt. Rfrac ) then
          eos_state % T = Tin
          u = 0
       else
          eos_state % T = Tin
          u(1) =   ( uin * cos(swrlang*M_PI/180.0) ) * ( x(2)-xcen(2) )/r
          u(2) = - ( uin * cos(swrlang*M_PI/180.0) ) * ( x(1)-xcen(1) )/r
          u(3) = uin * sin(swrlang*M_PI/180.0)
       endif
       u(3) = u(3)*(1.d0+pertmag*pert)

       call eos_tp(eos_state)

       if (rho_only .eqv. .true.) then
          u_ext(1) = eos_state % rho
       else
          u_ext(UMX:UMZ) = eos_state % rho  *  u(1:3)
          u_ext(URHO) = eos_state % rho
          u_ext(UEINT) = eos_state % rho  *  eos_state % e
          u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
          u_ext(UTEMP) = eos_state % T
          u_ext(UFS:UFS+nspec-1) = eos_state % rho  *  eos_state % massfrac(1:nspec)
       endif

       
    else ! ouflow: extrapolate (rho,Y,u), set p, compute (rho,e,T)
       if (rho_only .eqv. .true.) then
          ! Don't have species, assume it's consistent
          u_ext(1) = u_int(1)
       else
          eos_state % p = pamb
          eos_state % rho = u_int(URHO)
          rho_inv = 1.d0 / eos_state % rho
          eos_state % massfrac = u_int(UFS:UFS+nspec-1) * rho_inv
          u_ext(UMX:UMZ) = u_int(UMX:UMZ)
          u(1:3) = u_ext(UMX:UMZ) * rho_inv
          call eos_rp(eos_state)

          u_ext(URHO) = eos_state % rho
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

    integer :: n

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
  end subroutine pc_reactfill

end module bc_fill_module

