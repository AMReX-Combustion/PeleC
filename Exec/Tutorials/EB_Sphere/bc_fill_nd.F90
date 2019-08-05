
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


    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the proper
    ! upwinding.  However, this decision makes this routine look
    ! somewhat non-orthodox, in that we need to set external values in
    ! either case....how do we know it's Outflow?  We have to assume
    ! that the setup routines converted Outflow to FOEXTRAP.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR).and. adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,time)
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR).and. adv_hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, adv_hi(1)
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,time)
             end do
          end do
       end do
    end if

    if (dim .gt. 1) then
       !     YLO
       if ( (bc(2,1,1).eq.EXT_DIR).and. adv_lo(2).lt.domlo(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = adv_lo(2), domlo(2)-1
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,time)
                end do
             end do
          end do
       end if

       !     YHI
       if ( (bc(2,2,1).eq.EXT_DIR).and. adv_hi(2).gt.domhi(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = domhi(2)+1, adv_hi(2)
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,time)
                end do
             end do
          end do
       end if

       if (dim .gt. 2) then
          !     ZLO
          if ( (bc(3,1,1).eq.EXT_DIR).and. adv_lo(3).lt.domlo(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), domlo(3)-1
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,time)
                   end do
                end do
             end do
          end if

          !     ZHI
          if ( (bc(3,2,1).eq.EXT_DIR).and. adv_hi(3).gt.domhi(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = domhi(3)+1, adv_hi(3)
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,time)
                   end do
                end do
             end do
          end if
       end if
    end if


  end subroutine pc_hypfill

  subroutine bcnormal(x,u_int,u_ext,dir,sgn,time,bc_type,bc_params,bc_target)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, NVAR
    use network, only: nspecies, naux
    use prob_params_module, only : problo, probhi, dim, Interior, Inflow, Outflow, SlipWall, NoSlipWall
        
    use amrex_constants_module, only: M_PI
    
    implicit none

    include 'AMReX_bc_types.fi'

    double precision :: x(3)
    double precision :: u_int(NVAR),u_ext(NVAR)
    double precision :: time
    integer :: dir,sgn
    type (eos_t) :: eos_state

    integer, optional, intent(out) :: bc_type
    double precision, optional, intent(out) :: bc_params(6)
    double precision, optional, intent(out) :: bc_target(5)

    double precision :: relax_U, relax_V, relax_W, relax_T, beta, sigma_out
    integer :: flag_nscbc, which_bc_type


    double precision :: u(3), rho_inv
    double precision :: r, rmax, xcen(dim)
    integer :: i

    double precision :: holespace(dim), dist(dim), radius, holeRad
    integer :: nHoles, idholes(dim), d(dim)
    double precision :: pert

    flag_nscbc = 0

    ! When optional arguments are present, GC-NSCBC is activated
    ! Generic values are auto-filled for numerical parameters,
    ! but should be set by the user for each BC
    ! Note that in the impose_NSCBC_xD.f90 routine, not all parameters are used in same time
    if (present(bc_type).and.present(bc_params).and.present(bc_target)) then

      flag_nscbc = 1
      relax_U = 0.5d0 ! For inflow only, relax parameter for x_velocity
      relax_V = 0.5d0 ! For inflow only, relax parameter for y_velocity
      relax_W = 0.5d0 ! For inflow only, relax parameter for z_velocity
      relax_T = 0.2d0 ! For inflow only, relax parameter for temperature
      beta = 0.2d0  ! Control the contribution of transverse terms
      sigma_out = 0.25d0 ! For outflow only, relax parameter
      which_bc_type = Interior ! This is to ensure that nothing will be done if the user don't set anything
    endif

    call build(eos_state)

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

       eos_state % massfrac(1:nspecies) = 0.0
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

          u_ext(UMX:UMZ) = eos_state % rho  *  u(1:3)
          u_ext(URHO) = eos_state % rho
          u_ext(UEINT) = eos_state % rho  *  eos_state % e
          u_ext(UEDEN) = eos_state % rho  * (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2 + u(3)**2))
          u_ext(UTEMP) = eos_state % T
          u_ext(UFS:UFS+nspecies-1) = eos_state % rho  *  eos_state % massfrac(1:nspecies)

       

      ! Here the optional parameters are filled by the local variables if they were present
    if (flag_nscbc == 1) then
      bc_type = which_bc_type
      bc_params(1) = relax_T
      bc_params(2) = relax_U
      bc_params(3) = relax_V
      bc_params(4) = relax_W
      bc_params(5) = beta
      bc_params(6) = sigma_out
      bc_target(1) = U_ext(UMX)/U_ext(URHO)
      bc_target(2) = U_ext(UMY)/U_ext(URHO)
      bc_target(3) = U_ext(UMZ)/U_ext(URHO)
      bc_target(4) = U_ext(UTEMP)
      bc_target(5) = eos_state%p
    end if



  call destroy(eos_state)

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


