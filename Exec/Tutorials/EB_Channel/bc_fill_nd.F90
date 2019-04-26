module bc_fill_module

  implicit none

  public

contains

  subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_hypfill")

    use meth_params_module, only: NVAR
    use prob_params_module, only: dim
    use eos_type_module
    use eos_module

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
                do j = adv_lo(2), domlo(2)-1
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), adv_hi(3)
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
                do j = domhi(2)+1, adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), adv_hi(3)
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
    use network, only: nspec
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use amrex_constants_module, only: HALF, M_PI
    use prob_params_module, only : Interior, Inflow, Outflow, SlipWall, NoSlipWall, &
                                   problo, probhi, dim
    use eos_type_module
    use eos_module

    implicit none

    double precision :: x(3),time
    double precision :: u_int(*),u_ext(*)
   
    integer, optional, intent(out) :: bc_type
    double precision, optional, intent(out) :: bc_params(6)
    double precision, optional, intent(out) :: bc_target(5)

    integer :: dir,sgn
    type (eos_t) :: eos_state, eos_state_int
    double precision :: u, v, w
    double precision :: cost, sint, umag_int
    double precision :: relax_U, relax_V, relax_W, relax_T, beta, sigma_out
    integer :: flag_nscbc, which_bc_type

    ! flag_nscbc = 0

    ! if (present(bc_type).and.present(bc_params).and.present(bc_target)) then

    !   flag_nscbc = 1
    !   relax_U = 0.5d0 ! For inflow only, relax parameter for x_velocity
    !   relax_V = 0.5d0 ! For inflow only, relax parameter for y_velocity
    !   relax_W = 0.5d0 ! For inflow only, relax parameter for z_velocity
    !   relax_T = 0.2d0 ! For inflow only, relax parameter for temperature
    !   beta = 0.2d0  ! Control the contribution of transverse terms
    !   sigma_out = 0.25d0 ! For outflow only, relax parameter
    !   which_bc_type = Interior ! This is to ensure that nothing will be done if the user don't set anything
    ! endif

    ! call build(eos_state)

    ! if (dir == 1) then
    !    if (sgn == 1) then
    !       relax_U = 10.0d0
    !       relax_V = 10.0d0
    !       relax_T = - relax_V
    !       beta = 0.5d0
          
    !       which_bc_type = Inflow

    !       u = u0
    !       v = v0
    !       w = w0
    !       eos_state % rho = u_int(URHO)
    !       eos_state % T = u_int(UTEMP) 
    !       eos_state % massfrac = 0.d0
    !       eos_state % massfrac(nspec) = 1.d0
    !       call eos_rt(eos_state)
    !       ! eos_state % massfrac(1) = 1.d0
    !       ! eos_state % p = p0
    !       ! eos_state % T = T0
    !       ! call eos_tp(eos_state)
          
    !    else if (sgn == -1) then

    !       which_bc_type = Outflow
    !       sigma_out = 0.28d0
    !       beta = -0.60d0

    !       u = 0.d0
    !       v = 0.d0
    !       w = 0.d0          
    !       eos_state % massfrac(1) = 1.d0
    !       eos_state % rho = rho0
    !       eos_state % T = T0
    !       call eos_rt(eos_state)
          
    !    end if
    ! end if

    ! u_ext(URHO)            = eos_state % rho
    ! u_ext(UFS:UFS+nspec-1) = eos_state % massfrac * eos_state % rho
    ! u_ext(UMX)             = eos_state % rho * u
    ! u_ext(UMY)             = eos_state % rho * v
    ! u_ext(UMZ)             = eos_state % rho * w
    ! u_ext(UTEMP)           = eos_state % T
    ! u_ext(UEINT)           = eos_state % rho *  eos_state % e
    ! u_ext(UEDEN)           = eos_state % rho * (eos_state % e + 0.5d0 * (u**2 + v**2 + w**2))

    ! ! Here the optional parameters are filled by the local variables if they were present
    ! if (flag_nscbc == 1) then
    !    bc_type = which_bc_type
    !    bc_params(1) = relax_T !  For inflow only, relax parameter for temperature
    !    bc_params(2) = relax_U ! For inflow only, relax parameter for x_velocity
    !    bc_params(3) = relax_V ! For inflow only, relax parameter for y_velocity
    !    bc_params(4) = relax_W ! For inflow only, relax parameter for z_velocity
    !    bc_params(5) = beta  ! Control the contribution of transverse terms.
    !    bc_params(6) = sigma_out ! For outflow only, relax parameter
    !    bc_target(1) = u_ext(UMX)/u_ext(URHO)  ! Target for Inflow
    !    bc_target(2) = u_ext(UMY)/u_ext(URHO)  ! Target for Inflow
    !    bc_target(3) = u_ext(UMZ)/u_ext(URHO)  ! Target for Inflow
    !    bc_target(4) = u_ext(UTEMP)            ! Target for Inflow
    !    bc_target(5) = eos_state%p             ! Target for Outflow
    ! end if

    ! ==================================================================
    call build(eos_state)

    cost = cos(M_PI/180.d0 * angle)
    sint = sin(M_PI/180.d0 * angle)
    
    if (dir == 1) then
       if (sgn == 1) then
          eos_state % p = dpdx * (x(1) * cost + x(2) * sint) + p0
          eos_state % T = T0
          eos_state % massfrac    = 0.d0
          eos_state % massfrac(1) = 1.d0
          call eos_tp(eos_state)

          u = u_int(UMX) / u_int(URHO)
          v = u_int(UMY) / u_int(URHO)
          w = u_int(UMZ) / u_int(URHO)
          umag_int = sqrt(u**2 + v**2 + w**2)
          u = umag_int * cost
          v = umag_int * sint
          w = 0.d0

       else if (sgn == -1) then
          eos_state % p = dpdx * (x(1) * cost + x(2) * sint) + p0
          eos_state % T = T0
          eos_state % massfrac    = 0.d0
          eos_state % massfrac(1) = 1.d0
          call eos_tp(eos_state)

          u = u_int(UMX) / u_int(URHO)
          v = u_int(UMY) / u_int(URHO)
          w = u_int(UMZ) / u_int(URHO)
       end if
    end if

    u_ext(URHO)            = eos_state % rho
    u_ext(UFS:UFS+nspec-1) = eos_state % massfrac * eos_state % rho
    u_ext(UMX)             = eos_state % rho * u
    u_ext(UMY)             = eos_state % rho * v
    u_ext(UMZ)             = eos_state % rho * w
    u_ext(UTEMP)           = eos_state % T
    u_ext(UEINT)           = eos_state % rho *  eos_state % e
    u_ext(UEDEN)           = eos_state % rho * (eos_state % e + 0.5d0 * (u**2 + v**2 + w**2))

    ! ==================================================================
    ! eos_state % p = p0
    ! eos_state % T = T0
    ! eos_state % massfrac    = 0.d0
    ! eos_state % massfrac(1) = 1.d0
    ! call eos_tp(eos_state)

    ! u_ext(URHO)            = rho0
    ! u_ext(UFS:UFS+nspec-1) = rho0 * eos_state % massfrac(1:nspec)
    ! u_ext(UMX)             = rho0 * u0
    ! u_ext(UMY)             = rho0 * v0
    ! u_ext(UMZ)             = rho0 * w0
    ! u_ext(UEINT)           = rho0 * eint0
    ! u_ext(UEDEN)           = rho0 * (eint0 + HALF * (u0**2 + v0**2 + w0**2))
    ! u_ext(UTEMP)           = T0
    
    call destroy(eos_state)

  end subroutine bcnormal

end module bc_fill_module
