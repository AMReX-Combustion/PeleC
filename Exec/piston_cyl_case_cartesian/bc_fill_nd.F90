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
    logical rho_only

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the proper
    ! upwinding.  However, this decision makes this routine look
    ! somewhat non-orthodox, in that we need to set external values in
    ! either case....how do we know it's Outflow?  We have to assume
    ! that the setup routines converted Outflow to FOEXTRAP.

    !     Set flag for bc function
    rho_only = .FALSE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR).and. adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,rho_only)
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
                call bcnormal(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,rho_only)
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
                   call bcnormal(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,rho_only)
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
                   call bcnormal(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,rho_only)
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
                      call bcnormal(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,rho_only)
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
                      call bcnormal(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,rho_only)
                   end do
                end do
             end do
          end if
       end if
    end if

  end subroutine pc_hypfill


    subroutine pc_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_denfill")

    use prob_params_module, only: dim  
    use eos_type_module
    use eos_module

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

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    rho_only = .TRUE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR).and.adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k),adv(i,j,k),1,+1,rho_only)
             end do
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR).and.adv_hi(1).gt.domhi(1)) then
       do i = domhi(1)+1, adv_hi(1)
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domhi(1),j,k),adv(i,j,k),1,-1,rho_only)
             end do
          end do
       end do
    end if

    if (dim > 1) then
       !     YLO
       if ( (bc(2,1,1).eq.EXT_DIR).and.adv_lo(2).lt.domlo(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = adv_lo(2), domlo(2)-1
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domlo(2),k),adv(i,j,k),2,+1,rho_only)
                end do
             end do
          end do
       end if

       !     YHI
       if ( (bc(2,2,1).eq.EXT_DIR).and.adv_hi(2).gt.domhi(2)) then
          do i = adv_lo(1), adv_hi(1)
             x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
             do j = domhi(2)+1, adv_hi(2)
                x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                do k = adv_lo(3), adv_hi(3)
                   x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                   call bcnormal(x,adv(i,domhi(2),k),adv(i,j,k),2,-1,rho_only)
                end do
             end do
          end do
       end if

       if (dim > 2) then
          !     ZLO
          if ( (bc(3,1,1).eq.EXT_DIR).and.adv_lo(3).lt.domlo(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = adv_lo(2), domlo(2)-1
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), adv_hi(3)
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domlo(3)),adv(i,j,k),3,+1,rho_only)
                   end do
                end do
             end do
          end if

          !     ZHI
          if ( (bc(3,2,1).eq.EXT_DIR).and.adv_hi(3).gt.domhi(3)) then
             do i = adv_lo(1), adv_hi(1)
                x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
                do j = domhi(2)+1, adv_hi(2)
                   x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
                   do k = adv_lo(3), adv_hi(3)
                      x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                      call bcnormal(x,adv(i,j,domhi(3)),adv(i,j,k),3,-1,rho_only)
                   end do
                end do
             end do
          end if
       end if
    end if
  end subroutine pc_denfill


  subroutine bcnormal(x,u_int,u_ext,dir,sgn,rho_only)

    use probdata_module
    use network, only: nspec
    use bl_constants_module, only: M_PI
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use eos_type_module
    use eos_module

    implicit none

    double precision :: x(3)
    double precision :: u_int(*),u_ext(*)
    logical rho_only

    integer :: dir,sgn
    type (eos_t) :: eos_state_ext,eos_state_int

    double precision :: phi_angle,hole_cx,hole_cz,dphi_angle
    double precision :: vx_in, vz_in, vr_in, vy_in

    integer :: nh

    !default wall condition
    if(rho_only .eqv. .true.) then
        u_ext(URHO)=u_int(URHO)
    else
        u_ext(URHO)            =  u_int(URHO)
        u_ext(UMX)             = -u_int(UMX)
        u_ext(UMY)             = -u_int(UMY)
        u_ext(UMZ)             = -u_int(UMZ)
        u_ext(UEINT)           = u_int(UEINT)
        u_ext(UEDEN)           = u_int(UEDEN)
        u_ext(UFS:UFS+nspec-1) = u_int(UFS:UFS+nspec-1)
        u_ext(UTEMP)           = u_int(UTEMP)
    endif


    call build(eos_state_ext)
    call build(eos_state_int)


    !I am assuming this boundary is the top XZ plane
    ! and in-plane angle is measured from X axis
    dphi_angle=2.d0*M_PI/nholes

    do nh=1,nholes

        phi_angle = nh*dphi_angle
        hole_cx = centx + r_circ*cos(phi_angle)
        hole_cz = centz + r_circ*sin(phi_angle)

        if( ((x(1)-hole_cx)**2+(x(3)-hole_cz)**2) < r_hole**2) then

                eos_state_int % rho      = u_int(URHO)
                eos_state_int % T        = u_int(UTEMP)
                eos_state_int % massfrac = u_int(UFS:UFS+nspec-1)/u_int(URHO)
                call eos_rt(eos_state_int)
                
                !subsonic velocity inlet
                eos_state_ext % rho = dens_jet
                eos_state_ext % p   = eos_state_int % p !use interior pressure
                eos_state_ext % massfrac = 0.d0
                eos_state_ext % massfrac(fuel_id) = Yfuel_jet
                eos_state_ext % massfrac(ox_id)   = Yox_jet
                eos_state_ext % massfrac(N2_id)   = YN2_jet

                call eos_rp(eos_state_ext)

                !find velocity vector
                vy_in = -vel_jet*cos(cone_angle*M_PI/180.d0) !top XZ plane, remember!
                vr_in =  vel_jet*sin(cone_angle*M_PI/180.d0)
                vx_in = vr_in*cos(phi_angle)
                vz_in = vr_in*sin(phi_angle)

                if(rho_only .eqv. .true.) then
                        u_ext(URHO)   = eos_state_ext % rho
                else
                        u_ext(URHO)   = eos_state_ext % rho
                        u_ext(UMX)    = eos_state_ext % rho * vx_in
                        u_ext(UMY)    = eos_state_ext % rho * vy_in
                        u_ext(UMZ)    = eos_state_ext % rho * vz_in
                        u_ext(UEINT)  = eos_state_ext % rho * eos_state_ext % e
                        u_ext(UEDEN)  = u_ext(UEINT) + 0.5d0*eos_state_ext % rho*vel_jet**2
                        u_ext(UFS:UFS+nspec-1)  = eos_state_ext % rho * eos_state_ext % massfrac(:)
                        u_ext(UTEMP)            = eos_state_ext % T
                endif

        endif

    enddo
    
    call destroy(eos_state_int)
    call destroy(eos_state_ext)

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
