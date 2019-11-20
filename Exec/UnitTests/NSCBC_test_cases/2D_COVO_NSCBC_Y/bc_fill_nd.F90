module bc_fill_module

  implicit none

! This is a dedicated routine for imposing GC-NSCBC
! (Ghost-Cells Navier-Stokes Characteristic Boundary Conditions)
! For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp. 3399-3408, 2017. 

! PeleC requires ghost cells to be filled but it is not done in AMReX for EXT_DIR.
! So here we fill the corresponding ghost-cells with the target values provided by the user.
! These target values are temporary and ghost-cells will be recomputed with the NSCBC theory
! in the routine 'impose_NSCBC' located in Src_(dim)d


contains

  subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_hypfill")

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
    integer :: bc_type
    double precision :: bc_params(6)

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

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

    !!     ZLO
    !if ( (bc(2,1,1).eq.EXT_DIR).and. adv_lo(2).lt.domlo(2)) then
    !   do i = adv_lo(1), adv_hi(1)
    !      x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
    !      do j = adv_lo(2), domlo(2)-1
    !         x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
    !         do k = adv_lo(3), adv_hi(3)
    !            x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
    !            call bcnormal_nscbc(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,bc_type,bc_params)
    !         end do
    !      end do
    !   end do
    !end if
    !
    !!     ZHI
    !if ( (bc(2,2,1).eq.EXT_DIR).and. adv_hi(2).gt.domhi(2)) then
    !   do i = adv_lo(1), adv_hi(1)
    !      x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
    !      do j = domhi(2)+1, adv_hi(2)
    !         x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
    !         do k = adv_lo(3), adv_hi(3)
    !            x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
    !            call bcnormal_nscbc(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,bc_type,bc_params)
    !         end do
    !      end do
    !   end do
    !end if

  end subroutine pc_hypfill


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

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
  end subroutine pc_reactfill
  
  
    subroutine bcnormal(x,u_int,u_ext,dir,sgn,time,bc_type,bc_params,bc_target)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
    use network, only: nspecies, naux, molec_wt
    use prob_params_module, only : Interior, Inflow, Outflow, SlipWall, NoSlipWall, &
                                   problo, probhi
    
    
    use amrex_constants_module, only: M_PI
    
    implicit none

    double precision :: x(3), time
    double precision :: u_int(*),u_ext(*)
    integer :: dir,sgn
    integer, optional, intent(out) :: bc_type
    double precision, optional, intent(out) :: bc_params(6)
    double precision, optional, intent(out) :: bc_target(5)

    type (eos_t) :: eos_state
    double precision :: u(3)
    double precision :: y
    double precision :: relax_U, relax_V, relax_W, relax_T, beta, sigma_out
    integer :: flag_nscbc, which_bc_type

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
      sigma_out = 0.6d0 ! For outflow only, relax parameter
      which_bc_type = Interior ! This is to ensure that nothing will be done if the user don't set anything
    endif
    
    call build(eos_state)

! at low X
    if (dir == 1) then
      if (sgn == 1) then
      
        relax_U = 1.0d0
        relax_V = 1.0d0
        relax_T = - relax_V
        beta = 0.6d0  
           
        which_bc_type = Inflow
           
        u(1) = u_ref
        u(2) = 0.0d0
        u(3) = 0.0d0
        eos_state % massfrac(1) = 1.d0
        eos_state % p = p_ref
        eos_state % T = T_ref
        call eos_tp(eos_state)
    
      end if
    end if

! at hi X
    if (dir == 1) then
      if (sgn == -1) then
      
        ! Set outflow pressure
        which_bc_type = Outflow      
        sigma_out = 0.28d0
        beta = -1.0d0 ! This means that local mach number will be used
       
        u(1:3) = 0.d0
        eos_state % massfrac(1) = 1.d0
        eos_state % p = p_ref
        eos_state % T = T_ref
        call eos_tp(eos_state)
        
      end if
    end if
    
! at low Y
    if (dir == 2) then
      if (sgn == 1) then
      
        ! Set inflow      
        which_bc_type = Inflow
        relax_U = 1.0d0
        relax_V = 1.0d0
        relax_T = - relax_V
        beta = 0.6d0
           
        u(1) = 0.0d0 
        u(2) = u_ref
        u(3) = 0.0d0
        eos_state % massfrac(1) = 1.d0
        eos_state % p = p_ref
        eos_state % T = T_ref
        call eos_tp(eos_state)
    
      end if
    end if

! at hi Y
    if (dir == 2) then
      if (sgn == -1) then
      
       ! Set outflow pressure
       which_bc_type = Outflow
       sigma_out = 0.28d0
       beta = -1.0d0 ! This means that local mach number will be used
  
       u(1:3) = 0.d0
       eos_state % massfrac(1) = 1.d0
       eos_state % p = p_ref
       eos_state % T = T_ref
       call eos_tp(eos_state)

      end if
    end if


 
       u_ext(UFS:UFS+nspecies-1) = eos_state % massfrac * eos_state % rho
       u_ext(URHO)               = eos_state % rho
       u_ext(UMX)                = eos_state % rho  *  u(1)
       u_ext(UMY)                = eos_state % rho  *  u(2)
       u_ext(UMZ)                = eos_state % rho  *  u(3)
       u_ext(UTEMP)              = eos_state % T
       u_ext(UEINT)              = eos_state % rho  *   eos_state % e
       u_ext(UEDEN)              = eos_state % rho  *  (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2) + u(3)**2)

    
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
      bc_target(3) = 0.0d0
      bc_target(4) = U_ext(UTEMP)
      bc_target(5) = eos_state%p
    end if
 
    call destroy(eos_state)

  end subroutine bcnormal
  

end module bc_fill_module

