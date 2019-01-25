module bc_fill_module

  implicit none

  public

contains


  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.
  
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

    integer          :: n

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine pc_hypfill

  subroutine pc_reactfill(react,react_lo,react_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_reactfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: react_lo(3),react_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

    call filcc_nd(react,react_lo,react_hi,domlo,domhi,delta,xlo,bc)

  end subroutine pc_reactfill


!---------------
! This is a dedicated routine for imposing GC-NSCBC
! (Ghost-Cells Navier-Stokes Characteristic Boundary Conditions)
! For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp. 3399-3408, 2017. 

! This routine provides the boundary type (bc_type),
! the physical target values (bc_target)
! and numerical parameters associated to the BC (bc_params)

! Gghost-cells will be recomputed with the NSCBC theory
! in the routine 'impose_NSCBC' located in Src_(dim)d

  subroutine bcnormal(x,u_int,u_ext,dir,sgn,bc_type,bc_params,bc_target)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
    use network, only: nspec, naux, molec_wt
    use prob_params_module, only : Interior, Inflow, Outflow, SlipWall, NoSlipWall, &
                                   problo, probhi
    
    
    use bl_constants_module, only: M_PI
    
    implicit none
  
    double precision :: x(3)
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
  
  
    ! Below is an example where we impose outflow condition everywhere
    ! The user can impose any combination of BCs on all faces with the help of idir and isign
    which_bc_type = Outflow
    sigma_out = 0.28d0
    beta = -1.0d0 ! This means that local mach number will be used
  
    u(1) = u_int(UMX)
    u(2) = u_int(UMY)
    u(3) = u_int(UMZ) 
    eos_state % massfrac = u_int(UFS:UFS+nspec-1)
    eos_state % p = 101325.d0    ! This is a target value
    eos_state % T = u_int(UTEMP)
    call eos_tp(eos_state)
    

    u_ext(UFS:UFS+nspec-1) = eos_state % massfrac * eos_state % rho
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
      bc_target(3) = U_ext(UMZ)/U_ext(URHO)
      bc_target(4) = U_ext(UTEMP)
      bc_target(5) = eos_state%p
    end if
  
    call destroy(eos_state)

  end subroutine bcnormal

end module bc_fill_module
