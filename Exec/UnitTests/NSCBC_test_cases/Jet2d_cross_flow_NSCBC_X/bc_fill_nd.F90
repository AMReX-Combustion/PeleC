module bc_fill_module

  implicit none

  public :: bcnormal

! This is a dedicated module for GC-NSCBC
! (Ghost-Cells Navier-Stokes Characteristic Boundary Conditions)
! See Motheau et al. AIAA J. (Submitted) for the theory. 

! PeleC requires ghost cells to be filled but it is not done in AMReX for EXT_DIR.
! So here we fill the corresponding ghost-cells with the target values provided by the user.
! These target values are temporary and ghost-cells will be recomputed with the NSCBC theory
! in the routines 'impose_NSCBC_something' located in Src_(dim)d


contains

  subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_hypfill")

    use meth_params_module, only: NVAR, nb_nscbc_params
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
    integer :: bc_type
    double precision :: bc_params(nb_nscbc_params)

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! Set flag for bc function
    rho_only = .FALSE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR).and. adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,bc_type,bc_params,rho_only)
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
                call bcnormal(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,bc_type,bc_params,rho_only)
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
                call bcnormal(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,bc_type,bc_params,rho_only)
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
                call bcnormal(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,bc_type,bc_params,rho_only)
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
    !            call bcnormal(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,bc_type,bc_params,rho_only)
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
    !            call bcnormal(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,bc_type,bc_params,rho_only)
    !         end do
    !      end do
    !   end do
    !end if

  end subroutine pc_hypfill



  subroutine pc_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="pc_denfill")

    use prob_params_module, only: dim
    use meth_params_module, only: nb_nscbc_params

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
    integer :: bc_type
    double precision :: bc_params(nb_nscbc_params)

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

    rho_only = .TRUE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR).and.adv_lo(1).lt.domlo(1)) then
       do i = adv_lo(1), domlo(1)-1
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), adv_hi(2)
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(domlo(1),j,k),adv(i,j,k),1,+1,bc_type,bc_params,rho_only)
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
                call bcnormal(x,adv(domhi(1),j,k),adv(i,j,k),1,-1,bc_type,bc_params,rho_only)
             end do
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1,1).eq.EXT_DIR).and.adv_lo(2).lt.domlo(2)) then
       do i = adv_lo(1), adv_hi(1)
          x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
          do j = adv_lo(2), domlo(2)-1
             x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
             do k = adv_lo(3), adv_hi(3)
                x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
                call bcnormal(x,adv(i,domlo(2),k),adv(i,j,k),2,+1,bc_type,bc_params,rho_only)
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
                call bcnormal(x,adv(i,domhi(2),k),adv(i,j,k),2,-1,bc_type,bc_params,rho_only)
             end do
          end do
       end do
    end if

    !!     ZLO
    !if ( (bc(3,1,1).eq.EXT_DIR).and.adv_lo(3).lt.domlo(3)) then
    !   do i = adv_lo(1), adv_hi(1)
    !      x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
    !      do j = adv_lo(2), domlo(2)-1
    !         x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
    !         do k = adv_lo(3), adv_hi(3)
    !            x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
    !            call bcnormal(x,adv(i,j,domlo(3)),adv(i,j,k),3,+1,bc_type,bc_params,rho_only)
    !         end do
    !      end do
    !   end do
    !end if
    !
    !!     ZHI
    !if ( (bc(3,2,1).eq.EXT_DIR).and.adv_hi(3).gt.domhi(3)) then
    !   do i = adv_lo(1), adv_hi(1)
    !      x(1) = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5d0)
    !      do j = domhi(2)+1, adv_hi(2)
    !         x(2) = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5d0)
    !         do k = adv_lo(3), adv_hi(3)
    !            x(3) = xlo(3) + delta(3)*(dble(k-adv_lo(3)) + 0.5d0)
    !            call bcnormal(x,adv(i,j,domhi(3)),adv(i,j,k),3,-1,bc_type,bc_params,rho_only)
    !         end do
    !      end do
    !   end do
    !end if
    
  end subroutine pc_denfill



  subroutine bcnormal(x,u_int,u_ext,dir,sgn,bc_type,bc_params,rho_only)

    use probdata_module
    use eos_type_module
    use eos_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
    use meth_params_module, only: nb_nscbc_params
    use network, only: nspec, naux, molec_wt
    use prob_params_module, only : Interior, Inflow, Outflow, SlipWall, NoSlipWall, &
                                   problo, probhi
    
    
    use bl_constants_module, only: M_PI
    
    implicit none

    double precision :: x(3)
    double precision :: u_int(*),u_ext(*)
    logical rho_only
    integer :: dir,sgn
    integer, intent(out) :: bc_type
    double precision, intent(out) :: bc_params(nb_nscbc_params)

    type (eos_t) :: eos_state
    double precision :: u(3)
    double precision :: y
    double precision :: relax_U, relax_V, relax_T, beta, sigma_out

    call build(eos_state)

! at low Y
    if (dir == 2) then
      if (sgn == 1) then
       
         if ((x(1) <= 0.0175d0) .or. (x(1) >= 0.0225d0)) then
           bc_type = SlipWall 
           sigma_out = 0.28d0
           beta = .0d0
           bc_params(1) = sigma_out
           bc_params(4) = beta
           !
           u(1) = 0.0d0
           u(2) = 0.0d0
           u(3) = 0.0d0
           eos_state % massfrac(1) = 1.d0
           eos_state % p = p_ref
           eos_state % T = T_ref
           call eos_tp(eos_state)
         else
           
           relax_U = .0d0
           relax_V = 2.d0
           relax_T = - .0d0 !relax_V
           beta = 1.0d0
           
           bc_type = Inflow
           bc_params(1) = relax_U
           bc_params(2) = relax_V
           bc_params(3) = relax_T
           bc_params(4) = beta
           !
           y = x(1)-(probhi(1)/2.0d0)-0.0025d0 
           u(2) = 24.0d0*u_ref*((sin((M_PI/2.0d0)*(y/0.0025d0)))**2.0d0)
           !u(1) = 24.0d0*u_ref*y*(0.5d0-y)
           u(1) = 0.0d0
           u(3) = 0.0d0
           eos_state % massfrac(1) = 1.d0
           eos_state % p = p_ref
           eos_state % T = T_ref
           call eos_tp(eos_state)
         
         end if
    
      end if
    end if

! at hi Y
    if (dir == 2) then
      if (sgn == -1) then
      
       sigma_out = 0.28d0
       beta = 0.60d0

       ! Set outflow pressure
       bc_type = Outflow
       bc_params(1) = sigma_out
       bc_params(4) = beta

       eos_state % massfrac(1) = 1.d0
       eos_state % p = p_ref
       eos_state % T = T_ref
       call eos_tp(eos_state)

      end if
    end if
    
! at low X
    if (dir == 1) then
      if (sgn == 1) then
      
           relax_U = 0.50d0
           relax_V = 0.00d0
           relax_T = - 0.d0 !relax_V
           beta = 1.d0
           
           bc_type = Inflow
           bc_params(1) = relax_U
           bc_params(2) = relax_V
           bc_params(3) = relax_T
           bc_params(4) = beta
           u(1) = 2.0d0
           u(2) = 0.0d0 
           u(3) = 0.0d0
           eos_state % massfrac(1) = 1.d0
           eos_state % p = p_ref
           eos_state % T = T_ref
           call eos_tp(eos_state)
         
         !end if
    
      end if
    end if

! at hi X
    if (dir == 1) then
      if (sgn == -1) then
      
       sigma_out = 0.28d0
       beta = 0.6d0

       ! Set outflow pressure
       bc_type = Outflow
       bc_params(1) = sigma_out
       bc_params(4) = beta

       eos_state % massfrac(1) = 1.d0
       eos_state % p = p_ref
       eos_state % T = T_ref
       call eos_tp(eos_state)

      end if
    end if


    if (rho_only .EQV. .TRUE. ) then

       u_ext(1) = eos_state % rho

    else
 
       u_ext(UFS:UFS+nspec-1) = eos_state % massfrac * eos_state % rho
       u_ext(URHO)               = eos_state % rho
       u_ext(UMX)                = eos_state % rho  *  u(1)
       u_ext(UMY)                = eos_state % rho  *  u(2)
       u_ext(UMZ)                = eos_state % rho  *  u(3)
       u_ext(UTEMP)              = eos_state % T
       u_ext(UEINT)              = eos_state % rho  *   eos_state % e
       u_ext(UEDEN)              = eos_state % rho  *  (eos_state % e + 0.5d0 * (u(1)**2 + u(2)**2) + u(3)**2)

    endif
    
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

    double precision :: x(3)
    integer :: i, j, k, n
    logical rho_only

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
  end subroutine pc_reactfill

end module bc_fill_module

