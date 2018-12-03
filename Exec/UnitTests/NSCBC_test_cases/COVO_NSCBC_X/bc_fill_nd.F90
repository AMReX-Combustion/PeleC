module bc_fill_module

  implicit none

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
                call bcnormal_nscbc(x,adv(domlo(1),j,k,:),adv(i,j,k,:),1,+1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(domhi(1),j,k,:),adv(i,j,k,:),1,-1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(i,domlo(2),k,:),adv(i,j,k,:),2,+1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(i,domhi(2),k,:),adv(i,j,k,:),2,-1,bc_type,bc_params,rho_only)
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
    !            call bcnormal_nscbc(x,adv(i,j,domlo(3),:),adv(i,j,k,:),3,+1,bc_type,bc_params,rho_only)
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
    !            call bcnormal_nscbc(x,adv(i,j,domhi(3),:),adv(i,j,k,:),3,-1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(domlo(1),j,k),adv(i,j,k),1,+1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(domhi(1),j,k),adv(i,j,k),1,-1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(i,domlo(2),k),adv(i,j,k),2,+1,bc_type,bc_params,rho_only)
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
                call bcnormal_nscbc(x,adv(i,domhi(2),k),adv(i,j,k),2,-1,bc_type,bc_params,rho_only)
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
    !            call bcnormal_nscbc(x,adv(i,j,domlo(3)),adv(i,j,k),3,+1,bc_type,bc_params,rho_only)
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
    !            call bcnormal_nscbc(x,adv(i,j,domhi(3)),adv(i,j,k),3,-1,bc_type,bc_params,rho_only)
    !         end do
    !      end do
    !   end do
    !end if
    
  end subroutine pc_denfill





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

