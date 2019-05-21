module gc_nscbc_mod

  private
  public impose_NSCBC

contains

!------------------------------
! Imposing Ghost-Cells Navier-Stokes Characteristic BCs.
! For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp. 3399-3408, 2017. 
!
! Note that for the corner treatment, we depart from the AIAA paper, because
! we found out that the corner coupling method was superfluous and that providing
! transverse terms computed from one-sided derivative do the job.
!
!------------------------------

 subroutine impose_NSCBC(lo, hi, domlo, domhi, &
                         uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                         q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                         qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3, &
                         x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_l3, x_bcMask_h1, x_bcMask_h2, x_bcMask_h3, &
                         y_bcMask, y_bcMask_l1, y_bcMask_l2, y_bcMask_l3, y_bcMask_h1, y_bcMask_h2, y_bcMask_h3, &
                         z_bcMask, z_bcMask_l1, z_bcMask_l2, z_bcMask_l3, z_bcMask_h1, z_bcMask_h2, z_bcMask_h3, &
                         flag_nscbc_isAnyPerio, flag_nscbc_perio, &
                         time,delta,dt,verbose) bind(C, name="impose_NSCBC")
    
  use amrex_fort_module
  use network, only : nspec
  use eos_module
  use fundamental_constants_module, only: k_B, n_A

  use amrex_constants_module
  use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, UserBC
    
  use meth_params_module, only : NVAR, NQAUX,QVAR
  use bc_fill_module, only: bcnormal
 
  implicit none
    
  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) :: domlo(3), domhi(3)
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_l3, x_bcMask_h1, x_bcMask_h2, x_bcMask_h3
  integer, intent(in) :: y_bcMask_l1, y_bcMask_l2, y_bcMask_l3, y_bcMask_h1, y_bcMask_h2, y_bcMask_h3
  integer, intent(in) :: z_bcMask_l1, z_bcMask_l2, z_bcMask_l3, z_bcMask_h1, z_bcMask_h2, z_bcMask_h3
  integer, intent(in) :: flag_nscbc_isAnyPerio
  integer, intent(in) :: flag_nscbc_perio(3)
  
  double precision, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  double precision, intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  double precision, intent(inout) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
  integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2,x_bcMask_l3:x_bcMask_h3)
  integer, intent(inout) :: y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2,y_bcMask_l3:y_bcMask_h3)
  integer, intent(inout) :: z_bcMask(z_bcMask_l1:z_bcMask_h1,z_bcMask_l2:z_bcMask_h2,z_bcMask_l3:z_bcMask_h3)
  double precision, intent(in) :: delta(3), dt, time

  ! Local
  double precision dx, dy, dz
  double precision x, y, z
  
  double precision :: drhodx, dudx, dvdx, dwdx, dpdx
  double precision :: dpdy, dudy, dvdy, dwdy, drhody
  double precision :: dpdz, dudz, dvdz, dwdz, drhodz
  double precision :: L1, L2, L3, L4, L5
  double precision :: M1, M2, M3, M4, M5
  double precision :: N1, N2, N3, N4, N5
  
  double precision :: T1_X, T2_X, T3_X, T4_X, T5_X
  double precision :: T1_Y, T2_Y, T3_Y, T4_Y, T5_Y
  double precision :: T1_Z, T2_Z, T3_Z, T4_Z, T5_Z

  double precision :: U_dummy(NVAR)
  double precision :: U_ext(NVAR)
  double precision, parameter :: small = 1.d-8
  
  integer          :: q_lo(3), q_hi(3)
  integer          :: uin_lo(3),  uin_hi(3)
  integer          :: i, j, k
  integer          :: bc_type, x_bc_type, y_bc_type, z_bc_type
  integer          :: x_isign, y_isign, z_isign, x_idx_Mask, y_idx_Mask, z_idx_Mask
  integer          :: test_keyword_x, test_keyword_y, test_keyword_z
  double precision :: bc_params(6), x_bc_params(6), y_bc_params(6), z_bc_params(6)
  double precision :: bc_target(5), x_bc_target(5), y_bc_target(5), z_bc_target(5)
  
  type (eos_t) :: eos_state
  call build(eos_state)

  q_lo = [q_l1, q_l2, q_l3]
  q_hi = [q_h1, q_h2, q_h3]
  
  uin_lo  = [uin_l1, uin_l2, uin_l3]
  uin_hi  = [uin_h1, uin_h2, uin_h3]
 
  dx = delta(1)
  dy = delta(2)
  dz = delta(3)
  
  x_bcMask(:,:,:) = 0
  y_bcMask(:,:,:) = 0
  z_bcMask(:,:,:) = 0

  if ( flag_nscbc_isAnyPerio == 0) then

 !--------------------------------------------------------------------------   
 ! corners
 !--------------------------------------------------------------------------

   if       (((q_hi(1) > domhi(1)) .or. (q_lo(1) < domlo(1))) &
       .and. ((q_hi(2) > domhi(2)) .or. (q_lo(2) < domlo(2))) &
       .and. ((q_hi(3) > domhi(3)) .or. (q_lo(3) < domlo(3)))) then

    if (q_hi(1) > domhi(1)) then
      test_keyword_x = physbc_hi(1)
      i = domhi(1)
      x_isign = -1
      x_idx_Mask = i+1
    elseif (q_lo(1) < domlo(1)) then
      test_keyword_x = physbc_lo(1)
      i = domlo(1)
      x_isign = 1
      x_idx_Mask = i
    endif
 
    if (q_hi(2) > domhi(2)) then
      test_keyword_y = physbc_hi(2)
      j = domhi(2)
      y_isign = -1
      y_idx_Mask = j+1
    elseif (q_lo(2) < domlo(2)) then
      test_keyword_y = physbc_lo(2)
      j = domlo(2)
      y_isign = 1
      y_idx_Mask = j
    endif
 
    if (q_hi(3) > domhi(3)) then
      test_keyword_z = physbc_hi(3)
      k = domhi(3)
      z_isign = -1
      z_idx_Mask = k+1
    elseif (q_lo(3) < domlo(3)) then
      test_keyword_z = physbc_lo(3)
      k = domlo(3)
      z_isign = 1
      z_idx_Mask = k
    endif

    x   = (dble(i)+HALF)*dx
    y   = (dble(j)+HALF)*dy
    z   = (dble(k)+HALF)*dz
  
    ! Normal derivative along x
    call normal_derivative(i, j, k, 1, x_isign, dx, &
                           dpdx, dudx, dvdx, dwdx, drhodx, &
                           q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
                            
    ! Normal derivative along y
    call normal_derivative(i, j, k, 2, y_isign, dy, &
                           dpdy, dudy, dvdy, dwdy, drhody, &
                           q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
  
    ! Normal derivative along z
    call normal_derivative(i, j, k, 3, z_isign, dz, &
                           dpdz, dudz, dvdz, dwdz, drhodz, &
                           q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
                                                         
    ! Compute transverse terms for X
    call compute_transverse_terms(i, j, k, 1,  &
                                  T1_X, T2_X, T3_X, T4_X, T5_X, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                                  qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                                 
    ! Compute transverse terms for X
    call compute_transverse_terms(i, j, k, 2,  &
                                  T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                                  qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                                 
    ! Compute transverse terms for X
    call compute_transverse_terms(i, j, k, 3,  &
                                  T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                                  qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

    ! Calling user target BC values
    ! x face
    if (test_keyword_x == UserBC) then
      call bcnormal([x,y,z],U_dummy,U_ext,1,x_isign,time,x_bc_type,x_bc_params,x_bc_target)
    else
      x_bc_type = test_keyword_x
    endif
    x_bcMask(x_idx_Mask,j,k) = x_bc_type
    
    ! y face
    if (test_keyword_y == UserBC) then
      call bcnormal([x,y,z],U_dummy,U_ext,2,y_isign,time,y_bc_type,y_bc_params,y_bc_target)
    else
      y_bc_type = test_keyword_y
    endif
    y_bcMask(i,y_idx_Mask,k) = y_bc_type
   
    ! z face
    if (test_keyword_z == UserBC) then
      call bcnormal([x,y,z],U_dummy,U_ext,3,z_isign,time,z_bc_type,z_bc_params,z_bc_target)
    else
      z_bc_type = test_keyword_z
    endif
    z_bcMask(i,j,z_idx_Mask) = z_bc_type
    
    ! Computing the LODI system waves along X
    call compute_waves(i, j, k, 1, x_isign, &
                       x_bc_type, x_bc_params, x_bc_target, &
                       T1_X, T2_X, T3_X, T4_X, T5_X, &
                       L1, L2, L3, L4, L5, &
                       dpdx, dudx, dvdx, dwdx, drhodx, &
                       q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                       qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
   
    ! Computing the LODI system waves along Y
    call compute_waves(i, j, k, 2, y_isign, &
                       y_bc_type, y_bc_params, y_bc_target, &
                       T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                       M1, M2, M3, M4, M5, &
                       dpdy, dudy, dvdy, dwdy, drhody, &
                       q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                       qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
    
    ! Computing the LODI system waves along Z
    call compute_waves(i, j, k, 3, z_isign, &
                       z_bc_type, z_bc_params, z_bc_target, &
                       T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                       N1, N2, N3, N4, N5, &
                       dpdz, dudz, dvdz, dwdz, drhodz, &
                       q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                       qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

    ! Recomputing ghost-cells values with the LODI waves along X
    call update_ghost_cells(i, j, k, x_bc_type, 1, x_isign, dx, &
                            domlo, domhi, &
                            L1, L2, L3, L4, L5, &
                            uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                            q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                            qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                   
    ! Recomputing ghost-cells values with the LODI waves along Y
    call update_ghost_cells(i, j, k, y_bc_type, 2, y_isign, dy, &
                            domlo, domhi, &
                            M1, M2, M3, M4, M5, &
                            uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                            q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                            qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                           
    ! Recomputing ghost-cells values with the LODI waves along Y
    call update_ghost_cells(i, j, k, z_bc_type, 3, z_isign, dz, &
                            domlo, domhi, &
                            N1, N2, N3, N4, N5, &
                            uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                            q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                            qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                           


 endif
endif ! flag_nscbc_isAnyPerio ) 

 !--------------------------------------------------------------------------   
 ! lower X
 !--------------------------------------------------------------------------

 if ((q_lo(1) < domlo(1)) .and. (physbc_lo(1) == UserBC)) then

   i = domlo(1)
   x   = (dble(i)+HALF)*dx

   do j = q_lo(2)+1,q_hi(2)-1

     y   = (dble(j)+HALF)*dy
     
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
   
     do k = q_lo(3)+1,q_hi(3)-1

       z   = (dble(k)+HALF)*dz
     
       if ( flag_nscbc_isAnyPerio == 0) then
         if ((k== domlo(3)) .or. (k == domhi(3))) cycle !Doing that to avoid ghost cells already filled by corners
       endif
       
       ! Normal derivative along x
       call normal_derivative(i, j, k, 1, 1, dx, &
                              dpdx, dudx, dvdx, dwdx, drhodx, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
 
       ! Tangential derivative along y 
       call tangential_derivative(i, j, k, 2, dy, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
                               
       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 3, dz, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 1,  &
                               T1_X, T2_X, T3_X, T4_X, T5_X, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
            
       ! Calling user target BC values 
       call bcnormal([x,y,z],U_dummy,U_ext,1,1,time,bc_type,bc_params,bc_target)
     
       ! Filling bcMask with specific user defined BC type
       if ((j < q_lo(2)+3) .or. (j > q_hi(2)-3) .or. (k < q_lo(3)+3) .or. (k > q_hi(3)-3)) then
         continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
       else
         x_bcMask(i,j,k) = bc_type
       endif

       ! Computing the LODI system waves
       call compute_waves(i, j, k, 1, 1, &
                          bc_type, bc_params, bc_target, &
                          T1_X, T2_X, T3_X, T4_X, T5_X, &
                          L1, L2, L3, L4, L5, &
                          dpdx, dudx, dvdx, dwdx, drhodx, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Recomputing ghost-cells values with the LODI waves
       call update_ghost_cells(i, j, k, bc_type, 1, 1, dx, &
                               domlo, domhi, &
                               L1, L2, L3, L4, L5, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                                    
     enddo
   enddo
 end if
  
 !--------------------------------------------------------------------------   
 ! upper X
 !--------------------------------------------------------------------------
 
 if ((q_hi(1) > domhi(1)) .and. (physbc_hi(1) == UserBC)) then

   i = domhi(1)
   x   = (dble(i)+HALF)*dx
      
   do j = q_lo(2)+1,q_hi(2)-1
   
     y   = (dble(j)+HALF)*dy
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners 
     endif
     
     do k = q_lo(3)+1,q_hi(3)-1

       z   = (dble(k)+HALF)*dz

       if ( flag_nscbc_isAnyPerio == 0) then
         if ((k== domlo(3)) .or. (k == domhi(3))) cycle !Doing that to avoid ghost cells already filled by corners
       endif

     
       ! Normal derivative along x
       call normal_derivative(i, j, k, 1, -1, dx, &
                              dpdx, dudx, dvdx, dwdx, drhodx, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
       
       ! Tangential derivative along y 
       call tangential_derivative(i, j, k, 2, dy, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 3, dz, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
  
       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 1,  &
                               T1_X, T2_X, T3_X, T4_X, T5_X, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Filling bcMask with specific user defined BC type 
       call bcnormal([x,y,z],U_dummy,U_ext,1,-1,time,bc_type,bc_params,bc_target)
       if ((j < q_lo(2)+3) .or. (j > q_hi(2)-3) .or. (k < q_lo(3)+3) .or. (k > q_hi(3)-3)) then
         continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
       else
         x_bcMask(i+1,j,k) = bc_type
       endif
 
       ! Computing the LODI system waves
       call compute_waves(i, j, k, 1, -1, &
                          bc_type, bc_params, bc_target, &
                          T1_X, T2_X, T3_X, T4_X, T5_X, &
                          L1, L2, L3, L4, L5, &
                          dpdx, dudx, dvdx, dwdx, drhodx, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Recomputing ghost-cells values with the LODI waves
       call update_ghost_cells(i, j, k, bc_type, 1, -1, dx, &
                               domlo, domhi, &
                               L1, L2, L3, L4, L5, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     enddo
   enddo
 endif
 
 !--------------------------------------------------------------------------   
 ! lower Y
 !--------------------------------------------------------------------------
 
 if ((q_lo(2) < domlo(2)) .and. (physbc_lo(2) == UserBC)) then
 
   j = domlo(2)
   y   = (dble(j)+HALF)*dy
      
   do i = q_lo(1)+1,q_hi(1)-1
   
     x   = (dble(i)+HALF)*dx
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     do k = q_lo(3)+1,q_hi(3)-1

       z   = (dble(k)+HALF)*dz

       if ( flag_nscbc_isAnyPerio == 0) then
         if ((k== domlo(3)) .or. (k == domhi(3))) cycle !Doing that to avoid ghost cells already filled by corners
       endif
      
       ! Normal derivative along y
       call normal_derivative(i, j, k, 2, 1, dy, &
                              dpdy, dudy, dvdy, dwdy, drhody, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along x 
       call tangential_derivative(i, j, k, 1, dx, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 3, dz, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 2,  &
                               T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3) 
                               
       ! Filling bcMask with specific user defined BC type
       call bcnormal([x,y,z],U_dummy,U_ext,2,1,time,bc_type,bc_params,bc_target)
       if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3) .or. (k < q_lo(3)+3) .or. (k > q_hi(3)-3)) then
         continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
       else
         y_bcMask(i,j,k) = bc_type
       endif
     
       ! Computing the LODI system waves
       call compute_waves(i, j, k, 2, 1, &
                          bc_type, bc_params, bc_target, &
                          T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                          L1, L2, L3, L4, L5, &
                          dpdy, dudy, dvdy, dwdy, drhody, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Recomputing ghost-cells values with the LODI waves
       call update_ghost_cells(i, j, k, bc_type, 2, 1, dy, &
                               domlo, domhi, &
                               L1, L2, L3, L4, L5, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     enddo
   enddo
 end if
     
!--------------------------------------------------------------------------   
! upper Y
!--------------------------------------------------------------------------

 if ((q_hi(2) > domhi(2)) .and. (physbc_hi(2) == UserBC)) then
 
   j = domhi(2)
   y   = (dble(j)+HALF)*dy
      
   do i = q_lo(1)+1,q_hi(1)-1
     
     x   = (dble(i)+HALF)*dx
     
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     do k = q_lo(3)+1,q_hi(3)-1

       z   = (dble(k)+HALF)*dz

       if ( flag_nscbc_isAnyPerio == 0) then
         if ((k== domlo(3)) .or. (k == domhi(3))) cycle !Doing that to avoid ghost cells already filled by corners
       endif
     
        
       ! Normal derivative along y
       call normal_derivative(i, j, k, 2, -1, dy, &
                              dpdy, dudy, dvdy, dwdy, drhody, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along x 
       call tangential_derivative(i, j, k, 1, dx, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 3, dz, &
                                  dpdz, dudz, dvdz, dwdz, drhodz, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 2,  &
                               T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                              
     ! Filling bcMask with specific user defined BC type 
     call bcnormal([x,y,z],U_dummy,U_ext,2,-1,time,bc_type,bc_params,bc_target)
     if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3) .or. (k < q_lo(3)+3) .or. (k > q_hi(3)-3)) then
       continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
     else
       y_bcMask(i,j+1,k) = bc_type
     endif
    
     ! Computing the LODI system waves
     call compute_waves(i, j, k, 2, -1, &
                          bc_type, bc_params, bc_target, &
                          T1_Y, T2_Y, T3_Y, T4_Y, T5_Y, &
                          L1, L2, L3, L4, L5, &
                          dpdy, dudy, dvdy, dwdy, drhody, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     ! Recomputing ghost-cells values with the LODI waves
     call update_ghost_cells(i, j, k, bc_type, 2, -1, dy, &
                             domlo, domhi, &
                             L1, L2, L3, L4, L5, &
                             uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                             q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                             qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     enddo
   enddo
end if


 !--------------------------------------------------------------------------   
 ! lower Z
 !--------------------------------------------------------------------------
 
 if ((q_lo(3) < domlo(3)) .and. (physbc_lo(3) == UserBC)) then
 
   k = domlo(3)
   z   = (dble(k)+HALF)*dz
      
   do i = q_lo(1)+1,q_hi(1)-1
   
     x   = (dble(i)+HALF)*dx
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     do j = q_lo(2)+1,q_hi(2)-1

       y   = (dble(j)+HALF)*dy

       if ( flag_nscbc_isAnyPerio == 0) then
         if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners
       endif
      
       ! Normal derivative along y
       call normal_derivative(i, j, k, 3, 1, dz, &
                              dpdz, dudz, dvdz, dwdz, drhodz, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along x 
       call tangential_derivative(i, j, k, 1, dx, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 2, dy, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 3,  &
                               T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3) 
                               
       ! Filling bcMask with specific user defined BC type
       call bcnormal([x,y,z],U_dummy,U_ext,3,1,time,bc_type,bc_params,bc_target)
       if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3) .or. (j < q_lo(2)+3) .or. (j > q_hi(2)-3)) then
         continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
       else
         z_bcMask(i,j,k) = bc_type
       endif
     
       ! Computing the LODI system waves
       call compute_waves(i, j, k, 3, 1, &
                          bc_type, bc_params, bc_target, &
                          T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                          L1, L2, L3, L4, L5, &
                          dpdz, dudz, dvdz, dwdz, drhodz, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Recomputing ghost-cells values with the LODI waves
       call update_ghost_cells(i, j, k, bc_type, 3, 1, dz, &
                               domlo, domhi, &
                               L1, L2, L3, L4, L5, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     enddo
   enddo
 end if
     
!--------------------------------------------------------------------------   
! upper Z
!--------------------------------------------------------------------------

 if ((q_hi(3) > domhi(3)) .and. (physbc_hi(3) == UserBC)) then
 
   k = domhi(3)
   z   = (dble(k)+HALF)*dz
      
   do i = q_lo(1)+1,q_hi(1)-1
   
     x   = (dble(i)+HALF)*dx
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     do j = q_lo(2)+1,q_hi(2)-1

       y   = (dble(j)+HALF)*dy

       if ( flag_nscbc_isAnyPerio == 0) then
         if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners
       endif
      
       ! Normal derivative along y
       call normal_derivative(i, j, k, 3, -1, dz, &
                              dpdz, dudz, dvdz, dwdz, drhodz, &
                              q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along x 
       call tangential_derivative(i, j, k, 1, dx, &
                                  dpdx, dudx, dvdx, dwdx, drhodx, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Tangential derivative along z
       call tangential_derivative(i, j, k, 2, dy, &
                                  dpdy, dudy, dvdy, dwdy, drhody, &
                                  q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)

       ! Compute transverse terms
       call compute_transverse_terms(i, j, k, 3,  &
                               T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3) 
                               
       ! Filling bcMask with specific user defined BC type
       call bcnormal([x,y,z],U_dummy,U_ext,3,-1,time,bc_type,bc_params,bc_target)
       if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3) .or. (j < q_lo(2)+3) .or. (j > q_hi(2)-3)) then
         continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
       else
         z_bcMask(i,j,k) = bc_type
       endif
     
       ! Computing the LODI system waves
       call compute_waves(i, j, k, 3, -1, &
                          bc_type, bc_params, bc_target, &
                          T1_Z, T2_Z, T3_Z, T4_Z, T5_Z, &
                          L1, L2, L3, L4, L5, &
                          dpdz, dudz, dvdz, dwdz, drhodz, &
                          q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

       ! Recomputing ghost-cells values with the LODI waves
       call update_ghost_cells(i, j, k, bc_type, 3, -1, dz, &
                               domlo, domhi, &
                               L1, L2, L3, L4, L5, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)

     enddo
   enddo
 end if

call destroy(eos_state)

end subroutine impose_NSCBC

!-------------------------------------------------
! Generic routines below
!-------------------------------------------------

  subroutine normal_derivative(i, j, k, idir, isign, delta, &
                               dp, du, dv, dw, drho, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO
  
  integer, intent(in) :: i,j,k,idir,isign
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  double precision, intent(in) :: delta
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  
  double precision, intent(out) :: dp, du, dv, dw, drho
  
  
  if (idir == 1) then
    if (isign == 1) then
    
      ! 2nd order
      dp = ((-3.0d0/2.0d0)*q(i,j,k,QPRES)+2.0d0*q(i+1,j,k,QPRES)-0.5d0*q(i+2,j,k,QPRES))/delta
      du = ((-3.0d0/2.0d0)*q(i,j,k,QU)+2.0d0*q(i+1,j,k,QU)-0.5d0*q(i+2,j,k,QU))/delta
      dv = ((-3.0d0/2.0d0)*q(i,j,k,QV)+2.0d0*q(i+1,j,k,QV)-0.5d0*q(i+2,j,k,QV))/delta
      dw = ((-3.0d0/2.0d0)*q(i,j,k,QW)+2.0d0*q(i+1,j,k,QW)-0.5d0*q(i+2,j,k,QW))/delta
      drho = ((-3.0d0/2.0d0)*q(i,j,k,QRHO)+2.0d0*q(i+1,j,k,QRHO)-0.5d0*q(i+2,j,k,QRHO))/delta
      
    elseif (isign == -1) then
    
       !2nd order
       dp = ((3.0d0/2.0d0)*q(i,j,k,QPRES)-2.0d0*q(i-1,j,k,QPRES)+0.5d0*q(i-2,j,k,QPRES))/delta
       du = ((3.0d0/2.0d0)*q(i,j,k,QU)-2.0d0*q(i-1,j,k,QU)+0.5d0*q(i-2,j,k,QU))/delta
       dv = ((3.0d0/2.0d0)*q(i,j,k,QV)-2.0d0*q(i-1,j,k,QV)+0.5d0*q(i-2,j,k,QV))/delta
       dw = ((3.0d0/2.0d0)*q(i,j,k,QW)-2.0d0*q(i-1,j,k,QW)+0.5d0*q(i-2,j,k,QW))/delta
       drho = ((3.0d0/2.0d0)*q(i,j,k,QRHO)-2.0d0*q(i-1,j,k,QRHO)+0.5d0*q(i-2,j,k,QRHO))/delta
    
    else
      call bl_abort("Problem of isign in impose_NSCBC_3d:normal_derivative")
    end if
    
  elseif (idir == 2) then
  
    if (isign == 1) then
  
      !2nd order
      dp = ((-3.0d0/2.0d0)*q(i,j,k,QPRES)+2.0d0*q(i,j+1,k,QPRES)-0.5d0*q(i,j+2,k,QPRES))/delta
      du = ((-3.0d0/2.0d0)*q(i,j,k,QU)+2.0d0*q(i,j+1,k,QU)-0.5d0*q(i,j+2,k,QU))/delta
      dv = ((-3.0d0/2.0d0)*q(i,j,k,QV)+2.0d0*q(i,j+1,k,QV)-0.5d0*q(i,j+2,k,QV))/delta
      dw = ((-3.0d0/2.0d0)*q(i,j,k,QW)+2.0d0*q(i,j+1,k,QW)-0.5d0*q(i,j+2,k,QW))/delta
      drho = ((-3.0d0/2.0d0)*q(i,j,k,QRHO)+2.0d0*q(i,j+1,k,QRHO)-0.5d0*q(i,j+2,k,QRHO))/delta
     
    elseif (isign == -1) then
    
      !2nd order
      dp = ((3.0d0/2.0d0)*q(i,j,k,QPRES)-2.0d0*q(i,j-1,k,QPRES)+0.5d0*q(i,j-2,k,QPRES))/delta
      du = ((3.0d0/2.0d0)*q(i,j,k,QU)-2.0d0*q(i,j-1,k,QU)+0.5d0*q(i,j-2,k,QU))/delta
      dv = ((3.0d0/2.0d0)*q(i,j,k,QV)-2.0d0*q(i,j-1,k,QV)+0.5d0*q(i,j-2,k,QV))/delta
      dw = ((3.0d0/2.0d0)*q(i,j,k,QW)-2.0d0*q(i,j-1,k,QW)+0.5d0*q(i,j-2,k,QW))/delta
      drho = ((3.0d0/2.0d0)*q(i,j,k,QRHO)-2.0d0*q(i,j-1,k,QRHO)+0.5d0*q(i,j-2,k,QRHO))/delta
  
    else
      call bl_abort("Problem of isign in impose_NSCBC_3d:normal_derivative")
    end if
  
  elseif (idir == 3) then
  
    if (isign == 1) then
  
      !2nd order
      dp = ((-3.0d0/2.0d0)*q(i,j,k,QPRES)+2.0d0*q(i,j,k+1,QPRES)-0.5d0*q(i,j,k+2,QPRES))/delta
      du = ((-3.0d0/2.0d0)*q(i,j,k,QU)+2.0d0*q(i,j,k+1,QU)-0.5d0*q(i,j,k+2,QU))/delta
      dv = ((-3.0d0/2.0d0)*q(i,j,k,QV)+2.0d0*q(i,j,k+1,QV)-0.5d0*q(i,j,k+2,QV))/delta
      dw = ((-3.0d0/2.0d0)*q(i,j,k,QW)+2.0d0*q(i,j,k+1,QW)-0.5d0*q(i,j,k+2,QW))/delta
      drho = ((-3.0d0/2.0d0)*q(i,j,k,QRHO)+2.0d0*q(i,j,k+1,QRHO)-0.5d0*q(i,j,k+2,QRHO))/delta
     
    elseif (isign == -1) then
    
      !2nd order
      dp = ((3.0d0/2.0d0)*q(i,j,k,QPRES)-2.0d0*q(i,j,k-1,QPRES)+0.5d0*q(i,j,k-2,QPRES))/delta
      du = ((3.0d0/2.0d0)*q(i,j,k,QU)-2.0d0*q(i,j,k-1,QU)+0.5d0*q(i,j,k-2,QU))/delta
      dv = ((3.0d0/2.0d0)*q(i,j,k,QV)-2.0d0*q(i,j,k-1,QV)+0.5d0*q(i,j,k-2,QV))/delta
      dw = ((3.0d0/2.0d0)*q(i,j,k,QW)-2.0d0*q(i,j,k-1,QW)+0.5d0*q(i,j,k-2,QW))/delta
      drho = ((3.0d0/2.0d0)*q(i,j,k,QRHO)-2.0d0*q(i,j,k-1,QRHO)+0.5d0*q(i,j,k-2,QRHO))/delta
  
    else
      call bl_abort("Problem of isign in impose_NSCBC_3d:normal_derivative")
    end if
  
  else
      call bl_abort("Problem of idir in impose_NSCBC_3d:normal_derivative")
  end if
  
  end subroutine normal_derivative
  
  !----------------
  
  subroutine tangential_derivative(i, j, k, idir, delta, &
                               dp, du, dv, dw, drho, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO
  
  integer, intent(in) :: i,j,k,idir
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  double precision, intent(in) :: delta
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  
  double precision, intent(out) :: dp, du, dv, dw, drho
  
  ! Warning, idir means the tangential direction, this is different from 2D (sorry)
  if (idir == 1) then

    ! 2nd order Central
    dp = (q(i+1,j,k,QPRES)-q(i-1,j,k,QPRES))/(2.0d0*delta)
    du = (q(i+1,j,k,QU)-q(i-1,j,k,QU))/(2.0d0*delta)
    dv = (q(i+1,j,k,QV)-q(i-1,j,k,QV))/(2.0d0*delta)
    dw = (q(i+1,j,k,QW)-q(i-1,j,k,QW))/(2.0d0*delta)
    drho = (q(i+1,j,k,QRHO)-q(i-1,j,k,QRHO))/(2.0d0*delta)
        
  elseif (idir == 2) then
  
    ! 2nd order Central
    dp = (q(i,j+1,k,QPRES)-q(i,j-1,k,QPRES))/(2.0d0*delta)
    du = (q(i,j+1,k,QU)-q(i,j-1,k,QU))/(2.0d0*delta)
    dv = (q(i,j+1,k,QV)-q(i,j-1,k,QV))/(2.0d0*delta)
    dw = (q(i,j+1,k,QW)-q(i,j-1,k,QW))/(2.0d0*delta)
    drho = (q(i,j+1,k,QRHO)-q(i,j-1,k,QRHO))/(2.0d0*delta)
    
  elseif (idir == 3) then
  
    ! 2nd order Central
    dp = (q(i,j,k+1,QPRES)-q(i,j,k-1,QPRES))/(2.0d0*delta)
    du = (q(i,j,k+1,QU)-q(i,j,k-1,QU))/(2.0d0*delta)
    dv = (q(i,j,k+1,QV)-q(i,j,k-1,QV))/(2.0d0*delta)
    dw = (q(i,j,k+1,QW)-q(i,j,k-1,QW))/(2.0d0*delta)
    drho = (q(i,j,k+1,QRHO)-q(i,j,k-1,QRHO))/(2.0d0*delta)
  
  else
      call bl_abort("Problem of idir in impose_NSCBC_3d:tangential_derivative")
  end if
  
  end subroutine tangential_derivative
  
  !-----------------
  
  subroutine compute_transverse_terms(i, j, k, idir,  &
                               T1, T2, T3, T4, T5, &
                               dpdx, dudx, dvdx, dwdx, drhodx, &
                               dpdy, dudy, dvdy, dwdy, drhody, &
                               dpdz, dudz, dvdz, dwdz, drhodz, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO, NQAUX, QC, QGAMC
  
  integer, intent(in) :: i,j,k,idir
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  double precision, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  double precision, intent(in) :: dpdx, dudx, dvdx, dwdx, drhodx
  double precision, intent(in) :: dpdy, dudy, dvdy, dwdy, drhody
  double precision, intent(in) :: dpdz, dudz, dvdz, dwdz, drhodz
  
  double precision, intent(out) :: T1, T2, T3, T4, T5
  
  double precision :: inv_rho
  
  inv_rho = 1.0d0/q(i,j,k,QRHO)
  
  
  if (idir == 1) then
  
     T1 =  (q(i,j,k,QV)*(dpdy - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dudy)) &
         + (q(i,j,k,QW)*(dpdz - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dudz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dvdy + dwdz))
     
     T2 =  (q(i,j,k,QV)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhody)-dpdy))  &
         + (q(i,j,k,QW)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhodz)-dpdz))
       
     T3 = q(i,j,k,QV)*dvdy + q(i,j,k,QW)*dvdz + dpdy*inv_rho
     
     T4 = q(i,j,k,QV)*dwdy + q(i,j,k,QW)*dwdz + dpdz*inv_rho
     
     
     T5 =  (q(i,j,k,QV)*(dpdy + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dudy)) &
         + (q(i,j,k,QW)*(dpdz + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dudz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dvdy + dwdz))
  
  elseif (idir == 2) then
  
     T1 =  (q(i,j,k,QU)*(dpdx - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dvdx)) &
         + (q(i,j,k,QW)*(dpdz - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dvdz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dudx + dwdz))
     
     T2 = q(i,j,k,QU)*dudx + q(i,j,k,QW)*dudz + dpdx*inv_rho
     
     T3 =  (q(i,j,k,QU)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhodx)-dpdx))  &
         + (q(i,j,k,QW)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhodz)-dpdz))
     
     T4 = q(i,j,k,QU)*dwdx + q(i,j,k,QW)*dwdz + dpdz*inv_rho
     
     T5 =  (q(i,j,k,QU)*(dpdx + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dvdx)) &
         + (q(i,j,k,QW)*(dpdz + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dvdz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dudx + dwdz))
           
  elseif (idir == 3) then
  
     T1 =  (q(i,j,k,QU)*(dpdx - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dwdx)) &
         + (q(i,j,k,QV)*(dpdy - q(i,j,k,QRHO)*qaux(i,j,k,QC)*dwdz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dudx + dvdy))
     
     T2 =  q(i,j,k,QU)*dudx + q(i,j,k,QV)*dudy + dpdx*inv_rho
     
     T3 =  q(i,j,k,QU)*dvdx + q(i,j,k,QV)*dvdy + dpdy*inv_rho
     
     T4 =  (q(i,j,k,QU)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhodx)-dpdx))  &
         + (q(i,j,k,QV)*((qaux(i,j,k,QC)*qaux(i,j,k,QC)*drhody)-dpdy))
     
     T5 =  (q(i,j,k,QU)*(dpdx + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dwdx)) &
         + (q(i,j,k,QV)*(dpdy + q(i,j,k,QRHO)*qaux(i,j,k,QC)*dwdz)) &  
         + (qaux(i,j,k,QGAMC) * q(i,j,k,QPRES)*(dudx + dvdy))
     
  else
      call bl_abort("Problem of idir in impose_NSCBC_2d:compute_transverse_terms")
  end if

  end subroutine compute_transverse_terms
    
  !-------------------------
  
  subroutine update_ghost_cells(i, j, k, bc_type, idir, isign, delta, &
                                domlo, domhi, &
                                L1, L2, L3, L4, L5, &
                                uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                                qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                               
  use eos_module
  use amrex_constants_module, only : ONE
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP,&
                                 UFS, NQAUX, QC, QGAMC, QRSPEC, &
                                 QC, QDPDE, QDPDR, QCSML, QGAMC, &
                                 QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, &
                                 QFS, QFX, QGAME, NHYP
  use prob_params_module, only : SlipWall, NoSlipWall
  
  
  integer, intent(in) :: i,j,k,idir,isign,bc_type
  integer, intent(in) :: domlo(3), domhi(3)
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  double precision, intent(in) :: L1, L2, L3, L4, L5
  double precision, intent(in) :: delta
  double precision, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  double precision, intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  double precision, intent(inout) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    
  integer :: idx_gc1, idx_gc2, idx_gc3, idx_gc4, idx_int1, idx_int2, idx_int3
  integer :: idx_start, idx_end, hop, n, local_index
  double precision :: drho, du, dv, dw, dp, wall_sign
  double precision, parameter :: small = 1.d-8
  
  type (eos_t) :: eos_state
  
  call build(eos_state)
  
  if ((idir == 1) .or. (idir == 2) .or. (idir == 3)) then
    continue
  else
    call bl_abort("Problem of idir in impose_NSCBC_3d:update_ghost_cells")
  end if
  
  if ((isign == 1) .or. (isign == -1)) then
    continue
  else
    call bl_abort("Problem of isign in impose_NSCBC_3d:update_ghost_cells")
  end if
  
  ! Compute new spatial derivative
  
  if (idir == 1) then
    local_index = i
    drho = (L2 + 0.5d0*(L1 + L5))/(qaux(i,j,k,QC)**2.0d0)  
    du   = (L5-L1)/(2.0d0*qaux(i,j,k,QC)*q(i,j,k,QRHO))
    dv   = L3
    dw   = L4
    dp   = 0.5d0*(L1+L5)
  elseif (idir == 2) then
    local_index = j
    drho = (L3 + 0.5d0*(L1 + L5))/(qaux(i,j,k,QC)**2.0d0)
    du   = L2
    dv   = (L5-L1)/(2.0d0*qaux(i,j,k,QC)*q(i,j,k,QRHO))
    dw   = L4
    dp   = 0.5d0*(L1+L5)
  elseif (idir == 3) then
    local_index = k
    drho = (L4 + 0.5d0*(L1 + L5))/(qaux(i,j,k,QC)**2.0d0)
    du   = L2
    dv   = L3
    dw   = (L5-L1)/(2.0d0*qaux(i,j,k,QC)*q(i,j,k,QRHO))
    dp   = 0.5d0*(L1+L5)
  endif
    
  if (isign == 1) then
    idx_gc1 = local_index-1
    idx_gc2 = local_index-2
    idx_gc3 = local_index-3
    idx_gc4 = local_index-4
    idx_int1 = local_index+1
    idx_int2 = local_index+2
    idx_int3 = local_index+3
    idx_start = domlo(idir)-1
    idx_end   = domlo(idir)-4
  elseif (isign == -1) then
    idx_gc1 = local_index+1
    idx_gc2 = local_index+2
    idx_gc3 = local_index+3
    idx_gc4 = local_index+4
    idx_int1 = local_index-1
    idx_int2 = local_index-2
    idx_int3 = local_index-3
    idx_start = domhi(idir)+1
    idx_end   = domhi(idir)+4
  endif
    
  if (idir == 1) then
    
    ! Update ghost cells
    ! 2nd order
    q(idx_gc1,j,k,QU)    = q(idx_int1,j,k,QU) - 2.0d0*delta*du*isign
    q(idx_gc1,j,k,QV)    = q(idx_int1,j,k,QV) - 2.0d0*delta*dv*isign
    q(idx_gc1,j,k,QW)    = q(idx_int1,j,k,QW) - 2.0d0*delta*dw*isign
    q(idx_gc1,j,k,QRHO)  = q(idx_int1,j,k,QRHO)  - 2.0d0*delta*drho*isign
    q(idx_gc1,j,k,QPRES) = q(idx_int1,j,k,QPRES) - 2.0d0*delta*dp*isign
  
    !---------------- 
    q(idx_gc2,j,k,QU)    = -2.0d0*q(idx_int1,j,k,QU) - 3.0d0*q(i,j,k,QU) + 6.0d0*q(idx_gc1,j,k,QU) + 6.0d0*delta*du*isign
    q(idx_gc2,j,k,QV)    = -2.0d0*q(idx_int1,j,k,QV) - 3.0d0*q(i,j,k,QV) + 6.0d0*q(idx_gc1,j,k,QV) + 6.0d0*delta*dv*isign
    q(idx_gc2,j,k,QW)    = -2.0d0*q(idx_int1,j,k,QW) - 3.0d0*q(i,j,k,QW) + 6.0d0*q(idx_gc1,j,k,QW) + 6.0d0*delta*dw*isign
    q(idx_gc2,j,k,QRHO)  = -2.0d0*q(idx_int1,j,k,QRHO) - 3.0d0*q(i,j,k,QRHO) + 6.0d0*q(idx_gc1,j,k,QRHO) + 6.0d0*delta*drho*isign
    q(idx_gc2,j,k,QPRES) = -2.0d0*q(idx_int1,j,k,QPRES) - 3.0d0*q(i,j,k,QPRES) + 6.0d0*q(idx_gc1,j,k,QPRES) + 6.0d0*delta*dp*isign
 
    q(idx_gc3,j,k,QU)    = 3.0d0*q(idx_int1,j,k,QU) +10.0d0*q(i,j,k,QU) - 18.0d0*q(idx_gc1,j,k,QU) &
                         + 6.0d0*q(idx_gc2,j,k,QU) - 12.0d0*delta*du*isign
    q(idx_gc3,j,k,QV)    = 3.0d0*q(idx_int1,j,k,QV) +10.0d0*q(i,j,k,QV) - 18.0d0*q(idx_gc1,j,k,QV) &
                         + 6.0d0*q(idx_gc2,j,k,QV) - 12.0d0*delta*dv*isign
    q(idx_gc3,j,k,QW)    = 3.0d0*q(idx_int1,j,k,QW) +10.0d0*q(i,j,k,QW) - 18.0d0*q(idx_gc1,j,k,QW) &
                         + 6.0d0*q(idx_gc2,j,k,QW) - 12.0d0*delta*dw*isign
    q(idx_gc3,j,k,QRHO)    = 3.0d0*q(idx_int1,j,k,QRHO) +10.0d0*q(i,j,k,QRHO) - 18.0d0*q(idx_gc1,j,k,QRHO) &
                         + 6.0d0*q(idx_gc2,j,k,QRHO) - 12.0d0*delta*drho*isign
    q(idx_gc3,j,k,QPRES)    = 3.0d0*q(idx_int1,j,k,QPRES) +10.0d0*q(i,j,k,QPRES) - 18.0d0*q(idx_gc1,j,k,QPRES) &
                         + 6.0d0*q(idx_gc2,j,k,QPRES) - 12.0d0*delta*dp*isign
 
    q(idx_gc4,j,k,QU)    = -2.0d0*q(idx_int1,j,k,QU) - 13.0d0*q(i,j,k,QU) + 24.0d0*q(idx_gc1,j,k,QU) - 12.0d0*q(idx_gc2,j,k,QU)  &
                    + 4.0d0*q(idx_gc3,j,k,QU) + 12.0d0*delta*du*isign
    q(idx_gc4,j,k,QV)    = -2.0d0*q(idx_int1,j,k,QV) - 13.0d0*q(i,j,k,QV) + 24.0d0*q(idx_gc1,j,k,QV) - 12.0d0*q(idx_gc2,j,k,QV) &
                    + 4.0d0*q(idx_gc3,j,k,QV) + 12.0d0*delta*dv*isign
    q(idx_gc4,j,k,QW)    = -2.0d0*q(idx_int1,j,k,QW) - 13.0d0*q(i,j,k,QW) + 24.0d0*q(idx_gc1,j,k,QW) - 12.0d0*q(idx_gc2,j,k,QW) &
                    + 4.0d0*q(idx_gc3,j,k,QW) + 12.0d0*delta*dw*isign
    q(idx_gc4,j,k,QRHO)  = -2.0d0*q(idx_int1,j,k,QRHO) - 13.0d0*q(i,j,k,QRHO) + 24.0d0*q(idx_gc1,j,k,QRHO) - 12.0d0*q(idx_gc2,j,k,QRHO) &
                    + 4.0d0*q(idx_gc3,j,k,QRHO) + 12.0d0*delta*drho*isign
    q(idx_gc4,j,k,QPRES) = -2.0d0*q(idx_int1,j,k,QPRES) - 13.0d0*q(i,j,k,QPRES) + 24.0d0*q(idx_gc1,j,k,QPRES) - 12.0d0*q(idx_gc2,j,k,QPRES) &
                    + 4.0d0*q(idx_gc3,j,k,QPRES) + 12.0d0*delta*dp*isign
 
    if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
    
      if (bc_type .eq. NoSlipWall) then
        wall_sign = -1.0d0
      else if (bc_type .eq. SlipWall)  then
        wall_sign = 1.0d0
      end if
      
      q(idx_gc1,j,k,QU)    = -q(i,j,k,QU)
      q(idx_gc2,j,k,QU)    = -q(idx_int1,j,k,QU)
      q(idx_gc3,j,k,QU)    = -q(idx_int2,j,k,QU)
      q(idx_gc4,j,k,QU)    = -q(idx_int3,j,k,QU)
 
      q(idx_gc1,j,k,QV)    = wall_sign*q(i,j,k,QV)
      q(idx_gc2,j,k,QV)    = wall_sign*q(idx_int1,j,k,QV)
      q(idx_gc3,j,k,QV)    = wall_sign*q(idx_int2,j,k,QV)
      q(idx_gc4,j,k,QV)    = wall_sign*q(idx_int3,j,k,QV)
      
      q(idx_gc1,j,k,QW)    = wall_sign*q(i,j,k,QW)
      q(idx_gc2,j,k,QW)    = wall_sign*q(idx_int1,j,k,QW)
      q(idx_gc3,j,k,QW)    = wall_sign*q(idx_int2,j,k,QW)
      q(idx_gc4,j,k,QW)    = wall_sign*q(idx_int3,j,k,QW)
      
      q(idx_gc1,j,k,QRHO)  = q(i,j,k,QRHO)
      q(idx_gc2,j,k,QRHO)  = q(idx_int1,j,k,QRHO)
      q(idx_gc3,j,k,QRHO)  = q(idx_int2,j,k,QRHO)
      q(idx_gc4,j,k,QRHO)  = q(idx_int3,j,k,QRHO)
      
      q(idx_gc1,j,k,QPRES)  = q(i,j,k,QPRES)
      q(idx_gc2,j,k,QPRES)  = q(idx_int1,j,k,QPRES)
      q(idx_gc3,j,k,QPRES)  = q(idx_int2,j,k,QPRES)
      q(idx_gc4,j,k,QPRES)  = q(idx_int3,j,k,QPRES)
    
    end if
 
    ! Recompute missing values thanks to EOS
    do hop=idx_start,idx_end,-isign
    
      eos_state % p        = q(hop,j,k,QPRES )
      eos_state % rho      = q(hop,j,k,QRHO  )
      eos_state % massfrac = q(hop,j,k,QFS:QFS+nspec-1)
      eos_state % aux      = q(hop,j,k,QFX:QFX+naux-1)
 
      call eos_rp(eos_state)
      q(hop,j,k,QTEMP)  = eos_state % T
      q(hop,j,k,QREINT) = eos_state % e * q(hop,j,k,QRHO)
      q(hop,j,k,QGAME)  = q(hop,j,k,QPRES) / q(hop,j,k,QREINT) + ONE
      
      qaux(hop,j,k,QDPDR)  = eos_state % dpdr_e
      qaux(hop,j,k,QDPDE)  = eos_state % dpde
      qaux(hop,j,k,QGAMC)  = eos_state % gam1
      qaux(hop,j,k,QC   )  = eos_state % cs
      qaux(hop,j,k,QCSML)  = max(small, small * qaux(hop,j,k,QC))
 
      ! Here the update of the conservative variables uin seems to have only an impact
      ! on the application of artificial viscosity difmag. 
      uin(hop,j,k,URHO )  = eos_state % rho 
      uin(hop,j,k,UMX  )  = q(hop,j,k,QU ) * eos_state % rho 
      uin(hop,j,k,UMY  )  = q(hop,j,k,QV ) * eos_state % rho 
      uin(hop,j,k,UMZ  ) =  q(hop,j,k,QW ) * eos_state % rho
      uin(hop,j,k,UEINT) = eos_state % rho   *  eos_state % e
      uin(hop,j,k,UEDEN) = eos_state % rho  &
         * (eos_state % e + 0.5d0 * (q(hop,j,k,QU)**2 + q(hop,j,k,QV)**2 + q(hop,j,k,QW)**2))
      uin(hop,j,k,UTEMP) = eos_state % T
      do n=1, nspec
         uin(hop,j,k,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
      end do   
      
    enddo
  
  elseif (idir == 2) then
       
    ! Update ghost cells
    ! 2nd order
    q(i,idx_gc1,k,QU)    = q(i,idx_int1,k,QU) - 2.0d0*delta*du*isign
    q(i,idx_gc1,k,QV)    = q(i,idx_int1,k,QV) - 2.0d0*delta*dv*isign
    q(i,idx_gc1,k,QW)    = q(i,idx_int1,k,QW) - 2.0d0*delta*dw*isign
    q(i,idx_gc1,k,QRHO)  = q(i,idx_int1,k,QRHO)  - 2.0d0*delta*drho*isign
    q(i,idx_gc1,k,QPRES) = q(i,idx_int1,k,QPRES) - 2.0d0*delta*dp*isign
   
    !---------------- 
    q(i,idx_gc2,k,QU)    = -2.0d0*q(i,idx_int1,k,QU) - 3.0d0*q(i,j,k,QU) + 6.0d0*q(i,idx_gc1,k,QU) + 6.0d0*delta*du*isign
    q(i,idx_gc2,k,QV)    = -2.0d0*q(i,idx_int1,k,QV) - 3.0d0*q(i,j,k,QV) + 6.0d0*q(i,idx_gc1,k,QV) + 6.0d0*delta*dv*isign
    q(i,idx_gc2,k,QW)    = -2.0d0*q(i,idx_int1,k,QW) - 3.0d0*q(i,j,k,QW) + 6.0d0*q(i,idx_gc1,k,QW) + 6.0d0*delta*dw*isign
    q(i,idx_gc2,k,QRHO)  = -2.0d0*q(i,idx_int1,k,QRHO) - 3.0d0*q(i,j,k,QRHO) + 6.0d0*q(i,idx_gc1,k,QRHO) + 6.0d0*delta*drho*isign
    q(i,idx_gc2,k,QPRES) = -2.0d0*q(i,idx_int1,k,QPRES) - 3.0d0*q(i,j,k,QPRES) + 6.0d0*q(i,idx_gc1,k,QPRES) + 6.0d0*delta*dp*isign
  
    q(i,idx_gc3,k,QU)    = 3.0d0*q(i,idx_int1,k,QU) +10.0d0*q(i,j,k,QU) - 18.0d0*q(i,idx_gc1,k,QU) &
                       + 6.0d0*q(i,idx_gc2,k,QU) - 12.0d0*delta*du*isign
    q(i,idx_gc3,k,QV)    = 3.0d0*q(i,idx_int1,k,QV) +10.0d0*q(i,j,k,QV) - 18.0d0*q(i,idx_gc1,k,QV) &
                       + 6.0d0*q(i,idx_gc2,k,QV) - 12.0d0*delta*dv*isign
    q(i,idx_gc3,k,QW)    = 3.0d0*q(i,idx_int1,k,QW) +10.0d0*q(i,j,k,QW) - 18.0d0*q(i,idx_gc1,k,QW) &
                       + 6.0d0*q(i,idx_gc2,k,QW) - 12.0d0*delta*dw*isign
    q(i,idx_gc3,k,QRHO)  = 3.0d0*q(i,idx_int1,k,QRHO) +10.0d0*q(i,j,k,QRHO) - 18.0d0*q(i,idx_gc1,k,QRHO) &
                       + 6.0d0*q(i,idx_gc2,k,QRHO) - 12.0d0*delta*drho*isign
    q(i,idx_gc3,k,QPRES) = 3.0d0*q(i,idx_int1,k,QPRES) +10.0d0*q(i,j,k,QPRES) - 18.0d0*q(i,idx_gc1,k,QPRES) &
                       + 6.0d0*q(i,idx_gc2,k,QPRES) - 12.0d0*delta*dp*isign
  
    q(i,idx_gc4,k,QU)    = -2.0d0*q(i,idx_int1,k,QU) - 13.0d0*q(i,j,k,QU) + 24.0d0*q(i,idx_gc1,k,QU) - 12.0d0*q(i,idx_gc2,k,QU)  &
                    + 4.0d0*q(i,idx_gc3,k,QU) + 12.0d0*delta*du*isign
    q(i,idx_gc4,k,QV)    = -2.0d0*q(i,idx_int1,k,QV) - 13.0d0*q(i,j,k,QV) + 24.0d0*q(i,idx_gc1,k,QV) - 12.0d0*q(i,idx_gc2,k,QV) &
                    + 4.0d0*q(i,idx_gc3,k,QV) + 12.0d0*delta*dv*isign
    q(i,idx_gc4,k,QW)    = -2.0d0*q(i,idx_int1,k,QW) - 13.0d0*q(i,j,k,QW) + 24.0d0*q(i,idx_gc1,k,QW) - 12.0d0*q(i,idx_gc2,k,QW) &
                    + 4.0d0*q(i,idx_gc3,k,QW) + 12.0d0*delta*dw*isign
    q(i,idx_gc4,k,QRHO)  = -2.0d0*q(i,idx_int1,k,QRHO) - 13.0d0*q(i,j,k,QRHO) + 24.0d0*q(i,idx_gc1,k,QRHO) - 12.0d0*q(i,idx_gc2,k,QRHO) &
                    + 4.0d0*q(i,idx_gc3,k,QRHO) + 12.0d0*delta*drho*isign
    q(i,idx_gc4,k,QPRES) = -2.0d0*q(i,idx_int1,k,QPRES) - 13.0d0*q(i,j,k,QPRES) + 24.0d0*q(i,idx_gc1,k,QPRES) - 12.0d0*q(i,idx_gc2,k,QPRES) &
                    + 4.0d0*q(i,idx_gc3,k,QPRES) + 12.0d0*delta*dp*isign
  
    if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
      if (bc_type .eq. NoSlipWall) then
        wall_sign = -1.0d0
      else if (bc_type .eq. SlipWall)  then
        wall_sign = 1.0d0
      end if
       
      q(i,idx_gc1,k,QU)    = wall_sign*q(i,j,k,QU)
      q(i,idx_gc2,k,QU)    = wall_sign*q(i,idx_int1,k,QU)
      q(i,idx_gc3,k,QU)    = wall_sign*q(i,idx_int2,k,QU)
      q(i,idx_gc4,k,QU)    = wall_sign*q(i,idx_int3,k,QU)
       
      q(i,idx_gc1,k,QV)    = -q(i,j,k,QV)
      q(i,idx_gc2,k,QV)    = -q(i,idx_int1,k,QV)
      q(i,idx_gc3,k,QV)    = -q(i,idx_int2,k,QV)
      q(i,idx_gc4,k,QV)    = -q(i,idx_int3,k,QV)
       
      q(i,idx_gc1,k,QW)    = wall_sign*q(i,j,k,QW)
      q(i,idx_gc2,k,QW)    = wall_sign*q(i,idx_int1,k,QW)
      q(i,idx_gc3,k,QW)    = wall_sign*q(i,idx_int2,k,QW)
      q(i,idx_gc4,k,QW)    = wall_sign*q(i,idx_int3,k,QW)
       
      q(i,idx_gc1,k,QRHO)  = q(i,j,k,QRHO)
      q(i,idx_gc2,k,QRHO)  = q(i,idx_int1,k,QRHO)
      q(i,idx_gc3,k,QRHO)  = q(i,idx_int2,k,QRHO)
      q(i,idx_gc4,k,QRHO)  = q(i,idx_int3,k,QRHO)
       
      q(i,idx_gc1,k,QPRES)  = q(i,j,k,QPRES)
      q(i,idx_gc2,k,QPRES)  = q(i,idx_int1,k,QPRES)
      q(i,idx_gc3,k,QPRES)  = q(i,idx_int2,k,QPRES)
      q(i,idx_gc4,k,QPRES)  = q(i,idx_int3,k,QPRES)
     
    end if
      
    ! Recompute missing values thanks to EOS
    do hop=idx_start,idx_end,-isign
     
      eos_state % p        = q(i,hop,k,QPRES )
      eos_state % rho      = q(i,hop,k,QRHO  )
      eos_state % massfrac = q(i,hop,k,QFS:QFS+nspec-1)
      eos_state % aux      = q(i,hop,k,QFX:QFX+naux-1)
  
      call eos_rp(eos_state)
      q(i,hop,k,QTEMP)  = eos_state % T
      q(i,hop,k,QREINT) = eos_state % e * q(i,hop,k,QRHO)
      q(i,hop,k,QGAME)  = q(i,hop,k,QPRES) / q(i,hop,k,QREINT) + ONE
       
      qaux(i,hop,k,QDPDR)  = eos_state % dpdr_e
      qaux(i,hop,k,QDPDE)  = eos_state % dpde
      qaux(i,hop,k,QGAMC)  = eos_state % gam1
      qaux(i,hop,k,QC   )  = eos_state % cs
      qaux(i,hop,k,QCSML)  = max(small, small * qaux(i,hop,k,QC))
  
      ! Here the update of the conservative variables uin seems to have only an impact
      ! on the application of artificial viscosity difmag. 
      uin(i,hop,k,URHO )  = eos_state % rho 
      uin(i,hop,k,UMX  )  = q(i,hop,k,QU ) * eos_state % rho 
      uin(i,hop,k,UMY  )  = q(i,hop,k,QV ) * eos_state % rho 
      uin(i,hop,k,UMZ  )  = q(i,hop,k,QW ) * eos_state % rho 
      uin(i,hop,k,UEINT) = eos_state % rho   *  eos_state % e
      uin(i,hop,k,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (q(i,hop,k,QU)**2 + q(i,hop,k,QV)**2 + q(i,hop,k,QW)**2))
      uin(i,hop,k,UTEMP) = eos_state % T
      do n=1, nspec
         uin(i,hop,k,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
      end do   
       
    enddo
    
  elseif (idir == 3) then
       
    ! Update ghost cells
    ! 2nd order
    q(i,j,idx_gc1,QU)    = q(i,j,idx_int1,QU) - 2.0d0*delta*du*isign
    q(i,j,idx_gc1,QV)    = q(i,j,idx_int1,QV) - 2.0d0*delta*dv*isign
    q(i,j,idx_gc1,QW)    = q(i,j,idx_int1,QW) - 2.0d0*delta*dw*isign
    q(i,j,idx_gc1,QRHO)  = q(i,j,idx_int1,QRHO)  - 2.0d0*delta*drho*isign
    q(i,j,idx_gc1,QPRES) = q(i,j,idx_int1,QPRES) - 2.0d0*delta*dp*isign
   
    !---------------- 
    q(i,j,idx_gc2,QU)    = -2.0d0*q(i,j,idx_int1,QU) - 3.0d0*q(i,j,k,QU) + 6.0d0*q(i,j,idx_gc1,QU) + 6.0d0*delta*du*isign
    q(i,j,idx_gc2,QV)    = -2.0d0*q(i,j,idx_int1,QV) - 3.0d0*q(i,j,k,QV) + 6.0d0*q(i,j,idx_gc1,QV) + 6.0d0*delta*dv*isign
    q(i,j,idx_gc2,QW)    = -2.0d0*q(i,j,idx_int1,QW) - 3.0d0*q(i,j,k,QW) + 6.0d0*q(i,j,idx_gc1,QW) + 6.0d0*delta*dw*isign
    q(i,j,idx_gc2,QRHO)  = -2.0d0*q(i,j,idx_int1,QRHO) - 3.0d0*q(i,j,k,QRHO) + 6.0d0*q(i,j,idx_gc1,QRHO) + 6.0d0*delta*drho*isign
    q(i,j,idx_gc2,QPRES) = -2.0d0*q(i,j,idx_int1,QPRES) - 3.0d0*q(i,j,k,QPRES) + 6.0d0*q(i,j,idx_gc1,QPRES) + 6.0d0*delta*dp*isign
  
    q(i,j,idx_gc3,QU)    = 3.0d0*q(i,j,idx_int1,QU) +10.0d0*q(i,j,k,QU) - 18.0d0*q(i,j,idx_gc1,QU) &
                       + 6.0d0*q(i,j,idx_gc2,QU) - 12.0d0*delta*du*isign
    q(i,j,idx_gc3,QV)    = 3.0d0*q(i,j,idx_int1,QV) +10.0d0*q(i,j,k,QV) - 18.0d0*q(i,j,idx_gc1,QV) &
                       + 6.0d0*q(i,j,idx_gc2,QV) - 12.0d0*delta*dv*isign
    q(i,j,idx_gc3,QW)    = 3.0d0*q(i,j,idx_int1,QW) +10.0d0*q(i,j,k,QW) - 18.0d0*q(i,j,idx_gc1,QW) &
                       + 6.0d0*q(i,j,idx_gc2,QW) - 12.0d0*delta*dw*isign
    q(i,j,idx_gc3,QRHO)  = 3.0d0*q(i,j,idx_int1,QRHO) +10.0d0*q(i,j,k,QRHO) - 18.0d0*q(i,j,idx_gc1,QRHO) &
                       + 6.0d0*q(i,j,idx_gc2,QRHO) - 12.0d0*delta*drho*isign
    q(i,j,idx_gc3,QPRES) = 3.0d0*q(i,j,idx_int1,QPRES) +10.0d0*q(i,j,k,QPRES) - 18.0d0*q(i,j,idx_gc1,QPRES) &
                       + 6.0d0*q(i,j,idx_gc2,QPRES) - 12.0d0*delta*dp*isign
  
    q(i,j,idx_gc4,QU)    = -2.0d0*q(i,j,idx_int1,QU) - 13.0d0*q(i,j,k,QU) + 24.0d0*q(i,j,idx_gc1,QU) - 12.0d0*q(i,j,idx_gc2,QU)  &
                    + 4.0d0*q(i,j,idx_gc3,QU) + 12.0d0*delta*du*isign
    q(i,j,idx_gc4,QV)    = -2.0d0*q(i,j,idx_int1,QV) - 13.0d0*q(i,j,k,QV) + 24.0d0*q(i,j,idx_gc1,QV) - 12.0d0*q(i,j,idx_gc2,QV) &
                    + 4.0d0*q(i,j,idx_gc3,QV) + 12.0d0*delta*dv*isign
    q(i,j,idx_gc4,QW)    = -2.0d0*q(i,j,idx_int1,QW) - 13.0d0*q(i,j,k,QW) + 24.0d0*q(i,j,idx_gc1,QW) - 12.0d0*q(i,j,idx_gc2,QW) &
                    + 4.0d0*q(i,j,idx_gc3,QW) + 12.0d0*delta*dw*isign
    q(i,j,idx_gc4,QRHO)  = -2.0d0*q(i,j,idx_int1,QRHO) - 13.0d0*q(i,j,k,QRHO) + 24.0d0*q(i,j,idx_gc1,QRHO) - 12.0d0*q(i,j,idx_gc2,QRHO) &
                    + 4.0d0*q(i,j,idx_gc3,QRHO) + 12.0d0*delta*drho*isign
    q(i,j,idx_gc4,QPRES) = -2.0d0*q(i,j,idx_int1,QPRES) - 13.0d0*q(i,j,k,QPRES) + 24.0d0*q(i,j,idx_gc1,QPRES) - 12.0d0*q(i,j,idx_gc2,QPRES) &
                    + 4.0d0*q(i,j,idx_gc3,QPRES) + 12.0d0*delta*dp*isign
  
    if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
      if (bc_type .eq. NoSlipWall) then
        wall_sign = -1.0d0
      else if (bc_type .eq. SlipWall)  then
        wall_sign = 1.0d0
      end if
       
      q(i,j,idx_gc1,QU)    = wall_sign*q(i,j,k,QU)
      q(i,j,idx_gc2,QU)    = wall_sign*q(i,j,idx_int1,QU)
      q(i,j,idx_gc3,QU)    = wall_sign*q(i,j,idx_int2,QU)
      q(i,j,idx_gc4,QU)    = wall_sign*q(i,j,idx_int3,QU)
              
      q(i,j,idx_gc1,QV)    = wall_sign*q(i,j,k,QV)
      q(i,j,idx_gc2,QV)    = wall_sign*q(i,j,idx_int1,QV)
      q(i,j,idx_gc3,QV)    = wall_sign*q(i,j,idx_int2,QV)
      q(i,j,idx_gc4,QV)    = wall_sign*q(i,j,idx_int3,QV)
       
      q(i,j,idx_gc1,QW)    = -q(i,j,k,QW)
      q(i,j,idx_gc2,QW)    = -q(i,j,idx_int1,QW)
      q(i,j,idx_gc3,QW)    = -q(i,j,idx_int2,QW)
      q(i,j,idx_gc4,QW)    = -q(i,j,idx_int3,QW)
       
      q(i,j,idx_gc1,QRHO)  = q(i,j,k,QRHO)
      q(i,j,idx_gc2,QRHO)  = q(i,j,idx_int1,QRHO)
      q(i,j,idx_gc3,QRHO)  = q(i,j,idx_int2,QRHO)
      q(i,j,idx_gc4,QRHO)  = q(i,j,idx_int3,QRHO)
       
      q(i,j,idx_gc1,QPRES)  = q(i,j,k,QPRES)
      q(i,j,idx_gc2,QPRES)  = q(i,j,idx_int1,QPRES)
      q(i,j,idx_gc3,QPRES)  = q(i,j,idx_int2,QPRES)
      q(i,j,idx_gc4,QPRES)  = q(i,j,idx_int3,QPRES)
     
    end if
      
    ! Recompute missing values thanks to EOS
    do hop=idx_start,idx_end,-isign
     
      eos_state % p        = q(i,j,hop,QPRES )
      eos_state % rho      = q(i,j,hop,QRHO  )
      eos_state % massfrac = q(i,j,hop,QFS:QFS+nspec-1)
      eos_state % aux      = q(i,j,hop,QFX:QFX+naux-1)
  
      call eos_rp(eos_state)
      q(i,j,hop,QTEMP)  = eos_state % T
      q(i,j,hop,QREINT) = eos_state % e * q(i,j,hop,QRHO)
      q(i,j,hop,QGAME)  = q(i,j,hop,QPRES) / q(i,j,hop,QREINT) + ONE
       
      qaux(i,j,hop,QDPDR)  = eos_state % dpdr_e
      qaux(i,j,hop,QDPDE)  = eos_state % dpde
      qaux(i,j,hop,QGAMC)  = eos_state % gam1
      qaux(i,j,hop,QC   )  = eos_state % cs
      qaux(i,j,hop,QCSML)  = max(small, small * qaux(i,j,hop,QC))
  
      ! Here the update of the conservative variables uin seems to have only an impact
      ! on the application of artificial viscosity difmag. 
      uin(i,j,hop,URHO )  = eos_state % rho 
      uin(i,j,hop,UMX  )  = q(i,j,hop,QU ) * eos_state % rho 
      uin(i,j,hop,UMY  )  = q(i,j,hop,QV ) * eos_state % rho 
      uin(i,j,hop,UMZ  )  = q(i,j,hop,QW ) * eos_state % rho 
      uin(i,j,hop,UEINT) = eos_state % rho   *  eos_state % e
      uin(i,j,hop,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (q(i,j,hop,QU)**2 + q(i,j,hop,QV)**2 + q(i,j,hop,QW)**2))
      uin(i,j,hop,UTEMP) = eos_state % T
      do n=1, nspec
         uin(i,j,hop,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
      end do   
       
    enddo
  
  endif
  
  call destroy(eos_state)
  
  end subroutine update_ghost_cells
  
  !--------------------
  
  subroutine compute_waves(i, j, k, idir, isign, &
                           bc_type, bc_params, bc_target, &
                           T1, T2, T3, T4, T5, &
                           L1, L2, L3, L4, L5, &
                           dp, du, dv, dw, drho, &
                           q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                           qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3)
                                
  use meth_params_module, only : QVAR, QPRES, QU, QV, QW, QRHO, NQAUX, QC, QGAMC, QTEMP, QRSPEC
  use prob_params_module, only : probhi, Inflow, Outflow, SlipWall, NoSlipWall
  
  integer, intent(in) :: i, j, k, idir, isign
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
  double precision, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,qa_l3:qa_h3,NQAUX)
  double precision, intent(in) :: dp, du, dv, dw, drho
  
  integer, intent(in)          :: bc_type
  double precision, intent(in) :: bc_params(6)
  double precision, intent(in) :: bc_target(5)
  
  double precision, intent(in) :: T1, T2, T3, T4, T5
  
  ! Note: for convenience, all waves are called L because we have just 1 direction here
  ! Waves M and N will be employed for corners
  double precision, intent(out) :: L1, L2, L3, L4, L5
  
  ! Local
  double precision :: mach_local, TARGET_VX, TARGET_VY, TARGET_VZ, TARGET_TEMPERATURE, TARGET_PRESSURE
  double precision :: beta, relax_T, relax_U, relax_V, relax_W, sigma_out, Kout
  
  if ((idir == 1) .or. (idir == 2) .or. (idir == 3)) then
    continue
  else
    call bl_abort("Problem of idir in impose_NSCBC_3d:compute_waves")
  end if
  
  if ((isign == 1) .or. (isign == -1)) then
    continue
  else
    call bl_abort("Problem of isign in impose_NSCBC_3d:compute_waves")
  end if
  
  mach_local = dsqrt(q(i,j,k,QU)**2.0d0 + q(i,j,k,QV)**2.0d0 + q(i,j,k,QW)**2.0d0)/qaux(i,j,k,QC)

  !--------
  ! Recasting targets values and numerical parameters   
  TARGET_VX = bc_target(1)
  TARGET_VY = bc_target(2)
  TARGET_VZ = bc_target(3)
  TARGET_TEMPERATURE = bc_target(4)
  TARGET_PRESSURE = bc_target(5)
  
  relax_T = bc_params(1)      
  relax_U = bc_params(2)
  relax_V = bc_params(3)
  relax_W = bc_params(4)
  ! Here we have the abilities to set beta=local Mach number
  ! it may works better for outflow BCs
  if (bc_params(5) < 0.0d0) then 
    beta = mach_local
  else
    beta =  bc_params(5)
  endif
  sigma_out = bc_params(6)
  
  !--------
  ! Computing known numerical LODI waves 
  if (idir == 1) then
    ! Numerical LODI waves along X
    L1 = (q(i,j,k,QU)-qaux(i,j,k,QC))* (dp - (q(i,j,k,QRHO)*qaux(i,j,k,QC))*du)
    L2 =  q(i,j,k,QU) * ( ((qaux(i,j,k,QC)**2.0d0)*drho) - dp)
    L3 =  q(i,j,k,QU) * dv
    L4 =  q(i,j,k,QU) * dw
    L5 = (q(i,j,k,QU)+qaux(i,j,k,QC))* (dp + (q(i,j,k,QRHO)*qaux(i,j,k,QC))*du)
  elseif (idir == 2) then
    ! Numerical LODI waves along Y
    L1 = (q(i,j,k,QV)-qaux(i,j,k,QC))* (dp - (q(i,j,k,QRHO)*qaux(i,j,k,QC))*dv) 
    L2 =  q(i,j,k,QV) * du
    L3 =  q(i,j,k,QV) * ( ((qaux(i,j,k,QC)**2.0d0)*drho) - dp)
    L4 =  q(i,j,k,QV) * dw
    L5 = (q(i,j,k,QV)+qaux(i,j,k,QC))* (dp + (q(i,j,k,QRHO)*qaux(i,j,k,QC))*dv)
  elseif (idir == 3) then
    ! Numerical LODI waves along Z
    L1 = (q(i,j,k,QW)-qaux(i,j,k,QC))* (dp - (q(i,j,k,QRHO)*qaux(i,j,k,QC))*dw)
    L2 =  q(i,j,k,QW) * du
    L3 =  q(i,j,k,QW) * dv
    L4 =  q(i,j,k,QW) * ( ((qaux(i,j,k,QC)**2.0d0)*drho) - dp)
    L5 = (q(i,j,k,QW)+qaux(i,j,k,QC))* (dp + (q(i,j,k,QRHO)*qaux(i,j,k,QC))*dw)      
  endif
  
  !--------
  ! Computing missing LODI waves from BC model
  if (bc_type == Inflow) then

    if (idir == 1) then
      if (isign == 1) then      
        L5 = relax_U * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                     (q(i,j,k,QU) - TARGET_VX)  - ((1.0d0 - beta)*T5)
      elseif (isign == -1) then
        L1 = relax_U * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                        (q(i,j,k,QU) - TARGET_VX) - ((1.0d0 - beta)*T1)
      endif
    
      L2 = relax_T * (q(i,j,k,QRHO)*qaux(i,j,k,QC)*qaux(i,j,k,QRSPEC)/probhi(idir)) &
                 * (q(i,j,k,QTEMP) - TARGET_TEMPERATURE) - ((1.0d0 - beta)*T2)
      L3 = relax_V * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QV) - TARGET_VY) &
                     - ((1.0d0 - beta)*T3)
      L4 = relax_W * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QW) - TARGET_VZ) &
                     - ((1.0d0 - beta)*T4)

    elseif (idir == 2) then
    
      if (isign == 1) then      
         L5 = relax_V * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                       (q(i,j,k,QV) - TARGET_VY) - ((1.0d0 - beta)*T5)
      elseif (isign == -1) then
         L1 = relax_V * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                         (q(i,j,k,QV) - TARGET_VY)  -  ((1.0d0 - beta)*T1)
      endif
  
      L2 = relax_U * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QU) - TARGET_VX)  &
                     - ((1.0d0 - beta)*T2)
      L3 = relax_T * (q(i,j,k,QRHO)*qaux(i,j,k,QC)*qaux(i,j,k,QRSPEC)/probhi(idir)) &
                   * (q(i,j,k,QTEMP) - TARGET_TEMPERATURE) - ((1.0d0 - beta)*T3)
      L4 = relax_W * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QW) - TARGET_VZ) &
                     - ((1.0d0 - beta)*T4)

    elseif (idir == 3) then
    
      if (isign == 1) then      
         L5 = relax_W * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                       (q(i,j,k,QW) - TARGET_VZ) - ((1.0d0 - beta)*T5)
      elseif (isign == -1) then
         L1 = relax_W * ((q(i,j,k,QRHO)*qaux(i,j,k,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                         (q(i,j,k,QW) - TARGET_VZ)  -  ((1.0d0 - beta)*T1)
      endif
      
      L2 = relax_W * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QU) - TARGET_VX) &
                     - ((1.0d0 - beta)*T2)
      L3 = relax_W * (qaux(i,j,k,QC)/probhi(idir)) * (q(i,j,k,QV) - TARGET_VY) &
                     - ((1.0d0 - beta)*T3)
      L4 = relax_T * (q(i,j,k,QRHO)*qaux(i,j,k,QC)*qaux(i,j,k,QRSPEC)/probhi(idir)) &
                   * (q(i,j,k,QTEMP) - TARGET_TEMPERATURE) - ((1.0d0 - beta)*T4)
    
    else
      call bl_error("Error:: Wait, is this the fourth dimension?")
    endif
            
  elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
    ! Values long Y will be computed by mirror functions below
    ! but we set waves values to 0 to avoid undefined variables
    L1 = 0.0d0
    L2 = 0.0d0
    L3 = 0.0d0
    L4 = 0.0d0
    L5 = 0.0d0
           
  elseif (bc_type == Outflow) then
       
    Kout = sigma_out*(1.0d0 - (mach_local**2.0d0))*(qaux(i,j,k,QC)/probhi(idir))

    if (isign == 1) then
      L5 = (Kout*(q(i,j,k,QPRES) - TARGET_PRESSURE)) - ((1.0d0 - beta)*T5)
    elseif (isign == -1) then
      L1 = (Kout*(q(i,j,k,QPRES) - TARGET_PRESSURE)) - ((1.0d0 - beta)*T1)
    endif
    
  else
    call bl_error("Error:: This BC is not yet implemented for x dir in characteristic form")
  endif

  !--------
  ! Shaping the waves to be at the good dimension
  if (idir == 1) then
    L1 = L1 / (q(i,j,k,QU)-qaux(i,j,k,QC))
    L5 = L5 / (q(i,j,k,QU)+qaux(i,j,k,QC))
    if (q(i,j,k,QU) == 0.0d0) then
      L2 = 0.0d0
      L3 = 0.0d0
      L4 = 0.0d0
    else       
      L2 = L2 / q(i,j,k,QU)
      L3 = L3 / q(i,j,k,QU)
      L4 = L4 / q(i,j,k,QU)
    endif
  elseif (idir == 2) then
    L1 = L1 / (q(i,j,k,QV)-qaux(i,j,k,QC))
    L5 = L5 / (q(i,j,k,QV)+qaux(i,j,k,QC))
    if (q(i,j,k,QV) == 0.0d0) then
      L2 = 0.0d0
      L3 = 0.0d0
      L4 = 0.0d0
     else
      L2 = L2 / q(i,j,k,QV)
      L3 = L3 / q(i,j,k,QV)
      L4 = L4 / q(i,j,k,QV)
    endif
  elseif (idir == 3) then
    L1 = L1 / (q(i,j,k,QW)-qaux(i,j,k,QC))
    L5 = L5 / (q(i,j,k,QW)+qaux(i,j,k,QC))
    if (q(i,j,k,QW) == 0.0d0) then
      L2 = 0.0d0
      L3 = 0.0d0
      L4 = 0.0d0
     else
      L2 = L2 / q(i,j,k,QW)
      L3 = L3 / q(i,j,k,QW)
      L4 = L4 / q(i,j,k,QW)
    endif
  end if
  
  end subroutine compute_waves  
  
end module gc_nscbc_mod
