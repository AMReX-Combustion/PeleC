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
                         uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                         q, q_l1, q_l2, q_h1, q_h2, &
                         qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                         x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2, &
                         y_bcMask, y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2, &
                         flag_nscbc_isAnyPerio, flag_nscbc_perio, &
                         time,delta,dt,verbose) bind(C, name="impose_NSCBC")
    
  use amrex_fort_module
  use eos_module

  use amrex_constants_module
  use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, UserBC
    
  use meth_params_module, only : NVAR, NQAUX,QVAR
  use bc_fill_module, only: bcnormal
 
  implicit none
    
  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  integer, intent(in) :: uin_l1, uin_l2, uin_h1, uin_h2
  integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2
  integer, intent(in) :: y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2
  integer, intent(in) :: flag_nscbc_isAnyPerio
  integer, intent(in) :: flag_nscbc_perio(2)
  
  double precision, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  double precision, intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision, intent(inout) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2)
  integer, intent(inout) :: y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2)
  double precision, intent(in) :: delta(2), dt, time

  ! Local
  double precision dx, dy
  double precision x, y
  
  double precision :: drhodx, dudx, dvdx, dpdx
  double precision :: dpdy, dudy, dvdy, drhody
  double precision :: L1, L2, L3, L4
  double precision :: M1, M2, M3, M4
  double precision :: T1, T2, T3, T4
  double precision :: T1_X, T2_X, T3_X, T4_X
  double precision :: T1_Y, T2_Y, T3_Y, T4_Y

  double precision :: U_dummy(NVAR)
  double precision :: U_ext(NVAR)
  double precision, parameter :: small = 1.d-8
  
  integer          :: q_lo(2), q_hi(2)
  integer          :: uin_lo(2),  uin_hi(2)
  integer          :: i, j
  integer          :: bc_type, x_bc_type, y_bc_type
  integer          :: x_isign, y_isign, x_idx_Mask, y_idx_Mask
  integer          :: test_keyword_x, test_keyword_y
  double precision :: bc_params(6), x_bc_params(6), y_bc_params(6)
  double precision :: bc_target(5), x_bc_target(5), y_bc_target(5)
  
  type (eos_t) :: eos_state
  call build(eos_state)

  q_lo = [q_l1, q_l2]
  q_hi = [q_h1, q_h2]
  
  uin_lo  = [uin_l1, uin_l2]
  uin_hi  = [uin_h1, uin_h2]
 
  dx = delta(1)
  dy = delta(2)
  
  x_bcMask(:,:) = 0
  y_bcMask(:,:) = 0

  bc_params(:) = 0
  bc_target(:) = 0
  x_bc_params(:) = 0
  x_bc_target(:) = 0
  y_bc_params(:) = 0
  y_bc_target(:) = 0

  if ( flag_nscbc_isAnyPerio == 0) then

    if       (((q_hi(1) > domhi(1)) .or. (q_lo(1) < domlo(1))) &
       .and. ((q_hi(2) > domhi(2)) .or. (q_lo(2) < domlo(2)))) then
 
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
   
      x   = (dble(i)+HALF)*dx
      y   = (dble(j)+HALF)*dy
  
      ! Normal derivative along x
      call normal_derivative(i, j, 1, x_isign, dx, &
                             dpdx, dudx, dvdx, drhodx, &
                             q, q_l1, q_l2, q_h1, q_h2)
                            
      ! Normal derivative along y
      call normal_derivative(i, j, 2, y_isign, dy, &
                             dpdy, dudy, dvdy, drhody, &
                             q, q_l1, q_l2, q_h1, q_h2)
                            
      ! Compute transverse terms to x
      call compute_transverse_terms(i, j, 1,  &
                                    T1_X, T2_X, T3_X, T4_X, &
                                    dpdy, dudy, dvdy, drhody, &
                                    q, q_l1, q_l2, q_h1, q_h2, &
                                    qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                                   
      ! Compute transverse terms to y
      call compute_transverse_terms(i, j, 2,  &
                                    T1_Y, T2_Y, T3_Y, T4_Y, &
                                    dpdx, dudx, dvdx, drhodx, &
                                    q, q_l1, q_l2, q_h1, q_h2, &
                                    qaux, qa_l1, qa_l2, qa_h1, qa_h2)
   
      ! Calling user target BC values
      ! x face
      if (test_keyword_x == UserBC) then
        call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,x_isign,time,x_bc_type,x_bc_params,x_bc_target)
      else
        x_bc_type = test_keyword_x
      endif
      x_bcMask(x_idx_Mask,j) = x_bc_type
      
      ! y face
      if (test_keyword_y == UserBC) then
        call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,y_isign,time,y_bc_type,y_bc_params,y_bc_target)
      else
        y_bc_type = test_keyword_y
      endif
      y_bcMask(i,y_idx_Mask) = y_bc_type
     
      ! Computing the LODI system waves along X
      call compute_waves(i, j, 1, x_isign, &
                         x_bc_type, x_bc_params, x_bc_target, &
                         T1_X, T2_X, T3_X, T4_X, &
                         L1, L2, L3, L4, &
                         dpdx, dudx, dvdx, drhodx, &
                         q, q_l1, q_l2, q_h1, q_h2, &
                         qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                          
      ! Computing the LODI system waves along Y
      call compute_waves(i, j, 2, y_isign, &
                         y_bc_type, y_bc_params, y_bc_target, &
                         T1_Y, T2_Y, T3_Y, T4_Y, &
                         M1, M2, M3, M4, &
                         dpdy, dudy, dvdy, drhody, &
                         q, q_l1, q_l2, q_h1, q_h2, &
                         qaux, qa_l1, qa_l2, qa_h1, qa_h2)
     
      ! Along X
      ! Recomputing ghost-cells values with the LODI waves
      call update_ghost_cells(i, j, x_bc_type, 1, x_isign, dx, &
                              domlo, domhi, &
                              L1, L2, L3, L4, &
                              uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                              q, q_l1, q_l2, q_h1, q_h2, &
                              qaux, qa_l1, qa_l2, qa_h1, qa_h2)
   
      ! Along Y
      ! Recomputing ghost-cells values with the LODI waves
      call update_ghost_cells(i, j, y_bc_type, 2, y_isign, dy, &
                              domlo, domhi, &
                              M1, M2, M3, M4, &
                              uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                              q, q_l1, q_l2, q_h1, q_h2, &
                              qaux, qa_l1, qa_l2, qa_h1, qa_h2)
    
    endif
  endif ! flag_nscbc_isAnyPerio ) 

 !--------------------------------------------------------------------------   
 ! lower X
 !--------------------------------------------------------------------------

 if ((q_lo(1) < domlo(1)) .and. (physbc_lo(1) == UserBC)) then
   i = domlo(1)

   do j = q_lo(2)+1,q_hi(2)-1

      x   = (dble(i)+HALF)*dx
      y   = (dble(j)+HALF)*dy
     
      if ( flag_nscbc_isAnyPerio == 0) then
        if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners
      endif
           
     ! Normal derivative along x
     call normal_derivative(i, j, 1, 1, dx, &
                            dpdx, dudx, dvdx, drhodx, &
                            q, q_l1, q_l2, q_h1, q_h2)
 
     ! Tangential (to idir=x axis) derivative along y 
     call tangential_derivative(i, j, 1, dy, &
                               dpdy, dudy, dvdy, drhody, &
                               q, q_l1, q_l2, q_h1, q_h2)

     ! Compute transverse terms
     call compute_transverse_terms(i, j, 1,  &
                               T1, T2, T3, T4, &
                               dpdy, dudy, dvdy, drhody, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2)
     
     ! Calling user target BC values 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,1,time,bc_type,bc_params,bc_target)
     
     ! Filling bcMask with specific user defined BC type
     if ((j < q_lo(2)+3) .or. (j > q_hi(2)-3)) then
       continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
     else
       x_bcMask(i,j) = bc_type
     endif
          
     ! Computing the LODI system waves
     call compute_waves(i, j, 1, 1, &
                        bc_type, bc_params, bc_target, &
                        T1, T2, T3, T4, &
                        L1, L2, L3, L4, &
                        dpdx, dudx, dvdx, drhodx, &
                        q, q_l1, q_l2, q_h1, q_h2, &
                        qaux, qa_l1, qa_l2, qa_h1, qa_h2)
          
     ! Recomputing ghost-cells values with the LODI waves
     call update_ghost_cells(i, j, bc_type, 1, 1, dx, &
                             domlo, domhi, &
                             L1, L2, L3, L4, &
                             uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                             q, q_l1, q_l2, q_h1, q_h2, &
                             qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                               
   enddo
 end if
  
 !--------------------------------------------------------------------------   
 ! upper X
 !--------------------------------------------------------------------------
 
 if ((q_hi(1) > domhi(1)) .and. (physbc_hi(1) == UserBC)) then
   i = domhi(1)

   do j = q_lo(2)+1,q_hi(2)-1
   
     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners 
     endif
     
     ! Normal derivative along x
     call normal_derivative(i, j, 1, -1, dx, &
                            dpdx, dudx, dvdx, drhodx, &
                            q, q_l1, q_l2, q_h1, q_h2)
       
     ! Tangential (to idir=x axis) derivative along y 
     call tangential_derivative(i, j, 1, dy, &
                               dpdy, dudy, dvdy, drhody, &
                               q, q_l1, q_l2, q_h1, q_h2)
  
     ! Compute transverse terms
     call compute_transverse_terms(i, j, 1,  &
                               T1, T2, T3, T4, &
                               dpdy, dudy, dvdy, drhody, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2)
     
     ! Filling bcMask with specific user defined BC type 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,-1,time,bc_type,bc_params,bc_target)
     if ((j < q_lo(2)+3) .or. (j > q_hi(2)-3)) then
       continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
     else
       x_bcMask(i+1,j) = bc_type
     endif
 
     ! Computing the LODI system waves
     call compute_waves(i, j, 1, -1, &
                        bc_type, bc_params, bc_target, &
                        T1, T2, T3, T4, &
                        L1, L2, L3, L4, &
                        dpdx, dudx, dvdx, drhodx, &
                        q, q_l1, q_l2, q_h1, q_h2, &
                        qaux, qa_l1, qa_l2, qa_h1, qa_h2)    
     
     ! Recomputing ghost-cells values with the LODI waves
     call update_ghost_cells(i, j, bc_type, 1, -1, dx, &
                             domlo, domhi, &
                             L1, L2, L3, L4, &
                             uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                             q, q_l1, q_l2, q_h1, q_h2, &
                             qaux, qa_l1, qa_l2, qa_h1, qa_h2)

  enddo
endif

 !--------------------------------------------------------------------------   
 ! lower Y
 !--------------------------------------------------------------------------

 if ((q_lo(2) < domlo(2)) .and. (physbc_lo(2) == UserBC)) then
 
   j = domlo(2)
      
   do i = q_lo(1)+1,q_hi(1)-1

     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
          
     ! Normal derivative along y
     call normal_derivative(i, j, 2, 1, dy, &
                            dpdy, dudy, dvdy, drhody, &
                            q, q_l1, q_l2, q_h1, q_h2)

     ! Tangential (to idir=y axis) derivative along x 
     call tangential_derivative(i, j, 2, dx, &
                               dpdx, dudx, dvdx, drhodx, &
                               q, q_l1, q_l2, q_h1, q_h2)
      
     ! Compute transverse terms
     call compute_transverse_terms(i, j, 2,  &
                               T1, T2, T3, T4, &
                               dpdx, dudx, dvdx, drhodx, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                               
     ! Filling bcMask with specific user defined BC type
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,1,time,bc_type,bc_params,bc_target)
     if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3)) then
       continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
     else
       y_bcMask(i,j) = bc_type
     endif
     
     ! Computing the LODI system waves
     call compute_waves(i, j, 2, 1, &
                        bc_type, bc_params, bc_target, &
                        T1, T2, T3, T4, &
                        L1, L2, L3, L4, &
                        dpdy, dudy, dvdy, drhody, &
                        q, q_l1, q_l2, q_h1, q_h2, &
                        qaux, qa_l1, qa_l2, qa_h1, qa_h2)
 
     ! Recomputing ghost-cells values with the LODI waves
     call update_ghost_cells(i, j, bc_type, 2, 1, dy, &
                             domlo, domhi, &
                             L1, L2, L3, L4, &
                             uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                             q, q_l1, q_l2, q_h1, q_h2, &
                             qaux, qa_l1, qa_l2, qa_h1, qa_h2)
   
   enddo
 end if
     
!--------------------------------------------------------------------------   
! upper Y
!--------------------------------------------------------------------------

 if ((q_hi(2) > domhi(2)) .and. (physbc_hi(2) == UserBC)) then
 
   j = domhi(2)
      
   do i = q_lo(1)+1,q_hi(1)-1
     
     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
        
     ! Normal derivative along y
     call normal_derivative(i, j, 2, -1, dy, &
                            dpdy, dudy, dvdy, drhody, &
                            q, q_l1, q_l2, q_h1, q_h2)

     ! Tangential (to idir=y axis) derivative along x 
     call tangential_derivative(i, j, 2, dx, &
                               dpdx, dudx, dvdx, drhodx, &
                               q, q_l1, q_l2, q_h1, q_h2)
      
     ! Compute transverse terms
     call compute_transverse_terms(i, j, 2,  &
                               T1, T2, T3, T4, &
                               dpdx, dudx, dvdx, drhodx, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                               
     ! Filling bcMask with specific user defined BC type 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,-1,time,bc_type,bc_params,bc_target)
     if ((i < q_lo(1)+3) .or. (i > q_hi(1)-3)) then
       continue ! There is just 1 ghost-cell with bcMask because of the Riemann solver
     else
       y_bcMask(i,j+1) = bc_type
     endif

     ! Computing the LODI system waves
     call compute_waves(i, j, 2, -1, &
                        bc_type, bc_params, bc_target, &
                        T1, T2, T3, T4, &
                        L1, L2, L3, L4, &
                        dpdy, dudy, dvdy, drhody, &
                        q, q_l1, q_l2, q_h1, q_h2, &
                        qaux, qa_l1, qa_l2, qa_h1, qa_h2)

     ! Recomputing ghost-cells values with the LODI waves
     call update_ghost_cells(i, j, bc_type, 2, -1, dy, &
                             domlo, domhi, &
                             L1, L2, L3, L4, &
                             uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                             q, q_l1, q_l2, q_h1, q_h2, &
                             qaux, qa_l1, qa_l2, qa_h1, qa_h2)
   
   enddo
end if

call destroy(eos_state)

end subroutine impose_NSCBC


!-------------------------------------------------
! Generic routines below
!-------------------------------------------------


  subroutine normal_derivative(i, j, idir, isign, delta, &
                               dp, du, dv, drho, &
                               q, q_l1, q_l2, q_h1, q_h2)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QRHO

  implicit none
  
  integer, intent(in) :: i,j,idir,isign
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  double precision, intent(in) :: delta
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  
  double precision, intent(out) :: dp, du, dv, drho
  
  
  if (idir == 1) then
    if (isign == 1) then
    
      ! 2nd order
      dp = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i+1,j,QPRES)-0.5d0*q(i+2,j,QPRES))/delta
      du = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i+1,j,QU)-0.5d0*q(i+2,j,QU))/delta
      dv = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i+1,j,QV)-0.5d0*q(i+2,j,QV))/delta
      drho = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i+1,j,QRHO)-0.5d0*q(i+2,j,QRHO))/delta
      
    elseif (isign == -1) then
    
       !2nd order
       dp = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i-1,j,QPRES)+0.5d0*q(i-2,j,QPRES))/delta
       du = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i-1,j,QU)+0.5d0*q(i-2,j,QU))/delta
       dv = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i-1,j,QV)+0.5d0*q(i-2,j,QV))/delta
       drho = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i-1,j,QRHO)+0.5d0*q(i-2,j,QRHO))/delta
    
    else
      call bl_abort("Problem of isign in impose_NSCBC_2d:normal_derivative")
    end if
    
  elseif (idir == 2) then
  
    if (isign == 1) then
  
      !2nd order
      dp = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i,j+1,QPRES)-0.5d0*q(i,j+2,QPRES))/delta
      du = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i,j+1,QU)-0.5d0*q(i,j+2,QU))/delta
      dv = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i,j+1,QV)-0.5d0*q(i,j+2,QV))/delta
      drho = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i,j+1,QRHO)-0.5d0*q(i,j+2,QRHO))/delta
     
    elseif (isign == -1) then
    
      !2nd order
      dp = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i,j-1,QPRES)+0.5d0*q(i,j-2,QPRES))/delta
      du = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i,j-1,QU)+0.5d0*q(i,j-2,QU))/delta
      dv = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i,j-1,QV)+0.5d0*q(i,j-2,QV))/delta
      drho = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i,j-1,QRHO)+0.5d0*q(i,j-2,QRHO))/delta
  
    else
      call bl_abort("Problem of isign in impose_NSCBC_2d:normal_derivative")
    end if
  
  else
      call bl_abort("Problem of idir in impose_NSCBC_2d:normal_derivative")
  end if
  
  end subroutine normal_derivative
  
  !----------------
  
  subroutine tangential_derivative(i, j, idir, delta, &
                               dp, du, dv, drho, &
                               q, q_l1, q_l2, q_h1, q_h2)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QRHO

  implicit none
  
  integer, intent(in) :: i,j,idir
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  double precision, intent(in) :: delta
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  
  double precision, intent(out) :: dp, du, dv, drho
  
  ! Warning, idir means the normal direction and this routine compute the tangential derivative
  if (idir == 1) then

    ! 2nd order Central
    dp = (q(i,j+1,QPRES)-q(i,j-1,QPRES))/(2.0d0*delta)
    du = (q(i,j+1,QU)-q(i,j-1,QU))/(2.0d0*delta)
    dv = (q(i,j+1,QV)-q(i,j-1,QV))/(2.0d0*delta)
    drho = (q(i,j+1,QRHO)-q(i,j-1,QRHO))/(2.0d0*delta)
    
  elseif (idir == 2) then
  
    ! 2nd order Central
    dp = (q(i+1,j,QPRES)-q(i-1,j,QPRES))/(2.0d0*delta)
    du = (q(i+1,j,QU)-q(i-1,j,QU))/(2.0d0*delta)
    dv = (q(i+1,j,QV)-q(i-1,j,QV))/(2.0d0*delta)
    drho = (q(i+1,j,QRHO)-q(i-1,j,QRHO))/(2.0d0*delta)
    
  else
      call bl_abort("Problem of idir in impose_NSCBC_2d:tangential_derivative")
  end if
  
  end subroutine tangential_derivative
  
  !-----------------
  
  subroutine compute_transverse_terms(i, j, idir,  &
                               T1, T2, T3, T4, &
                               dp, du, dv, drho, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                               
                               
  use meth_params_module, only : QVAR, QPRES, QU, QV, QRHO, NQAUX, QC, QGAMC

  implicit none
  
  integer, intent(in) :: i,j,idir
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  double precision, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision, intent(in) :: dp, du, dv, drho
  
  double precision, intent(out) :: T1, T2, T3, T4
  
  if (idir == 1) then
  
     T1 = (q(i,j,QV)*(dp - q(i,j,QRHO)*qaux(i,j,QC)*du)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dv)
     T2 = (q(i,j,QV)*((qaux(i,j,QC)*qaux(i,j,QC)*drho)-dp)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dv) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dv)
     T3 = ((q(i,j,QV)*dv))+(dp/q(i,j,QRHO))
     T4 = (q(i,j,QV)*(dp + q(i,j,QRHO)*qaux(i,j,QC)*du)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dv)
  
  
  elseif (idir == 2) then
  
     T1 = (q(i,j,QU)*(dp - q(i,j,QRHO)*qaux(i,j,QC)*dv)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*du)
     T2 = ((q(i,j,QU)*du))+(dp/q(i,j,QRHO))
     T3 = (q(i,j,QU)*((qaux(i,j,QC)*qaux(i,j,QC)*drho)-dp)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*du) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*du)
     T4 = (q(i,j,QU)*(dp + q(i,j,QRHO)*qaux(i,j,QC)*dv)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*du)
           
  else
      call bl_abort("Problem of idir in impose_NSCBC_2d:compute_transverse_terms")
  end if

  end subroutine compute_transverse_terms
    
  !-------------------------
  
  subroutine update_ghost_cells(i, j, bc_type, idir, isign, delta, &
                                domlo, domhi, &
                                L1, L2, L3, L4, &
                                uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                                q, q_l1, q_l2, q_h1, q_h2, &
                                qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                               
  use eos_module
  use fuego_chemistry, only : nspecies
  use amrex_constants_module, only : ONE
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP,&
                                 UFS, NQAUX, QC, QGAMC, QRSPEC, &
                                 QC, QDPDE, QDPDR, QCSML, QGAMC, &
                                 QVAR, QRHO, QU, QV, QREINT, QPRES, QTEMP, &
                                 QFS, QFX, QGAME, NHYP
  use prob_params_module, only : SlipWall, NoSlipWall

  implicit none
  
  integer, intent(in) :: i,j,idir,isign,bc_type
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  integer, intent(in) :: uin_l1, uin_l2, uin_h1, uin_h2
  double precision, intent(in) :: L1, L2, L3, L4
  double precision, intent(in) :: delta
  double precision, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  double precision, intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision, intent(inout) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    
  integer :: idx_gc1, idx_gc2, idx_gc3, idx_gc4, idx_int1, idx_int2, idx_int3
  integer :: idx_start, idx_end, hop, n, local_index, bndy_loc
  double precision :: drho, du, dv, dp, wall_sign
  double precision, parameter :: small = 1.d-8
  
  type (eos_t) :: eos_state
  
  call build(eos_state)
  
  if ((idir == 1) .or. (idir == 2)) then
    continue
  else
    call bl_abort("Problem of idir in impose_NSCBC_2d:compute_waves")
  end if
  
  if ((isign == 1) .or. (isign == -1)) then
    continue
  else
    call bl_abort("Problem of isign in impose_NSCBC_2d:compute_waves")
  end if
  
  ! Compute new spatial derivative
  
  if (idir == 1) then
    local_index = i
    drho = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)  
    du   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
    dv   = L3
    dp   = 0.5d0*(L1+L4) 
  elseif (idir == 2) then
    local_index = j
    drho = (L3 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)
    du   = L2
    dv   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
    dp   = 0.5d0*(L1+L4)
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
      bndy_loc = domlo(idir)
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
      bndy_loc = domhi(idir)
    endif
    
    if (idir == 1) then
    
     ! Update ghost cells
     ! 2nd order
     q(idx_gc1,j,QU)    = q(idx_int1,j,QU) - 2.0d0*delta*du*isign
     q(idx_gc1,j,QV)    = q(idx_int1,j,QV) - 2.0d0*delta*dv*isign 
     q(idx_gc1,j,QRHO)  = q(idx_int1,j,QRHO)  - 2.0d0*delta*drho*isign
     q(idx_gc1,j,QPRES) = q(idx_int1,j,QPRES) - 2.0d0*delta*dp*isign
   
     !---------------- 
     q(idx_gc2,j,QU)    = -2.0d0*q(idx_int1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(idx_gc1,j,QU) + 6.0d0*delta*du*isign
     q(idx_gc2,j,QV)    = -2.0d0*q(idx_int1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(idx_gc1,j,QV) + 6.0d0*delta*dv*isign
     q(idx_gc2,j,QRHO)    = -2.0d0*q(idx_int1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(idx_gc1,j,QRHO) + 6.0d0*delta*drho*isign
     q(idx_gc2,j,QPRES)    = -2.0d0*q(idx_int1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(idx_gc1,j,QPRES) + 6.0d0*delta*dp*isign
  
     q(idx_gc3,j,QU)    = 3.0d0*q(idx_int1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(idx_gc1,j,QU) + 6.0d0*q(idx_gc2,j,QU) - 12.0d0*delta*du*isign
     q(idx_gc3,j,QV)    = 3.0d0*q(idx_int1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(idx_gc1,j,QV) + 6.0d0*q(idx_gc2,j,QV) - 12.0d0*delta*dv*isign
     q(idx_gc3,j,QRHO)    = 3.0d0*q(idx_int1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(idx_gc1,j,QRHO) + 6.0d0*q(idx_gc2,j,QRHO) - 12.0d0*delta*drho*isign
     q(idx_gc3,j,QPRES)    = 3.0d0*q(idx_int1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(idx_gc1,j,QPRES) + 6.0d0*q(idx_gc2,j,QPRES) - 12.0d0*delta*dp*isign
  
     q(idx_gc4,j,QU)    = -2.0d0*q(idx_int1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(idx_gc1,j,QU) - 12.0d0*q(idx_gc2,j,QU)  &
                     + 4.0d0*q(idx_gc3,j,QU) + 12.0d0*delta*du*isign
     q(idx_gc4,j,QV)    = -2.0d0*q(idx_int1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(idx_gc1,j,QV) - 12.0d0*q(idx_gc2,j,QV) &
                     + 4.0d0*q(idx_gc3,j,QV) + 12.0d0*delta*dv*isign
     q(idx_gc4,j,QRHO)  = -2.0d0*q(idx_int1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(idx_gc1,j,QRHO) - 12.0d0*q(idx_gc2,j,QRHO) &
                     + 4.0d0*q(idx_gc3,j,QRHO) + 12.0d0*delta*drho*isign
     q(idx_gc4,j,QPRES) = -2.0d0*q(idx_int1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(idx_gc1,j,QPRES) - 12.0d0*q(idx_gc2,j,QPRES) &
                     + 4.0d0*q(idx_gc3,j,QPRES) + 12.0d0*delta*dp*isign
  
     if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(idx_gc1,j,QU)    = -q(i,j,QU)
       q(idx_gc2,j,QU)    = -q(idx_int1,j,QU)
       q(idx_gc3,j,QU)    = -q(idx_int2,j,QU)
       q(idx_gc4,j,QU)    = -q(idx_int3,j,QU)
  
       q(idx_gc1,j,QV)    = wall_sign*q(i,j,QV)
       q(idx_gc2,j,QV)    = wall_sign*q(idx_int1,j,QV)
       q(idx_gc3,j,QV)    = wall_sign*q(idx_int2,j,QV)
       q(idx_gc4,j,QV)    = wall_sign*q(idx_int3,j,QV)
       
       q(idx_gc1,j,QRHO)  = q(i,j,QRHO)
       q(idx_gc2,j,QRHO)  = q(idx_int1,j,QRHO)
       q(idx_gc3,j,QRHO)  = q(idx_int2,j,QRHO)
       q(idx_gc4,j,QRHO)  = q(idx_int3,j,QRHO)
       
       q(idx_gc1,j,QPRES)  = q(i,j,QPRES)
       q(idx_gc2,j,QPRES)  = q(idx_int1,j,QPRES)
       q(idx_gc3,j,QPRES)  = q(idx_int2,j,QPRES)
       q(idx_gc4,j,QPRES)  = q(idx_int3,j,QPRES)
     
     end if
  
     ! Recompute missing values thanks to EOS
     do hop=idx_start,idx_end,-isign
     
       eos_state % p        = q(hop,j,QPRES )
       eos_state % rho      = q(hop,j,QRHO  )
       eos_state % massfrac = q(bndy_loc,j,QFS:QFS+nspecies-1)
       eos_state % aux      = q(bndy_loc,j,QFX:QFX+naux-1)
  
       call eos_rp(eos_state)
       q(hop,j,QTEMP)  = eos_state % T
       q(hop,j,QREINT) = eos_state % e * q(hop,j,QRHO)
       q(hop,j,QGAME)  = q(hop,j,QPRES) / q(hop,j,QREINT) + ONE
       
       qaux(hop,j,QDPDR)  = eos_state % dpdr_e
       qaux(hop,j,QDPDE)  = eos_state % dpde
       qaux(hop,j,QGAMC)  = eos_state % gam1
       qaux(hop,j,QC   )  = eos_state % cs
       qaux(hop,j,QCSML)  = max(small, small * qaux(hop,j,QC))
  
       ! Here the update of the conservative variables uin seems to have only an impact
       ! on the application of artificial viscosity difmag. 
       uin(hop,j,URHO )  = eos_state % rho 
       uin(hop,j,UMX  )  = q(hop,j,QU ) * eos_state % rho 
       uin(hop,j,UMY  )  = q(hop,j,QV ) * eos_state % rho 
       uin(hop,j,UMZ  ) = 0.0d0
       uin(hop,j,UEINT) = eos_state % rho   *  eos_state % e
       uin(hop,j,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (q(hop,j,QU)**2 + q(hop,j,QV)**2))
       uin(hop,j,UTEMP) = eos_state % T
       do n=1, nspecies
          uin(hop,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do   
       
     enddo
  
  elseif (idir == 2) then
       
     ! Update ghost cells
     ! 2nd order
     q(i,idx_gc1,QU)    = q(i,idx_int1,QU) - 2.0d0*delta*du*isign
     q(i,idx_gc1,QV)    = q(i,idx_int1,QV) - 2.0d0*delta*dv*isign 
     q(i,idx_gc1,QRHO)  = q(i,idx_int1,QRHO)  - 2.0d0*delta*drho*isign
     q(i,idx_gc1,QPRES) = q(i,idx_int1,QPRES) - 2.0d0*delta*dp*isign
   
     !---------------- 
     q(i,idx_gc2,QU)    = -2.0d0*q(i,idx_int1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,idx_gc1,QU) + 6.0d0*delta*du*isign
     q(i,idx_gc2,QV)    = -2.0d0*q(i,idx_int1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,idx_gc1,QV) + 6.0d0*delta*dv*isign
     q(i,idx_gc2,QRHO)  = -2.0d0*q(i,idx_int1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,idx_gc1,QRHO) + 6.0d0*delta*drho*isign
     q(i,idx_gc2,QPRES) = -2.0d0*q(i,idx_int1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,idx_gc1,QPRES) + 6.0d0*delta*dp*isign
  
     q(i,idx_gc3,QU)    = 3.0d0*q(i,idx_int1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,idx_gc1,QU) + 6.0d0*q(i,idx_gc2,QU) - 12.0d0*delta*du*isign
     q(i,idx_gc3,QV)    = 3.0d0*q(i,idx_int1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,idx_gc1,QV) + 6.0d0*q(i,idx_gc2,QV) - 12.0d0*delta*dv*isign
     q(i,idx_gc3,QRHO)  = 3.0d0*q(i,idx_int1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,idx_gc1,QRHO) + 6.0d0*q(i,idx_gc2,QRHO) - 12.0d0*delta*drho*isign
     q(i,idx_gc3,QPRES) = 3.0d0*q(i,idx_int1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,idx_gc1,QPRES) + 6.0d0*q(i,idx_gc2,QPRES) - 12.0d0*delta*dp*isign
  
     q(i,idx_gc4,QU)    = -2.0d0*q(i,idx_int1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,idx_gc1,QU) - 12.0d0*q(i,idx_gc2,QU)  &
                     + 4.0d0*q(i,idx_gc3,QU) + 12.0d0*delta*du*isign
     q(i,idx_gc4,QV)    = -2.0d0*q(i,idx_int1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,idx_gc1,QV) - 12.0d0*q(i,idx_gc2,QV) &
                     + 4.0d0*q(i,idx_gc3,QV) + 12.0d0*delta*dv*isign
     q(i,idx_gc4,QRHO)  = -2.0d0*q(i,idx_int1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,idx_gc1,QRHO) - 12.0d0*q(i,idx_gc2,QRHO) &
                     + 4.0d0*q(i,idx_gc3,QRHO) + 12.0d0*delta*drho*isign
     q(i,idx_gc4,QPRES) = -2.0d0*q(i,idx_int1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,idx_gc1,QPRES) - 12.0d0*q(i,idx_gc2,QPRES) &
                     + 4.0d0*q(i,idx_gc3,QPRES) + 12.0d0*delta*dp*isign
  
     if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i,idx_gc1,QU)    = wall_sign*q(i,j,QU)
       q(i,idx_gc2,QU)    = wall_sign*q(i,idx_int1,QU)
       q(i,idx_gc3,QU)    = wall_sign*q(i,idx_int2,QU)
       q(i,idx_gc4,QU)    = wall_sign*q(i,idx_int3,QU)
       
       q(i,idx_gc1,QV)    = -q(i,j,QV)
       q(i,idx_gc2,QV)    = -q(i,idx_int1,QV)
       q(i,idx_gc3,QV)    = -q(i,idx_int2,QV)
       q(i,idx_gc4,QV)    = -q(i,idx_int3,QV)
       
       q(i,idx_gc1,QRHO)  = q(i,j,QRHO)
       q(i,idx_gc2,QRHO)  = q(i,idx_int1,QRHO)
       q(i,idx_gc3,QRHO)  = q(i,idx_int2,QRHO)
       q(i,idx_gc4,QRHO)  = q(i,idx_int3,QRHO)
       
       q(i,idx_gc1,QPRES)  = q(i,j,QPRES)
       q(i,idx_gc2,QPRES)  = q(i,idx_int1,QPRES)
       q(i,idx_gc3,QPRES)  = q(i,idx_int2,QPRES)
       q(i,idx_gc4,QPRES)  = q(i,idx_int3,QPRES)
     
     end if
      
     ! Recompute missing values thanks to EOS
     do hop=idx_start,idx_end,-isign
     
       eos_state % p        = q(i,hop,QPRES )
       eos_state % rho      = q(i,hop,QRHO  )
       eos_state % massfrac = q(i,bndy_loc,QFS:QFS+nspecies-1)
       eos_state % aux      = q(i,bndy_loc,QFX:QFX+naux-1)
  
       call eos_rp(eos_state)
       q(i,hop,QTEMP)  = eos_state % T
       q(i,hop,QREINT) = eos_state % e * q(i,hop,QRHO)
       q(i,hop,QGAME)  = q(i,hop,QPRES) / q(i,hop,QREINT) + ONE
       
       qaux(i,hop,QDPDR)  = eos_state % dpdr_e
       qaux(i,hop,QDPDE)  = eos_state % dpde
       qaux(i,hop,QGAMC)  = eos_state % gam1
       qaux(i,hop,QC   )  = eos_state % cs
       qaux(i,hop,QCSML)  = max(small, small * qaux(i,hop,QC))
  
       ! Here the update of the conservative variables uin seems to have only an impact
       ! on the application of artificial viscosity difmag. 
       uin(i,hop,URHO )  = eos_state % rho 
       uin(i,hop,UMX  )  = q(i,hop,QU ) * eos_state % rho 
       uin(i,hop,UMY  )  = q(i,hop,QV ) * eos_state % rho 
       uin(i,hop,UMZ  ) = 0.0d0
       uin(i,hop,UEINT) = eos_state % rho   *  eos_state % e
       uin(i,hop,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (q(i,hop,QU)**2 + q(i,hop,QV)**2))
       uin(i,hop,UTEMP) = eos_state % T
       do n=1, nspecies
          uin(i,hop,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do   
       
     enddo
    
  endif
  
  call destroy(eos_state)
  
  end subroutine update_ghost_cells
  
  !--------------------
  
  subroutine compute_waves(i, j, idir, isign, &
                           bc_type, bc_params, bc_target, &
                           T1, T2, T3, T4, &
                           L1, L2, L3, L4, &
                           dp, du, dv, drho, &
                           q, q_l1, q_l2, q_h1, q_h2, &
                           qaux, qa_l1, qa_l2, qa_h1, qa_h2)
                                                    
  use meth_params_module, only : QVAR, QPRES, QU, QV, QRHO, NQAUX, QC, QGAMC, QTEMP, QRSPEC
  use prob_params_module, only : probhi, UserBC, Inflow, Outflow, SlipWall, NoSlipWall

  implicit none

  integer, intent(in) :: i,j,idir,isign
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
  double precision, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision, intent(in) :: dp, du, dv, drho
  
  integer, intent(in)          :: bc_type
  double precision, intent(in) :: bc_params(6)
  double precision, intent(in) :: bc_target(5)
  
  double precision, intent(in) :: T1, T2, T3, T4
  
  double precision, intent(out) :: L1, L2, L3, L4
  
  ! Local
  double precision :: mach_local, TARGET_VX, TARGET_VY, TARGET_TEMPERATURE, TARGET_PRESSURE
  double precision :: beta, relax_T, relax_U, relax_V, sigma_out, Kout
  
  if ((idir == 1) .or. (idir == 2)) then
    continue
  else
    call bl_abort("Problem of idir in impose_NSCBC_2d:compute_waves")
  end if
  
  if ((isign == 1) .or. (isign == -1)) then
    continue
  else
    call bl_abort("Problem of isign in impose_NSCBC_2d:compute_waves")
  end if
  
  mach_local = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)
     
  !--------
  ! Recasting targets values and numerical parameters   
  TARGET_VX = bc_target(1)
  TARGET_VY = bc_target(2)
  TARGET_TEMPERATURE = bc_target(4)
  TARGET_PRESSURE = bc_target(5)
  
  relax_T = bc_params(1)      
  relax_U = bc_params(2)
  relax_V = bc_params(3)
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
    L1 = (q(i,j,QU)-qaux(i,j,QC))* (dp - (q(i,j,QRHO)*qaux(i,j,QC))*du)
    L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drho) - dp)
    L3 = q(i,j,QU) * dv
    L4 = (q(i,j,QU)+qaux(i,j,QC))* (dp + (q(i,j,QRHO)*qaux(i,j,QC))*du)
  elseif (idir == 2) then
    ! Numerical LODI waves along Y
    L1 = (q(i,j,QV)-qaux(i,j,QC))* (dp - (q(i,j,QRHO)*qaux(i,j,QC))*dv) 
    L2 = q(i,j,QV) * du
    L3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drho) - dp)
    L4 = (q(i,j,QV)+qaux(i,j,QC))* (dp + (q(i,j,QRHO)*qaux(i,j,QC))*dv)
  endif
     
  !--------
  ! Computing missing LODI waves from BC model
  if (bc_type == Inflow) then
        
    if (idir == 1) then
      if (isign == 1) then      
        L4 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                     (q(i,j,QU) - TARGET_VX)  - ((1.0d0 - beta)*T4)
      elseif (isign == -1) then
        L1 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                        (q(i,j,QU) - TARGET_VX) - ((1.0d0 - beta)*T1)
      endif
    
      L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(idir)) &
                 * (q(i,j,QTEMP) - TARGET_TEMPERATURE) - ((1.0d0 - beta)*T2)
      L3 = relax_V * (qaux(i,j,QC)/probhi(idir)) * (q(i,j,QV) - TARGET_VY) &
                     - ((1.0d0 - beta)*T3)

    elseif (idir == 2) then
    
      if (isign == 1) then      
         L4 = relax_V * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                       (q(i,j,QV) - TARGET_VY) - ((1.0d0 - beta)*T4)
      elseif (isign == -1) then
         L1 = relax_V * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local*mach_local)/probhi(idir)) * &
                         (q(i,j,QV) - TARGET_VY)  -  ((1.0d0 - beta)*T1)
      endif
  
      L2 = relax_U * (qaux(i,j,QC)/probhi(idir)) * (q(i,j,QU) - TARGET_VX)  &
                     - ((1.0d0 - beta)*T2)
      L3 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(idir)) &
                   * (q(i,j,QTEMP) - TARGET_TEMPERATURE) - ((1.0d0 - beta)*T3)

    endif
            
  elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
    ! Values long Y will be computed by mirror functions below
    ! but we set waves values to 0 to avoid undefined variables
    L1 = 0.0d0
    L2 = 0.0d0
    L3 = 0.0d0
    L4 = 0.0d0
       
  elseif (bc_type == Outflow) then
       
    Kout = sigma_out*(1.0d0 - (mach_local**2.0d0))*(qaux(i,j,QC)/probhi(idir))
    
    if (isign == 1) then
      L4 = (Kout*(q(i,j,QPRES) - TARGET_PRESSURE)) - ((1.0d0 - beta)*T4)
    elseif (isign == -1) then
      L1 = (Kout*(q(i,j,QPRES) - TARGET_PRESSURE)) - ((1.0d0 - beta)*T1)
    endif

  else
    call bl_error("Error:: This BC is not yet implemented in characteristic form")
  endif
 
  !--------
  ! Shaping the waves to be at the good dimension
  if (idir == 1) then
    L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
    L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
    if (q(i,j,QU) == 0.0d0) then  
      L2 = 0.0d0
      L3 = 0.0d0
    else       
      L2 = L2 / q(i,j,QU)
      L3 = L3 / q(i,j,QU)
    endif
  elseif (idir == 2) then
    L1 = L1 / (q(i,j,QV)-qaux(i,j,QC))
    L4 = L4 / (q(i,j,QV)+qaux(i,j,QC))
    if (q(i,j,QV) == 0.0d0) then
      L2 = 0.0d0
      L3 = 0.0d0
    else
      L2 = L2 / q(i,j,QV)
      L3 = L3 / q(i,j,QV)
    endif
  end if
  
  end subroutine compute_waves
  
end module gc_nscbc_mod
