module gc_nscbc_mod

  private
  public impose_NSCBC

contains

 subroutine impose_NSCBC(lo, hi, domlo, domhi, &
                                    uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                                   q, q_l1, q_l2, q_h1, q_h2, &
                                   qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                                   x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2, &
                                   y_bcMask, y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2, &
                                   flag_nscbc_isAnyPerio, flag_nscbc_perio, &
                                   time,delta,dt,verbose) bind(C, name="impose_NSCBC")
    
 
    use bl_error_module
    use network, only : nspec
    use eos_module
    use fundamental_constants_module, only: k_B, n_A

    use bl_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, &
                                   Interior, Inflow, Outflow, SlipWall, NoSlipWall
    
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP,&
                                   UFS, NQAUX, QC, QGAMC, QRSPEC, &
                                   QC, QDPDE, QDPDR, QCSML, QGAMC, &
                                   QVAR, QRHO, QU, QV, QREINT, QPRES, QTEMP, &
                                   QFS, QFX, QGAME, NHYP
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
  double precision :: L1, L2, L3, L4, T1, T2, T3, T4
  double precision :: S1, S2, M1, M2, M3, M4
  double precision :: Kout, sigma_out, beta
  double precision :: relax_U, relax_V, relax_T
  double precision :: mach_local_hi_x, mach_local_lo_x
  double precision :: mach_local_hi_y, mach_local_lo_y
  double precision :: mach_local_corner
  double precision :: INLET_PRESSURE, INLET_VX, INLET_VY, INLET_TEMPERATURE
  double precision :: INLET_PRESSURE_Xdir, INLET_VX_Xdir, INLET_VY_Xdir, INLET_TEMPERATURE_Xdir
  double precision :: INLET_PRESSURE_Ydir, INLET_VX_Ydir, INLET_VY_Ydir, INLET_TEMPERATURE_Ydir
  
  double precision :: U_dummy(NVAR)
  double precision :: U_ext(NVAR)
  double precision, parameter :: small = 1.d-8
  
  integer          :: q_lo(2), q_hi(2)
  integer          ::  uin_lo(2),  uin_hi(2)
  integer          :: i, j, n, hop
  integer          :: bc_type
  double precision :: bc_params(6)
  double precision :: wall_sign, beta_X, beta_Y
  
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
write(*,*) 'IN ROUTINE NSCBC, flag_nscbc_isAnyPerio ',flag_nscbc_isAnyPerio
write(*,*) 'IN ROUTINE NSCBC, flag_nscbc_perio',flag_nscbc_perio
  if ( flag_nscbc_isAnyPerio == 0) then
write(*,*) 'WE DONT HAVE PERIO, SO WE HAVE CORNERS'


 !--------------------------------------------------------------------------   
 ! upper right corner (Lx,Ly)
 ! phi = 1, psi = 1
 !--------------------------------------------------------------------------

 if ((q_hi(1) > domhi(1)) .and. (q_hi(2) > domhi(2))) then

   i = domhi(1)
   j = domhi(2)

   x   = (dble(i)+HALF)*dx
   y   = (dble(j)+HALF)*dy
   
   ! Calling user target BC values
   ! right face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,-1,.false.,bc_type,bc_params)
   x_bcMask(i+1,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i+1,j+1,1) = bc_type
   
   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta = bc_params(5)
    
   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Xdir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Xdir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Xdir = U_ext(UTEMP)
   INLET_PRESSURE_Xdir = eos_state%p
 
   ! top face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,-1,.false.,bc_type,bc_params)
   y_bcMask(i,j+1) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i+1,j+1,1) = bc_type
   
   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta = bc_params(5)
   
   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Ydir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Ydir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Ydir = U_ext(UTEMP)
   INLET_PRESSURE_Ydir = eos_state%p
 
   mach_local_corner = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)

   ! Derivatives along Y
   dpdy = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i,j-1,QPRES)+0.5d0*q(i,j-2,QPRES))/dy
   dudy = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i,j-1,QU)+0.5d0*q(i,j-2,QU))/dy
   dvdy = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i,j-1,QV)+0.5d0*q(i,j-2,QV))/dy
   drhody = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i,j-1,QRHO)+0.5d0*q(i,j-2,QRHO))/dy

   ! Derivative along X
   dpdx = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i-1,j,QPRES)+0.5d0*q(i-2,j,QPRES))/dx
   dudx = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i-1,j,QU)+0.5d0*q(i-2,j,QU))/dx
   dvdx = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i-1,j,QV)+0.5d0*q(i-2,j,QV))/dx
   drhodx = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i-1,j,QRHO)+0.5d0*q(i-2,j,QRHO))/dx
   
   ! Test BCs
   if ((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j+1) == Outflow)) then
   ! This is the case when the upper right corner is outflow/outflow

     ! LODI waves along X
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx
     L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx)

     ! LODI waves along Y
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy) 
     M4 = (q(i,j,QV)+qaux(i,j,QC))* (dpdy + (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 

     ! Compute transverse terms
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     T1 = -0.5d0*M4 + (q(i,j,QRHO)*qaux(i,j,QC))*M2
     S1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))

     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(2))
     T2 = -0.5d0*L4 + (q(i,j,QRHO)*qaux(i,j,QC))*L3
     S2 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))
 
     S2 = S2 + (1.0d0-beta)*T2 
     S1 = S1 + (1.0d0-beta)*T1
   
     L1 = -2.0d0*(S2*beta + 2.0d0*S1 - S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     M1 = -2.0d0*(S1*(beta-1.0d0) + 2.0d0*S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)

   elseif (((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j+1) == NoSlipWall)) .or. &
           ((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j+1) ==   SlipWall))) then

   ! This is the case when the upper right corner is wall(y)/outflow(x)
     ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j,QPRES)-q(i,j-1,QPRES))/(2.0d0*dy)
     dudy = (-q(i,j,QU)-q(i,j-1,QU))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == NoSlipWall) dvdy = (-q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == SlipWall) dvdy = (q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     drhody = (q(i,j,QRHO)-q(i,j-1,QRHO))/(2.0d0*dy)
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     
     ! LODI waves along X
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     L1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir)) - ((1.0d0 - beta)*T1)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx) 
     L3 = q(i,j,QU) * dvdx 
     L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx) 
     
   elseif (((x_bcMask(i+1,j)  == NoSlipWall) .and. (y_bcMask(i,j+1) == Outflow)) .or. &
           ((x_bcMask(i+1,j)  == SlipWall) .and. (y_bcMask(i,j+1) == Outflow))) then

   ! This is the case when the upper right corner is outflow(y)/wall(x)
      call bl_error("NSCBC not implemented for upper right corner")
      
   elseif (((x_bcMask(i+1,j)  == SlipWall) .or. (x_bcMask(i+1,j)  == NoSlipWall)) .and.  & 
           ((y_bcMask(i,j+1) == SlipWall) .or. (y_bcMask(i,j+1) == NoSlipWall))) then
   ! This is the case when the upper right corner is wall(y)/wall(x)
     ! Values long Y will be computed by mirror functions below

   else
     call bl_error("NSCBC not implemented for upper right corner")
   endif

   if (q(i,j,QV) == 0.0d0) then
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = 0.0d0
       M3 = 0.0d0
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   else
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = M2 / q(i,j,QV)
       M3 = M3 / q(i,j,QV)
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   endif   

   if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   else
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   endif

   ! Along X
   ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)

     ! 2nd order
     q(i+1,j,QU)    = q(i-1,j,QU) + 2.0d0*dx*dudx
     q(i+1,j,QV)    = q(i-1,j,QV) + 2.0d0*dx*dvdx
     q(i+1,j,QRHO)  = q(i-1,j,QRHO)  + 2.0d0*dx*drhodx
     q(i+1,j,QPRES) = q(i-1,j,QPRES) + 2.0d0*dx*dpdx

     !----------------
     q(i+2,j,QU)    = -2.0d0*q(i-1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i+1,j,QU) - 6.0d0*dx*dudx
     q(i+2,j,QV)    = -2.0d0*q(i-1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i+1,j,QV) - 6.0d0*dx*dvdx
     q(i+2,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i+1,j,QRHO) - 6.0d0*dx*drhodx
     q(i+2,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i+1,j,QPRES) - 6.0d0*dx*dpdx

     q(i+3,j,QU)    = 3.0d0*q(i-1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i+1,j,QU) + 6.0d0*q(i+2,j,QU) + 12.0d0*dx*dudx
     q(i+3,j,QV)    = 3.0d0*q(i-1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i+1,j,QV) + 6.0d0*q(i+2,j,QV) + 12.0d0*dx*dvdx
     q(i+3,j,QRHO)    = 3.0d0*q(i-1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i+1,j,QRHO) + 6.0d0*q(i+2,j,QRHO) + 12.0d0*dx*drhodx
     q(i+3,j,QPRES)    = 3.0d0*q(i-1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i+1,j,QPRES) + 6.0d0*q(i+2,j,QPRES) + 12.0d0*dx*dpdx

     q(i+4,j,QU)    = -2.0d0*q(i-1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i+1,j,QU) - 12.0d0*q(i+2,j,QU)  &
                     + 4.0d0*q(i+3,j,QU) - 12.0d0*dx*dudx
     q(i+4,j,QV)    = -2.0d0*q(i-1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i+1,j,QV) - 12.0d0*q(i+2,j,QV) &
                       + 4.0d0*q(i+3,j,QV) - 12.0d0*dx*dvdx
     q(i+4,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i+1,j,QRHO) - 12.0d0*q(i+2,j,QRHO) &
                     + 4.0d0*q(i+3,j,QRHO) - 12.0d0*dx*drhodx
     q(i+4,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i+1,j,QPRES) - 12.0d0*q(i+2,j,QPRES) &
                     + 4.0d0*q(i+3,j,QPRES) - 12.0d0*dx*dpdx

     if ((x_bcMask(i+1,j) .eq. NoSlipWall).or.(x_bcMask(i+1,j) .eq. SlipWall)) then
     
       if (x_bcMask(i+1,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (x_bcMask(i+1,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i+1,j,QU)    = -q(i,j,QU)
       q(i+2,j,QU)    = -q(i-1,j,QU)
       q(i+3,j,QU)    = -q(i-2,j,QU)
       q(i+4,j,QU)    = -q(i-3,j,QU)

       q(i+1,j,QV)    = wall_sign*q(i,j,QV)
       q(i+2,j,QV)    = wall_sign*q(i-1,j,QV)
       q(i+3,j,QV)    = wall_sign*q(i-2,j,QV)
       q(i+4,j,QV)    = wall_sign*q(i-3,j,QV)

       q(i+1,j,QRHO)  = q(i,j,QRHO)
       q(i+2,j,QRHO)  = q(i-1,j,QRHO)
       q(i+3,j,QRHO)  = q(i-2,j,QRHO)
       q(i+4,j,QRHO)  = q(i-3,j,QRHO)

       q(i+1,j,QPRES)  = q(i,j,QPRES)
       q(i+2,j,QPRES)  = q(i-1,j,QPRES)
       q(i+3,j,QPRES)  = q(i-2,j,QPRES)
       q(i+4,j,QPRES)  = q(i-3,j,QPRES)
     
     end if

   ! Along Y
   ! Compute new spatial derivative
     drhody = (M3 + 0.5d0*(M1 + M4))/(qaux(i,j,QC)**2.0d0)
     dvdy   = (M4-M1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = M2
     dpdy   = 0.5d0*(M1+M4)

     ! Update ghost cells

     q(i,j+1,QU)    = q(i,j-1,QU) + 2.0d0*dy*dudy
     q(i,j+1,QV)    = q(i,j-1,QV) + 2.0d0*dy*dvdy
     q(i,j+1,QRHO)  = q(i,j-1,QRHO)  + 2.0d0*dy*drhody
     q(i,j+1,QPRES) = q(i,j-1,QPRES) + 2.0d0*dy*dpdy

     !----------------
     q(i,j+2,QU)    = -2.0d0*q(i,j-1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j+1,QU) - 6.0d0*dy*dudy
     q(i,j+2,QV)    = -2.0d0*q(i,j-1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j+1,QV) - 6.0d0*dy*dvdy
     q(i,j+2,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j+1,QRHO) - 6.0d0*dy*drhody
     q(i,j+2,QPRES) = -2.0d0*q(i,j-1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j+1,QPRES) - 6.0d0*dy*dpdy

     q(i,j+3,QU)    = 3.0d0*q(i,j-1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j+1,QU) + 6.0d0*q(i,j+2,QU) + 12.0d0*dy*dudy
     q(i,j+3,QV)    = 3.0d0*q(i,j-1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j+1,QV) + 6.0d0*q(i,j+2,QV) + 12.0d0*dy*dvdy
     q(i,j+3,QRHO)  = 3.0d0*q(i,j-1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j+1,QRHO) + 6.0d0*q(i,j+2,QRHO) + 12.0d0*dy*drhody
     q(i,j+3,QPRES) = 3.0d0*q(i,j-1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j+1,QPRES) + 6.0d0*q(i,j+2,QPRES) + 12.0d0*dy*dpdy

     q(i,j+4,QU)    = -2.0d0*q(i,j-1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j+1,QU) - 12.0d0*q(i,j+2,QU)  &
                     + 4.0d0*q(i,j+3,QU) - 12.0d0*dy*dudy
     q(i,j+4,QV)    = -2.0d0*q(i,j-1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j+1,QV) - 12.0d0*q(i,j+2,QV) &
                     + 4.0d0*q(i,j+3,QV) - 12.0d0*dy*dvdy
     q(i,j+4,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j+1,QRHO) - 12.0d0*q(i,j+2,QRHO) &
                     + 4.0d0*q(i,j+3,QRHO) - 12.0d0*dy*drhody
     q(i,j+4,QPRES) = -2.0d0*q(i,j-1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j+1,QPRES) - 12.0d0*q(i,j+2,QPRES) &
                     + 4.0d0*q(i,j+3,QPRES) - 12.0d0*dy*dpdy
 
     if ((y_bcMask(i,j+1) .eq. NoSlipWall).or.(y_bcMask(i,j+1) .eq. SlipWall)) then
     
       if (y_bcMask(i,j+1) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (y_bcMask(i,j+1) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i,j+1,QU)    = wall_sign*q(i,j,QU)
       q(i,j+2,QU)    = wall_sign*q(i,j-1,QU)
       q(i,j+3,QU)    = wall_sign*q(i,j-2,QU)
       q(i,j+4,QU)    = wall_sign*q(i,j-3,QU)

       q(i,j+1,QV)    = -q(i,j,QV)
       q(i,j+2,QV)    = -q(i,j-1,QV)
       q(i,j+3,QV)    = -q(i,j-2,QV)
       q(i,j+4,QV)    = -q(i,j-3,QV)

       q(i,j+1,QRHO)  = q(i,j,QRHO)
       q(i,j+2,QRHO)  = q(i,j-1,QRHO)
       q(i,j+3,QRHO)  = q(i,j-2,QRHO)
       q(i,j+4,QRHO)  = q(i,j-3,QRHO)

       q(i,j+1,QPRES)  = q(i,j,QPRES)
       q(i,j+2,QPRES)  = q(i,j-1,QPRES)
       q(i,j+3,QPRES)  = q(i,j-2,QPRES)
       q(i,j+4,QPRES)  = q(i,j-3,QPRES)
     
     end if
 
     do j = domhi(2)+1,domhi(2)+4,1
       do i = domhi(1)+1,domhi(1)+4,1

         eos_state % p        = q(i,j,QPRES )
         eos_state % rho      = q(i,j,QRHO  )
         eos_state % massfrac = q(i,j,QFS:QFS+nspec-1)
         eos_state % aux      = q(i,j,QFX:QFX+naux-1)

         call eos_rp(eos_state)
         q(i,j,QTEMP)  = eos_state % T
         q(i,j,QREINT) = eos_state % e * q(i,j,QRHO)
         q(i,j,QGAME)  = q(i,j,QPRES) / q(i,j,QREINT) + ONE
       
         qaux(i,j,QDPDR)  = eos_state % dpdr_e
         qaux(i,j,QDPDE)  = eos_state % dpde
         qaux(i,j,QGAMC)  = eos_state % gam1
         qaux(i,j,QC   )  = eos_state % cs
         qaux(i,j,QCSML)  = max(small, small * qaux(i,j,QC))

         uin(i,j,URHO )  = eos_state % rho 
         uin(i,j,UMX  )  = q(i,j,QU ) * eos_state % rho 
         uin(i,j,UMY  )  = q(i,j,QV ) * eos_state % rho 
         uin(i,j,UMZ  ) = 0.0d0
         uin(i,j,UEINT) = eos_state % rho   *  eos_state % e
         uin(i,j,UEDEN) = eos_state % rho  &
            * (eos_state % e + 0.5d0 * (uin(i,j,UMX)**2 + uin(i,j,UMY)**2))
         uin(i,j,UTEMP) = eos_state % T
         do n=1, nspec
           uin(i,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
         end do
         
       enddo
     enddo

 endif
 
 
 !--------------------------------------------------------------------------   
 ! Bottom right corner (Lx,0)
 ! phi = 1, psi = 4
 !--------------------------------------------------------------------------

 if ((q_hi(1) > domhi(1)) .and. (q_lo(2) < domlo(2))) then

   i = domhi(1)
   j = domlo(2)
  
   x   = (dble(i)+HALF)*dx
   y   = (dble(j)+HALF)*dy
  
   ! Calling user target BC values
   ! right face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,-1,.false.,bc_type,bc_params)
   x_bcMask(i+1,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i+1,j-1,1) = bc_type

   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta = bc_params(5)

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Xdir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Xdir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Xdir = U_ext(UTEMP)
   INLET_PRESSURE_Xdir = eos_state%p

   ! bottom face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,1,.false.,bc_type,bc_params)
   y_bcMask(i,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i+1,j-1,1) = bc_type
   
   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta = bc_params(5)

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Ydir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Ydir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Ydir = U_ext(UTEMP)
   INLET_PRESSURE_Ydir = eos_state%p

   mach_local_corner = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)

   ! Derivatives along Y
   dpdy = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i,j+1,QPRES)-0.5d0*q(i,j+2,QPRES))/dy
   dudy = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i,j+1,QU)-0.5d0*q(i,j+2,QU))/dy
   dvdy = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i,j+1,QV)-0.5d0*q(i,j+2,QV))/dy
   drhody = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i,j+1,QRHO)-0.5d0*q(i,j+2,QRHO))/dy

   ! Derivative along X
   dpdx = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i-1,j,QPRES)+0.5d0*q(i-2,j,QPRES))/dx
   dudx = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i-1,j,QU)+0.5d0*q(i-2,j,QU))/dx
   dvdx = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i-1,j,QV)+0.5d0*q(i-2,j,QV))/dx
   drhodx = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i-1,j,QRHO)+0.5d0*q(i-2,j,QRHO))/dx
   
   ! Test BCs
   if ((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j) == Outflow)) then
   ! This is the case when the bottom right corner is outflow/outflow

     ! LODI waves along X
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx
     L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx)

     ! LODI waves along Y
     M1 = (q(i,j,QV)-qaux(i,j,QC))* (dpdy - (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy) 

     ! Compute transverse terms
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     T1 = -0.5d0*M1 + (q(i,j,QRHO)*qaux(i,j,QC))*M2
     S1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Ydir))

     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(2))
     T2 = -0.5d0*L4 - (q(i,j,QRHO)*qaux(i,j,QC))*L3
     S2 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Ydir))
 
     S2 = S2 + (1.0d0-beta)*T2
     S1 = S1 + (1.0d0-beta)*T1
   
     L1 = -2.0d0*(S2*beta + 2.0d0*S1 - S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     M4 = -2.0d0*(S1*(beta-1.0d0) + 2.0d0*S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     
   elseif (((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j) == NoSlipWall)) .or. &
           ((x_bcMask(i+1,j) == Outflow) .and. (y_bcMask(i,j) == SlipWall))) then
              
   ! This is the case when the bottom right corner is wall(y)/outflow(x)
     ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j+1,QPRES)-q(i,j,QPRES))/(2.0d0*dy)
     dudy = (q(i,j+1,QU)+q(i,j,QU))/(2.0d0*dy)
     if (y_bcMask(i,j) == NoSlipWall) dvdy = (q(i,j+1,QV)+q(i,j,QV))/(2.0d0*dy)
     if (y_bcMask(i,j) == SlipWall) dvdy = (q(i,j+1,QV)-q(i,j,QV))/(2.0d0*dy)
     drhody = (q(i,j+1,QRHO)-q(i,j,QRHO))/(2.0d0*dy)
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     
     ! LODI waves along X
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     L1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir)) - ((1.0d0 - beta)*T1)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx) 
     L3 = q(i,j,QU) * dvdx 
     L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx) 
     
   elseif (((x_bcMask(i+1,j) == NoSlipWall) .and. (y_bcMask(i,j) == Outflow)) .or. &
           ((x_bcMask(i+1,j) == SlipWall) .and. (y_bcMask(i,j) == Outflow))) then
              
   ! This is the case when the bottom right corner is outflow(y)/wall(x)
     call bl_error("NSCBC not implemented for bottom right corner")

   elseif (((x_bcMask(i+1,j) == SlipWall) .or. (x_bcMask(i+1,j) == NoSlipWall)) .and.  & 
           ((y_bcMask(i,j) == SlipWall) .or. (y_bcMask(i,j) == NoSlipWall))) then
     ! This is the case when the bottom right corner is wall(y)/wall(x)
      ! Values long Y will be computed by mirror functions below

   else
     call bl_error("NSCBC not implemented for bottom right corner")
   endif  
     
   if (q(i,j,QV) == 0.0d0) then
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = 0.0d0
       M3 = 0.0d0
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   else
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = M2 / q(i,j,QV)
       M3 = M3 / q(i,j,QV)
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   endif   

   if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   else
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   endif

   ! Along X
   ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)

     ! 2nd order
     q(i+1,j,QU)    = q(i-1,j,QU) + 2.0d0*dx*dudx
     q(i+1,j,QV)    = q(i-1,j,QV) + 2.0d0*dx*dvdx
     q(i+1,j,QRHO)  = q(i-1,j,QRHO)  + 2.0d0*dx*drhodx
     q(i+1,j,QPRES) = q(i-1,j,QPRES) + 2.0d0*dx*dpdx

     !----------------
     q(i+2,j,QU)    = -2.0d0*q(i-1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i+1,j,QU) - 6.0d0*dx*dudx
     q(i+2,j,QV)    = -2.0d0*q(i-1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i+1,j,QV) - 6.0d0*dx*dvdx
     q(i+2,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i+1,j,QRHO) - 6.0d0*dx*drhodx
     q(i+2,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i+1,j,QPRES) - 6.0d0*dx*dpdx

     q(i+3,j,QU)    = 3.0d0*q(i-1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i+1,j,QU) + 6.0d0*q(i+2,j,QU) + 12.0d0*dx*dudx
     q(i+3,j,QV)    = 3.0d0*q(i-1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i+1,j,QV) + 6.0d0*q(i+2,j,QV) + 12.0d0*dx*dvdx
     q(i+3,j,QRHO)    = 3.0d0*q(i-1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i+1,j,QRHO) + 6.0d0*q(i+2,j,QRHO) + 12.0d0*dx*drhodx
     q(i+3,j,QPRES)    = 3.0d0*q(i-1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i+1,j,QPRES) + 6.0d0*q(i+2,j,QPRES) + 12.0d0*dx*dpdx

     q(i+4,j,QU)    = -2.0d0*q(i-1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i+1,j,QU) - 12.0d0*q(i+2,j,QU)  &
                     + 4.0d0*q(i+3,j,QU) - 12.0d0*dx*dudx
     q(i+4,j,QV)    = -2.0d0*q(i-1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i+1,j,QV) - 12.0d0*q(i+2,j,QV) &
                       + 4.0d0*q(i+3,j,QV) - 12.0d0*dx*dvdx
     q(i+4,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i+1,j,QRHO) - 12.0d0*q(i+2,j,QRHO) &
                     + 4.0d0*q(i+3,j,QRHO) - 12.0d0*dx*drhodx
     q(i+4,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i+1,j,QPRES) - 12.0d0*q(i+2,j,QPRES) &
                     + 4.0d0*q(i+3,j,QPRES) - 12.0d0*dx*dpdx

     if ((x_bcMask(i+1,j) .eq. NoSlipWall).or.(x_bcMask(i+1,j) .eq. SlipWall)) then

       if (x_bcMask(i+1,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (x_bcMask(i+1,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i+1,j,QU)    = -q(i,j,QU)
       q(i+2,j,QU)    = -q(i-1,j,QU)
       q(i+3,j,QU)    = -q(i-2,j,QU)
       q(i+4,j,QU)    = -q(i-3,j,QU)

       q(i+1,j,QV)    = wall_sign*q(i,j,QV)
       q(i+2,j,QV)    = wall_sign*q(i-1,j,QV)
       q(i+3,j,QV)    = wall_sign*q(i-2,j,QV)
       q(i+4,j,QV)    = wall_sign*q(i-3,j,QV)

       q(i+1,j,QRHO)  = q(i,j,QRHO)
       q(i+2,j,QRHO)  = q(i-1,j,QRHO)
       q(i+3,j,QRHO)  = q(i-2,j,QRHO)
       q(i+4,j,QRHO)  = q(i-3,j,QRHO)

       q(i+1,j,QPRES)  = q(i,j,QPRES)
       q(i+2,j,QPRES)  = q(i-1,j,QPRES)
       q(i+3,j,QPRES)  = q(i-2,j,QPRES)
       q(i+4,j,QPRES)  = q(i-3,j,QPRES)

     end if


   ! Along Y
   ! Compute new spatial derivative
     drhody = (M3 + 0.5d0*(M1 + M4))/(qaux(i,j,QC)**2.0d0)
     dvdy   = (M4-M1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = M2
     dpdy   = 0.5d0*(M1+M4)

     ! Update ghost cells

     q(i,j-1,QU)    = q(i,j+1,QU) - 2.0d0*dy*dudy
     q(i,j-1,QV)    = q(i,j+1,QV) - 2.0d0*dy*dvdy 
     q(i,j-1,QRHO)  = q(i,j+1,QRHO)  - 2.0d0*dy*drhody
     q(i,j-1,QPRES) = q(i,j+1,QPRES) - 2.0d0*dy*dpdy

     !----------------
     q(i,j-2,QU)    = -2.0d0*q(i,j+1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j-1,QU) + 6.0d0*dy*dudy
     q(i,j-2,QV)    = -2.0d0*q(i,j+1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j-1,QV) + 6.0d0*dy*dvdy
     q(i,j-2,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j-1,QRHO) + 6.0d0*dy*drhody
     q(i,j-2,QPRES) = -2.0d0*q(i,j+1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j-1,QPRES) + 6.0d0*dy*dpdy

     q(i,j-3,QU)    = 3.0d0*q(i,j+1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j-1,QU) + 6.0d0*q(i,j-2,QU) - 12.0d0*dy*dudy
     q(i,j-3,QV)    = 3.0d0*q(i,j+1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j-1,QV) + 6.0d0*q(i,j-2,QV) - 12.0d0*dy*dvdy
     q(i,j-3,QRHO)  = 3.0d0*q(i,j+1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j-1,QRHO) + 6.0d0*q(i,j-2,QRHO) - 12.0d0*dy*drhody
     q(i,j-3,QPRES) = 3.0d0*q(i,j+1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j-1,QPRES) + 6.0d0*q(i,j-2,QPRES) - 12.0d0*dy*dpdy

     q(i,j-4,QU)    = -2.0d0*q(i,j+1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j-1,QU) - 12.0d0*q(i,j-2,QU)  &
                     + 4.0d0*q(i,j-3,QU) + 12.0d0*dy*dudy
     q(i,j-4,QV)    = -2.0d0*q(i,j+1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j-1,QV) - 12.0d0*q(i,j-2,QV) &
                     + 4.0d0*q(i,j-3,QV) + 12.0d0*dy*dvdy
     q(i,j-4,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j-1,QRHO) - 12.0d0*q(i,j-2,QRHO) &
                     + 4.0d0*q(i,j-3,QRHO) + 12.0d0*dy*drhody
     q(i,j-4,QPRES) = -2.0d0*q(i,j+1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j-1,QPRES) - 12.0d0*q(i,j-2,QPRES) &
                     + 4.0d0*q(i,j-3,QPRES) + 12.0d0*dy*dpdy
     
     if ((y_bcMask(i,j) .eq. NoSlipWall).or.(y_bcMask(i,j) .eq. SlipWall)) then

       if (y_bcMask(i,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (y_bcMask(i,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i,j-1,QU)    = wall_sign*q(i,j,QU)
       q(i,j-2,QU)    = wall_sign*q(i,j+1,QU)
       q(i,j-3,QU)    = wall_sign*q(i,j+2,QU)
       q(i,j-4,QU)    = wall_sign*q(i,j+3,QU)

       q(i,j-1,QV)    = -q(i,j,QV)
       q(i,j-2,QV)    = -q(i,j+1,QV)
       q(i,j-3,QV)    = -q(i,j+2,QV)
       q(i,j-4,QV)    = -q(i,j+3,QV)

       q(i,j-1,QRHO)  = q(i,j,QRHO)
       q(i,j-2,QRHO)  = q(i,j+1,QRHO)
       q(i,j-3,QRHO)  = q(i,j+2,QRHO)
       q(i,j-4,QRHO)  = q(i,j+3,QRHO)

       q(i,j-1,QPRES)  = q(i,j,QPRES)
       q(i,j-2,QPRES)  = q(i,j+1,QPRES)
       q(i,j-3,QPRES)  = q(i,j+2,QPRES)
       q(i,j-4,QPRES)  = q(i,j+3,QPRES)

     end if


     do j = domlo(2)-1,domlo(2)-4,-1
       do i = domhi(1)+1,domhi(1)+4,1
       
         eos_state % p        = q(i,j,QPRES )
         eos_state % rho      = q(i,j,QRHO  )
         eos_state % massfrac = q(i,j,QFS:QFS+nspec-1)
         eos_state % aux      = q(i,j,QFX:QFX+naux-1)

         call eos_rp(eos_state)
         q(i,j,QTEMP)  = eos_state % T
         q(i,j,QREINT) = eos_state % e * q(i,j,QRHO)
         q(i,j,QGAME)  = q(i,j,QPRES) / q(i,j,QREINT) + ONE

         qaux(i,j,QDPDR)  = eos_state % dpdr_e
         qaux(i,j,QDPDE)  = eos_state % dpde
         qaux(i,j,QGAMC)  = eos_state % gam1
         qaux(i,j,QC   )  = eos_state % cs
         qaux(i,j,QCSML)  = max(small, small * qaux(i,j,QC))

         uin(i,j,URHO )  = eos_state % rho
         uin(i,j,UMX  )  = q(i,j,QU ) * eos_state % rho
         uin(i,j,UMY  )  = q(i,j,QV ) * eos_state % rho
         uin(i,j,UMZ  ) = 0.0d0
         uin(i,j,UEINT) = eos_state % rho   *  eos_state % e
         uin(i,j,UEDEN) = eos_state % rho  &
            * (eos_state % e + 0.5d0 * (uin(i,j,UMX)**2 + uin(i,j,UMY)**2))
         uin(i,j,UTEMP) = eos_state % T
         do n=1, nspec
           uin(i,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
         end do

       enddo
     enddo

 endif
 
 !--------------------------------------------------------------------------   
 ! upper left corner (0,Ly)
 ! phi = 4, psi = 1
 !--------------------------------------------------------------------------

 if ((q_lo(1) < domlo(1)) .and. (q_hi(2) > domhi(2))) then

   i = domlo(1)
   j = domhi(2)
   
   x   = (dble(i)+HALF)*dx
   y   = (dble(j)+HALF)*dy
   
   ! Calling user target BC values
   ! left face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,1,.false.,bc_type,bc_params)
   x_bcMask(i,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i-1,j+1,1) = bc_type

   if (bc_type == Outflow) sigma_out = bc_params(6)
   relax_T = bc_params(1)
   relax_U = bc_params(2)
   relax_V = bc_params(3)
   beta =  bc_params(5)

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Xdir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Xdir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Xdir = U_ext(UTEMP)
   INLET_PRESSURE_Xdir = eos_state%p

   ! top face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,-1,.false.,bc_type,bc_params)
   y_bcMask(i,j+1) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i-1,j+1,1) = bc_type
   
   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta = bc_params(5)

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Ydir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Ydir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Ydir = U_ext(UTEMP)
   INLET_PRESSURE_Ydir = eos_state%p

   mach_local_corner = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)
 
   ! Derivatives along Y
   dpdy = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i,j-1,QPRES)+0.5d0*q(i,j-2,QPRES))/dy
   dudy = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i,j-1,QU)+0.5d0*q(i,j-2,QU))/dy
   dvdy = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i,j-1,QV)+0.5d0*q(i,j-2,QV))/dy
   drhody = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i,j-1,QRHO)+0.5d0*q(i,j-2,QRHO))/dy

   ! Derivative along X
   dpdx = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i+1,j,QPRES)-0.5d0*q(i+2,j,QPRES))/dx
   dudx = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i+1,j,QU)-0.5d0*q(i+2,j,QU))/dx
   dvdx = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i+1,j,QV)-0.5d0*q(i+2,j,QV))/dx
   drhodx = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i+1,j,QRHO)-0.5d0*q(i+2,j,QRHO))/dx
   
   ! Test BCs
   if ((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j+1) == Outflow)) then
   ! This is the case when the upper left corner is outflow/outflow

     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx

     ! LODI waves along Y
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy) 
     M4 = (q(i,j,QV)+qaux(i,j,QC))* (dpdy + (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 

     ! Compute transverse terms
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     T1 = -0.5d0*M4 - (q(i,j,QRHO)*qaux(i,j,QC))*M2
     S1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))

     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(2))
     T2 = -0.5d0*L1 + (q(i,j,QRHO)*qaux(i,j,QC))*L3
     S2 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))
 
     S2 = S2 + (1.0d0-beta)*T2
     S1 = S1 + (1.0d0-beta)*T1
   
     L4 = -2.0d0*(S2*beta + 2.0d0*S1 - S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     M1 = -2.0d0*(S1*(beta-1.0d0) + 2.0d0*S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     
   elseif ((x_bcMask(i,j) == Inflow) .and. (y_bcMask(i,j+1) == Outflow)) then
   ! This is the case when the upper left corner is inflow(x)/outflow(y)

     ! LODI waves along Y
     M1 = 0.0d0 ! Compatibility condition
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(2))
     
       M1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Ydir)) !- ((1.0d0 - beta )*T1)
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy)
     M4 = (q(i,j,QV)+qaux(i,j,QC))* (dpdy + (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 
     
     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
     L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE_Xdir) !- M3
     L3 = (relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY_Xdir))  !- (M4/(2.0d0*q(i,j,QRHO)*qaux(i,j,QC)))
     L4 = (relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_corner*mach_local_corner)/probhi(1)) * &
            (q(i,j,QU) - INLET_VX_Xdir))  ! - (q(i,j,QRHO)*qaux(i,j,QC)*M2) - (M4/2.0d0)
    
   
   elseif (((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j+1) == NoSlipWall)) .or. &
           ((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j+1) == SlipWall))) then

   ! This is the case when the upper right corner is wall(y)/outflow(x)
     ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j,QPRES)-q(i,j-1,QPRES))/(2.0d0*dy)
     dudy = (-q(i,j,QU)-q(i,j-1,QU))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == NoSlipWall) dvdy = (-q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == SlipWall) dvdy = (q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     drhody = (q(i,j,QRHO)-q(i,j-1,QRHO))/(2.0d0*dy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     ! LODI waves along X
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx
     L4 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir)) - ((1.0d0 - beta)*T4)
     
   elseif ((x_bcMask(i,j) == Inflow) .and. &
           ((y_bcMask(i,j+1) == SlipWall).or.(y_bcMask(i,j+1) == NoSlipWall))) then
           
   ! This is the case when the upper left corner is wall(y)/inflow(x)
   ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j,QPRES)-q(i,j-1,QPRES))/(2.0d0*dy)
     dudy = (-q(i,j,QU)-q(i,j-1,QU))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == NoSlipWall) dvdy = (-q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     if (y_bcMask(i,j+1) == SlipWall) dvdy = (q(i,j,QV)-q(i,j-1,QV))/(2.0d0*dy)
     drhody = (q(i,j,QRHO)-q(i,j-1,QRHO))/(2.0d0*dy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
     L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE_Xdir) - ((1.0d0 - beta)*T2)
     L3 = relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY_Xdir)  - ((1.0d0 - beta)*T3)
     L4 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_corner*mach_local_corner)/probhi(1)) * &
          (q(i,j,QU) - INLET_VX_Xdir)  - ((1.0d0 - beta)*T4)

   elseif (((x_bcMask(i,j) == SlipWall) .or. (x_bcMask(i,j) == NoSlipWall)) .and. &
           ((y_bcMask(i,j+1) == SlipWall) .or. (y_bcMask(i,j+1) == NoSlipWall))) then
     ! This is the case when the upper left corner is wall(y)/wall(x)
       ! Values long Y will be computed by mirror functions below
 
   else
     call bl_error("NSCBC not implemented for upper left corner")
   endif 

   if (q(i,j,QV) == 0.0d0) then
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = 0.0d0
       M3 = 0.0d0
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   else
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = M2 / q(i,j,QV)
       M3 = M3 / q(i,j,QV)
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   endif   

   if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   else
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   endif

   ! Along X
   ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)

     ! 2nd order
     q(i-1,j,QU)    = q(i+1,j,QU) - 2.0d0*dx*dudx
     q(i-1,j,QV)    = q(i+1,j,QV) - 2.0d0*dx*dvdx 
     q(i-1,j,QRHO)  = q(i+1,j,QRHO)  - 2.0d0*dx*drhodx
     q(i-1,j,QPRES) = q(i+1,j,QPRES) - 2.0d0*dx*dpdx
   
     !---------------- 
     q(i-2,j,QU)    = -2.0d0*q(i+1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i-1,j,QU) + 6.0d0*dx*dudx
     q(i-2,j,QV)    = -2.0d0*q(i+1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i-1,j,QV) + 6.0d0*dx*dvdx
     q(i-2,j,QRHO)    = -2.0d0*q(i+1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i-1,j,QRHO) + 6.0d0*dx*drhodx
     q(i-2,j,QPRES)    = -2.0d0*q(i+1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i-1,j,QPRES) + 6.0d0*dx*dpdx

     q(i-3,j,QU)    = 3.0d0*q(i+1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i-1,j,QU) + 6.0d0*q(i-2,j,QU) - 12.0d0*dx*dudx
     q(i-3,j,QV)    = 3.0d0*q(i+1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i-1,j,QV) + 6.0d0*q(i-2,j,QV) - 12.0d0*dx*dvdx
     q(i-3,j,QRHO)    = 3.0d0*q(i+1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i-1,j,QRHO) + 6.0d0*q(i-2,j,QRHO) - 12.0d0*dx*drhodx
     q(i-3,j,QPRES)    = 3.0d0*q(i+1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i-1,j,QPRES) + 6.0d0*q(i-2,j,QPRES) - 12.0d0*dx*dpdx

     q(i-4,j,QU)    = -2.0d0*q(i+1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i-1,j,QU) - 12.0d0*q(i-2,j,QU)  &
                     + 4.0d0*q(i-3,j,QU) + 12.0d0*dx*dudx
     q(i-4,j,QV)    = -2.0d0*q(i+1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i-1,j,QV) - 12.0d0*q(i-2,j,QV) &
                     + 4.0d0*q(i-3,j,QV) + 12.0d0*dx*dvdx
     q(i-4,j,QRHO)  = -2.0d0*q(i+1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i-1,j,QRHO) - 12.0d0*q(i-2,j,QRHO) &
                     + 4.0d0*q(i-3,j,QRHO) + 12.0d0*dx*drhodx
     q(i-4,j,QPRES) = -2.0d0*q(i+1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i-1,j,QPRES) - 12.0d0*q(i-2,j,QPRES) &
                     + 4.0d0*q(i-3,j,QPRES) + 12.0d0*dx*dpdx

     if ((x_bcMask(i,j) .eq. NoSlipWall).or.(x_bcMask(i,j) .eq. SlipWall)) then
     
       if (x_bcMask(i,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (x_bcMask(i,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i-1,j,QU)    = -q(i,j,QU)
       q(i-2,j,QU)    = -q(i+1,j,QU)
       q(i-3,j,QU)    = -q(i+2,j,QU)
       q(i-4,j,QU)    = -q(i+3,j,QU)

       q(i-1,j,QV)    = wall_sign*q(i,j,QV)
       q(i-2,j,QV)    = wall_sign*q(i+1,j,QV)
       q(i-3,j,QV)    = wall_sign*q(i+2,j,QV)
       q(i-4,j,QV)    = wall_sign*q(i+3,j,QV)
       
       q(i-1,j,QRHO)  = q(i,j,QRHO)
       q(i-2,j,QRHO)  = q(i+1,j,QRHO)
       q(i-3,j,QRHO)  = q(i+2,j,QRHO)
       q(i-4,j,QRHO)  = q(i+3,j,QRHO)
       
       q(i-1,j,QPRES)  = q(i,j,QPRES)
       q(i-2,j,QPRES)  = q(i+1,j,QPRES)
       q(i-3,j,QPRES)  = q(i+2,j,QPRES)
       q(i-4,j,QPRES)  = q(i+3,j,QPRES)
     
     end if

   ! Along Y
   ! Compute new spatial derivative
     drhody = (M3 + 0.5d0*(M1 + M4))/(qaux(i,j,QC)**2.0d0)
     dvdy   = (M4-M1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = M2
     dpdy   = 0.5d0*(M1+M4)

     ! Update ghost cells

     q(i,j+1,QU)    = q(i,j-1,QU) + 2.0d0*dy*dudy
     q(i,j+1,QV)    = q(i,j-1,QV) + 2.0d0*dy*dvdy
     q(i,j+1,QRHO)  = q(i,j-1,QRHO)  + 2.0d0*dy*drhody
     q(i,j+1,QPRES) = q(i,j-1,QPRES) + 2.0d0*dy*dpdy

     !----------------
     q(i,j+2,QU)    = -2.0d0*q(i,j-1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j+1,QU) - 6.0d0*dy*dudy
     q(i,j+2,QV)    = -2.0d0*q(i,j-1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j+1,QV) - 6.0d0*dy*dvdy
     q(i,j+2,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j+1,QRHO) - 6.0d0*dy*drhody
     q(i,j+2,QPRES) = -2.0d0*q(i,j-1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j+1,QPRES) - 6.0d0*dy*dpdy

     q(i,j+3,QU)    = 3.0d0*q(i,j-1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j+1,QU) + 6.0d0*q(i,j+2,QU) + 12.0d0*dy*dudy
     q(i,j+3,QV)    = 3.0d0*q(i,j-1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j+1,QV) + 6.0d0*q(i,j+2,QV) + 12.0d0*dy*dvdy
     q(i,j+3,QRHO)  = 3.0d0*q(i,j-1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j+1,QRHO) + 6.0d0*q(i,j+2,QRHO) + 12.0d0*dy*drhody
     q(i,j+3,QPRES) = 3.0d0*q(i,j-1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j+1,QPRES) + 6.0d0*q(i,j+2,QPRES) + 12.0d0*dy*dpdy

     q(i,j+4,QU)    = -2.0d0*q(i,j-1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j+1,QU) - 12.0d0*q(i,j+2,QU)  &
                     + 4.0d0*q(i,j+3,QU) - 12.0d0*dy*dudy
     q(i,j+4,QV)    = -2.0d0*q(i,j-1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j+1,QV) - 12.0d0*q(i,j+2,QV) &
                     + 4.0d0*q(i,j+3,QV) - 12.0d0*dy*dvdy
     q(i,j+4,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j+1,QRHO) - 12.0d0*q(i,j+2,QRHO) &
                     + 4.0d0*q(i,j+3,QRHO) - 12.0d0*dy*drhody
     q(i,j+4,QPRES) = -2.0d0*q(i,j-1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j+1,QPRES) - 12.0d0*q(i,j+2,QPRES) &
                     + 4.0d0*q(i,j+3,QPRES) - 12.0d0*dy*dpdy
                     
     if ((y_bcMask(i,j+1) .eq. NoSlipWall).or.(y_bcMask(i,j+1) .eq. SlipWall)) then
     
       if (y_bcMask(i,j+1) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (y_bcMask(i,j+1) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i,j+1,QU)    = wall_sign*q(i,j,QU)
       q(i,j+2,QU)    = wall_sign*q(i,j-1,QU)
       q(i,j+3,QU)    = wall_sign*q(i,j-2,QU)
       q(i,j+4,QU)    = wall_sign*q(i,j-3,QU)

       q(i,j+1,QV)    = -q(i,j,QV)
       q(i,j+2,QV)    = -q(i,j-1,QV)
       q(i,j+3,QV)    = -q(i,j-2,QV)
       q(i,j+4,QV)    = -q(i,j-3,QV)

       q(i,j+1,QRHO)  = q(i,j,QRHO)
       q(i,j+2,QRHO)  = q(i,j-1,QRHO)
       q(i,j+3,QRHO)  = q(i,j-2,QRHO)
       q(i,j+4,QRHO)  = q(i,j-3,QRHO)

       q(i,j+1,QPRES)  = q(i,j,QPRES)
       q(i,j+2,QPRES)  = q(i,j-1,QPRES)
       q(i,j+3,QPRES)  = q(i,j-2,QPRES)
       q(i,j+4,QPRES)  = q(i,j-3,QPRES)
     
     end if
 
     do j = domhi(2)+1,domhi(2)+4,1
       do i = domlo(1)-1,domlo(1)-4,-1
         
         eos_state % p        = q(i,j,QPRES )
         eos_state % rho      = q(i,j,QRHO  )
         eos_state % massfrac = q(i,j,QFS:QFS+nspec-1)
         eos_state % aux      = q(i,j,QFX:QFX+naux-1)

         call eos_rp(eos_state)
         q(i,j,QTEMP)  = eos_state % T
         q(i,j,QREINT) = eos_state % e * q(i,j,QRHO)
         q(i,j,QGAME)  = q(i,j,QPRES) / q(i,j,QREINT) + ONE

         qaux(i,j,QDPDR)  = eos_state % dpdr_e
         qaux(i,j,QDPDE)  = eos_state % dpde
         qaux(i,j,QGAMC)  = eos_state % gam1
         qaux(i,j,QC   )  = eos_state % cs
         qaux(i,j,QCSML)  = max(small, small * qaux(i,j,QC))

         uin(i,j,URHO )  = eos_state % rho
         uin(i,j,UMX  )  = q(i,j,QU ) * eos_state % rho
         uin(i,j,UMY  )  = q(i,j,QV ) * eos_state % rho
         uin(i,j,UMZ  ) = 0.0d0
         uin(i,j,UEINT) = eos_state % rho   *  eos_state % e
         uin(i,j,UEDEN) = eos_state % rho  &
            * (eos_state % e + 0.5d0 * (uin(i,j,UMX)**2 + uin(i,j,UMY)**2))
         uin(i,j,UTEMP) = eos_state % T
         do n=1, nspec
           uin(i,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
         end do
         
       enddo
     enddo

 endif
 
 !--------------------------------------------------------------------------   
 ! Bottom left corner (0,0)
 ! phi = 4, psi = 4
 !--------------------------------------------------------------------------

 if ((q_lo(1) < domlo(1)) .and. (q_lo(2) < domlo(2))) then

   i = domlo(1)
   j = domlo(2)

   x   = (dble(i)+HALF)*dx
   y   = (dble(j)+HALF)*dy

   ! Calling user target BC values
   ! left face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,1,.false.,bc_type,bc_params)
   x_bcMask(i,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i-1,j-1,1) = bc_type

   if (bc_type == Outflow) sigma_out = bc_params(6)
   relax_T = bc_params(1)
   relax_U = bc_params(2)
   relax_V = bc_params(3)
   beta_X =  bc_params(5)

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Xdir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Xdir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Xdir = U_ext(UTEMP)
   INLET_PRESSURE_Xdir = eos_state%p

   ! bottom face
   call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,1,.false.,bc_type,bc_params)
   y_bcMask(i,j) = bc_type
   !if ((bc_type == NoSlipWall).or.(bc_type == SlipWall)) bcMask(i-1,j-1,1) = bc_type
   
   if (bc_type == Outflow) sigma_out = bc_params(6)
   beta_Y = bc_params(5)
   beta = beta_X

   eos_state %  T = U_ext(UTEMP)
   eos_state %  rho = U_ext(URHO)
   eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
   call eos_rt(eos_state)
   INLET_VX_Ydir = U_ext(UMX)/U_ext(URHO)
   INLET_VY_Ydir = U_ext(UMY)/U_ext(URHO)
   INLET_TEMPERATURE_Ydir = U_ext(UTEMP)
   INLET_PRESSURE_Ydir = eos_state%p

   mach_local_corner = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)

   ! Derivatives along Y
   dpdy = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i,j+1,QPRES)-0.5d0*q(i,j+2,QPRES))/dy
   dudy = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i,j+1,QU)-0.5d0*q(i,j+2,QU))/dy
   dvdy = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i,j+1,QV)-0.5d0*q(i,j+2,QV))/dy
   drhody = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i,j+1,QRHO)-0.5d0*q(i,j+2,QRHO))/dy

   ! Derivative along X
   dpdx = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i+1,j,QPRES)-0.5d0*q(i+2,j,QPRES))/dx
   dudx = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i+1,j,QU)-0.5d0*q(i+2,j,QU))/dx
   dvdx = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i+1,j,QV)-0.5d0*q(i+2,j,QV))/dx
   drhodx = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i+1,j,QRHO)-0.5d0*q(i+2,j,QRHO))/dx
   
   ! Test BCs
   if ((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j) == Outflow)) then
   ! This is the case when the bottom left corner is outflow/outflow

     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx

     ! LODI waves along Y
     M1 = (q(i,j,QV)-qaux(i,j,QC))* (dpdy - (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy)  

     ! Compute transverse terms
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))

     T1 = -0.5d0*M1 - (q(i,j,QRHO)*qaux(i,j,QC))*M2
     S1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))

     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(2))
     T2 = -0.5d0*L1 - (q(i,j,QRHO)*qaux(i,j,QC))*L3
     S2 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir))
 
     S2 = S2 + (1.0d0-beta)*T2
     S1 = S1 + (1.0d0-beta)*T1
   
     L4 = -2.0d0*(S2*beta + 2.0d0*S1 - S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     M4 = -2.0d0*(S1*(beta-1.0d0) + 2.0d0*S2)/ ((beta*beta) - (2.0d0*beta) - 3.0d0)
     
   elseif ((x_bcMask(i,j) == Inflow) .and. (y_bcMask(i,j) == Outflow)) then
   ! This is the case when the bottom left corner is inflow(x)/outflow(y)

     ! LODI waves along Y
     M1 = (q(i,j,QV)-qaux(i,j,QC))* (dpdy - (q(i,j,QRHO)*qaux(i,j,QC))*dvdy) 
     M2 = q(i,j,QV) * dudy
     M3 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy)
     M4 = 0.0d0 ! Compatibility condition
     
     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
     L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE_Xdir) - M3
     L3 = (relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY_Xdir))  + (M1/(2.0d0*q(i,j,QRHO)*qaux(i,j,QC)))
     L4 = (relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_corner*mach_local_corner)/probhi(1)) * &
            (q(i,j,QU) - INLET_VX_Xdir))  - (M1/2.0d0) - (q(i,j,QRHO)*qaux(i,j,QC)*M2)
            
   elseif (((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j) == NoSlipWall)) .or. &
           ((x_bcMask(i,j) == Outflow) .and. (y_bcMask(i,j) == SlipWall))) then
           
   ! This is the case when the bottom left corner is wall(y)/outflow(x)
     ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j+1,QPRES)-q(i,j,QPRES))/(2.0d0*dy)
     dudy = (q(i,j+1,QU)+q(i,j,QU))/(2.0d0*dy)
     if (y_bcMask(i,j) == NoSlipWall) dvdy = (q(i,j+1,QV)+q(i,j,QV))/(2.0d0*dy)
     if (y_bcMask(i,j) == SlipWall) dvdy = (q(i,j+1,QV)-q(i,j,QV))/(2.0d0*dy)
     drhody = (q(i,j+1,QRHO)-q(i,j,QRHO))/(2.0d0*dy)
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     ! LODI waves along X
     Kout = sigma_out*(1.0d0 - (mach_local_corner**2.0d0))*qaux(i,j,QC)/(probhi(1))
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)
     L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
     L3 = q(i,j,QU) * dvdx
     L4 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE_Xdir)) - ((1.0d0 - beta)*T4)
     
   elseif ((x_bcMask(i,j) == Inflow) .and. &
           ((y_bcMask(i,j) == SlipWall).or.(y_bcMask(i,j) == NoSlipWall))) then

   ! This is the case when the bottom left corner is wall(y)/inflow(x)
   ! Values long Y will be computed by mirror functions below
     dpdy = (q(i,j+1,QPRES)-q(i,j,QPRES))/(2.0d0*dy)
     dvdy = (q(i,j+1,QV)+q(i,j,QV))/(2.0d0*dy)
     drhody = (q(i,j+1,QRHO)+q(i,j,QRHO))/(2.0d0*dy)
     if (y_bcMask(i,j) == NoSlipWall) dudy = (q(i,j+1,QU)+q(i,j,QU))/(2.0d0*dy)
     if (y_bcMask(i,j) == SlipWall) dudy = (q(i,j+1,QU)-q(i,j,QU))/(2.0d0*dy)
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T3 = ((q(i,j,QV)*dvdy))+(dpdy/q(i,j,QRHO))
     T2 = (q(i,j,QV)*((qaux(i,j,QC)*qaux(i,j,QC)*drhody)-dpdy)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dvdy) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     beta_X = mach_local_corner

     ! LODI waves along X
     L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
     L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE_Xdir) - ((1.0d0 - beta_X)*T2)

     L3 = relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY_Xdir)  !- ((1.0d0 - beta_X)*T3)
     L4 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_corner*mach_local_corner)/probhi(1)) * &
          (q(i,j,QU) - INLET_VX_Xdir)  - ((1.0d0 - beta_X)*T4)
            

  elseif (((x_bcMask(i,j) == SlipWall) .or. (x_bcMask(i,j) == NoSlipWall)) .and.  & 
           ((y_bcMask(i,j) == SlipWall) .or. (y_bcMask(i,j) == NoSlipWall))) then
   ! This is the case when the bottom left corner is wall(y)/wall(x)
     ! Values long Y will be computed by mirror functions below

   else
     call bl_error("NSCBC not implemented for bottom left corner")
   endif 

   if (q(i,j,QV) == 0.0d0) then
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = 0.0d0
       M3 = 0.0d0
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   else
       M1 = M1 / (q(i,j,QV)-qaux(i,j,QC))
       M2 = M2 / q(i,j,QV)
       M3 = M3 / q(i,j,QV)
       M4 = M4 / (q(i,j,QV)+qaux(i,j,QC))
   endif   

   if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   else
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
   endif

   ! Along X
   ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)

     ! 2nd order
     q(i-1,j,QU)    = q(i+1,j,QU) - 2.0d0*dx*dudx
     q(i-1,j,QV)    = q(i+1,j,QV) - 2.0d0*dx*dvdx 
     q(i-1,j,QRHO)  = q(i+1,j,QRHO)  - 2.0d0*dx*drhodx
     q(i-1,j,QPRES) = q(i+1,j,QPRES) - 2.0d0*dx*dpdx

     !---------------- 
     q(i-2,j,QU)    = -2.0d0*q(i+1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i-1,j,QU) + 6.0d0*dx*dudx
     q(i-2,j,QV)    = -2.0d0*q(i+1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i-1,j,QV) + 6.0d0*dx*dvdx
     q(i-2,j,QRHO)    = -2.0d0*q(i+1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i-1,j,QRHO) + 6.0d0*dx*drhodx
     q(i-2,j,QPRES)    = -2.0d0*q(i+1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i-1,j,QPRES) + 6.0d0*dx*dpdx

     q(i-3,j,QU)    = 3.0d0*q(i+1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i-1,j,QU) + 6.0d0*q(i-2,j,QU) - 12.0d0*dx*dudx
     q(i-3,j,QV)    = 3.0d0*q(i+1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i-1,j,QV) + 6.0d0*q(i-2,j,QV) - 12.0d0*dx*dvdx
     q(i-3,j,QRHO)    = 3.0d0*q(i+1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i-1,j,QRHO) + 6.0d0*q(i-2,j,QRHO) - 12.0d0*dx*drhodx
     q(i-3,j,QPRES)    = 3.0d0*q(i+1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i-1,j,QPRES) + 6.0d0*q(i-2,j,QPRES) - 12.0d0*dx*dpdx

     q(i-4,j,QU)    = -2.0d0*q(i+1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i-1,j,QU) - 12.0d0*q(i-2,j,QU)  &
                     + 4.0d0*q(i-3,j,QU) + 12.0d0*dx*dudx
     q(i-4,j,QV)    = -2.0d0*q(i+1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i-1,j,QV) - 12.0d0*q(i-2,j,QV) &
                     + 4.0d0*q(i-3,j,QV) + 12.0d0*dx*dvdx
     q(i-4,j,QRHO)  = -2.0d0*q(i+1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i-1,j,QRHO) - 12.0d0*q(i-2,j,QRHO) &
                     + 4.0d0*q(i-3,j,QRHO) + 12.0d0*dx*drhodx
     q(i-4,j,QPRES) = -2.0d0*q(i+1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i-1,j,QPRES) - 12.0d0*q(i-2,j,QPRES) &
                     + 4.0d0*q(i-3,j,QPRES) + 12.0d0*dx*dpdx

     if ((x_bcMask(i,j) .eq. NoSlipWall).or.(x_bcMask(i,j) .eq. SlipWall)) then

       if (x_bcMask(i,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (x_bcMask(i,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i-1,j,QU)    = -q(i,j,QU)
       q(i-2,j,QU)    = -q(i+1,j,QU)
       q(i-3,j,QU)    = -q(i+2,j,QU)
       q(i-4,j,QU)    = -q(i+3,j,QU)

       q(i-1,j,QV)    = wall_sign*q(i,j,QV)
       q(i-2,j,QV)    = wall_sign*q(i+1,j,QV)
       q(i-3,j,QV)    = wall_sign*q(i+2,j,QV)
       q(i-4,j,QV)    = wall_sign*q(i+3,j,QV)

       q(i-1,j,QRHO)  = q(i,j,QRHO)
       q(i-2,j,QRHO)  = q(i+1,j,QRHO)
       q(i-3,j,QRHO)  = q(i+2,j,QRHO)
       q(i-4,j,QRHO)  = q(i+3,j,QRHO)

       q(i-1,j,QPRES)  = q(i,j,QPRES)
       q(i-2,j,QPRES)  = q(i+1,j,QPRES)
       q(i-3,j,QPRES)  = q(i+2,j,QPRES)
       q(i-4,j,QPRES)  = q(i+3,j,QPRES)

     end if


   ! Along Y
   ! Compute new spatial derivative
     drhody = (M3 + 0.5d0*(M1 + M4))/(qaux(i,j,QC)**2.0d0)
     dvdy   = (M4-M1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = M2
     dpdy   = 0.5d0*(M1+M4)

     ! Update ghost cells

     q(i,j-1,QU)    = q(i,j+1,QU) - 2.0d0*dy*dudy
     q(i,j-1,QV)    = q(i,j+1,QV) - 2.0d0*dy*dvdy 
     q(i,j-1,QRHO)  = q(i,j+1,QRHO)  - 2.0d0*dy*drhody
     q(i,j-1,QPRES) = q(i,j+1,QPRES) - 2.0d0*dy*dpdy

     !----------------
     q(i,j-2,QU)    = -2.0d0*q(i,j+1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j-1,QU) + 6.0d0*dy*dudy
     q(i,j-2,QV)    = -2.0d0*q(i,j+1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j-1,QV) + 6.0d0*dy*dvdy
     q(i,j-2,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j-1,QRHO) + 6.0d0*dy*drhody
     q(i,j-2,QPRES) = -2.0d0*q(i,j+1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j-1,QPRES) + 6.0d0*dy*dpdy

     q(i,j-3,QU)    = 3.0d0*q(i,j+1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j-1,QU) + 6.0d0*q(i,j-2,QU) - 12.0d0*dy*dudy
     q(i,j-3,QV)    = 3.0d0*q(i,j+1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j-1,QV) + 6.0d0*q(i,j-2,QV) - 12.0d0*dy*dvdy
     q(i,j-3,QRHO)  = 3.0d0*q(i,j+1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j-1,QRHO) + 6.0d0*q(i,j-2,QRHO) - 12.0d0*dy*drhody
     q(i,j-3,QPRES) = 3.0d0*q(i,j+1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j-1,QPRES) + 6.0d0*q(i,j-2,QPRES) - 12.0d0*dy*dpdy

     q(i,j-4,QU)    = -2.0d0*q(i,j+1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j-1,QU) - 12.0d0*q(i,j-2,QU)  &
                     + 4.0d0*q(i,j-3,QU) + 12.0d0*dy*dudy
     q(i,j-4,QV)    = -2.0d0*q(i,j+1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j-1,QV) - 12.0d0*q(i,j-2,QV) &
                     + 4.0d0*q(i,j-3,QV) + 12.0d0*dy*dvdy
     q(i,j-4,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j-1,QRHO) - 12.0d0*q(i,j-2,QRHO) &
                     + 4.0d0*q(i,j-3,QRHO) + 12.0d0*dy*drhody
     q(i,j-4,QPRES) = -2.0d0*q(i,j+1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j-1,QPRES) - 12.0d0*q(i,j-2,QPRES) &
                     + 4.0d0*q(i,j-3,QPRES) + 12.0d0*dy*dpdy

     if ((y_bcMask(i,j) .eq. NoSlipWall).or.(y_bcMask(i,j) .eq. SlipWall)) then

       if (y_bcMask(i,j) .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (y_bcMask(i,j) .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i,j-1,QU)    = wall_sign*q(i,j,QU)
       q(i,j-2,QU)    = wall_sign*q(i,j+1,QU)
       q(i,j-3,QU)    = wall_sign*q(i,j+2,QU)
       q(i,j-4,QU)    = wall_sign*q(i,j+3,QU)

       q(i,j-1,QV)    = -q(i,j,QV)
       q(i,j-2,QV)    = -q(i,j+1,QV)
       q(i,j-3,QV)    = -q(i,j+2,QV)
       q(i,j-4,QV)    = -q(i,j+3,QV)

       q(i,j-1,QRHO)  = q(i,j,QRHO)
       q(i,j-2,QRHO)  = q(i,j+1,QRHO)
       q(i,j-3,QRHO)  = q(i,j+2,QRHO)
       q(i,j-4,QRHO)  = q(i,j+3,QRHO)

       q(i,j-1,QPRES)  = q(i,j,QPRES)
       q(i,j-2,QPRES)  = q(i,j+1,QPRES)
       q(i,j-3,QPRES)  = q(i,j+2,QPRES)
       q(i,j-4,QPRES)  = q(i,j+3,QPRES)

     end if

     do j = domlo(2)-1,domlo(2)-4,-1
       do i = domlo(1)-1,domlo(1)-4,-1
       
         eos_state % p        = q(i,j,QPRES )
         eos_state % rho      = q(i,j,QRHO  )
         eos_state % massfrac = q(i,j,QFS:QFS+nspec-1)
         eos_state % aux      = q(i,j,QFX:QFX+naux-1)

         call eos_rp(eos_state)
         q(i,j,QTEMP)  = eos_state % T
         q(i,j,QREINT) = eos_state % e * q(i,j,QRHO)
         q(i,j,QGAME)  = q(i,j,QPRES) / q(i,j,QREINT) + ONE

         qaux(i,j,QDPDR)  = eos_state % dpdr_e
         qaux(i,j,QDPDE)  = eos_state % dpde
         qaux(i,j,QGAMC)  = eos_state % gam1
         qaux(i,j,QC   )  = eos_state % cs
         qaux(i,j,QCSML)  = max(small, small * qaux(i,j,QC))

         uin(i,j,URHO )  = eos_state % rho
         uin(i,j,UMX  )  = q(i,j,QU ) * eos_state % rho
         uin(i,j,UMY  )  = q(i,j,QV ) * eos_state % rho
         uin(i,j,UMZ  ) = 0.0d0
         uin(i,j,UEINT) = eos_state % rho   *  eos_state % e
         uin(i,j,UEDEN) = eos_state % rho  &
            * (eos_state % e + 0.5d0 * (uin(i,j,UMX)**2 + uin(i,j,UMY)**2))
         uin(i,j,UTEMP) = eos_state % T
         do n=1, nspec
           uin(i,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
         end do

       enddo
     enddo

 endif

endif ! flag_nscbc_isAnyPerio ) 

  
 !--------------------------------------------------------------------------   
 ! lower X
 !--------------------------------------------------------------------------
 !write(*,*) 'DEBUG IN THE NSCBC ROUTINE'
 if ((q_lo(1) < domlo(1)) .and. (physbc_lo(1) /= Interior)) then
   i = domlo(1)
!write(*,*) 'DEBUG IN THE LOWER X '
   do j = q_lo(2)+1,q_hi(2)-1

      if ( flag_nscbc_isAnyPerio == 0) then
        if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners
      endif
      
      x   = (dble(i)+HALF)*dx
      y   = (dble(j)+HALF)*dy

     !2nd order
     dpdx = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i+1,j,QPRES)-0.5d0*q(i+2,j,QPRES))/dx
     dudx = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i+1,j,QU)-0.5d0*q(i+2,j,QU))/dx
     dvdx = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i+1,j,QV)-0.5d0*q(i+2,j,QV))/dx
     drhodx = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i+1,j,QRHO)-0.5d0*q(i+2,j,QRHO))/dx
 
     ! Derivative along y
     ! 2nd order Central
     dpdy = (q(i,j+1,QPRES)-q(i,j-1,QPRES))/(2.0d0*dy)
     dudy = (q(i,j+1,QU)-q(i,j-1,QU))/(2.0d0*dy)
     dvdy = (q(i,j+1,QV)-q(i,j-1,QV))/(2.0d0*dy)
     drhody = (q(i,j+1,QRHO)-q(i,j-1,QRHO))/(2.0d0*dy)

     ! Compute transverse terms
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T2 = (q(i,j,QV)*((qaux(i,j,QC)*qaux(i,j,QC)*drhody)-dpdy)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dvdy) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T3 = ((q(i,j,QV)*dvdy))+(dpdy/q(i,j,QRHO))
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     ! Calling user target BC values 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,1,.false.,bc_type,bc_params)
     !write(*,*) 'DEBUG CALLING BCNORMAL ',i,j,bc_type
     if ((j < domlo(2)-1) .or. (j > domhi(2)+1)) then
       continue
     else
       x_bcMask(i,j) = bc_type
       !write(*,*) 'DEBUG CALLING x_bcMask(i,j) ',i,j,x_bcMask(i,j)
     endif
     
     !write(*,*) 'DEBUG bc_type,  ',i,j,bc_type
     !write(*,*) 'DEBUG bc_params ',i,j,bc_params
     
     eos_state %  T = U_ext(UTEMP)
     eos_state %  rho = U_ext(URHO)
     eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
     call eos_rt(eos_state)
     INLET_VX = U_ext(UMX)/U_ext(URHO)
     INLET_VY = U_ext(UMY)/U_ext(URHO)
     INLET_TEMPERATURE = U_ext(UTEMP)
     INLET_PRESSURE = eos_state%p

     mach_local_lo_x = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)
     
     ! Compute LODI equations
     if (bc_type == Inflow) then
 
       relax_T = bc_params(1)      
       relax_U = bc_params(2)
       relax_V = bc_params(3)
       beta =  bc_params(5)
           
       !write(*,*) 'I AM IN INFLOW',relax_T,relax_U,relax_V,beta
                  
       L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
       L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE) - ((1.0d0 - beta)*T2)
       L3 = relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY)  !- ((1.0d0 - beta)*T3)
       L4 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_lo_x*mach_local_lo_x)/probhi(1)) * &
            (q(i,j,QU) - INLET_VX)  - ((1.0d0 - beta)*T4)
            
     elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
       ! Values long Y will be computed by mirror functions below
       
     elseif (bc_type == Outflow) then
       
     ! We find that using a local Mach number gives better results for high Mach nb.
     ! This is in contradiction with Granet AIAA 2010
     ! However for low Mach number a surface averaged Mach number is much more better
     ! as reported in the paper of Granet
       sigma_out = bc_params(6)
       beta =  mach_local_lo_x 
       Kout = sigma_out*(1.0d0 - (mach_local_lo_x**2.0d0))*qaux(i,j,QC)/(probhi(1))

       L1 = (q(i,j,QU)-qaux(i,j,QC))* (dpdx - (q(i,j,QRHO)*qaux(i,j,QC))*dudx)
       L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx)
       L3 = q(i,j,QU) * dvdx
       L4 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE)) - ((1.0d0 - beta)*T4)
       
     else
       call bl_error("Error:: This BC is not yet implemented for lo_x in characteristic form")
     endif
 
     if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
     else       
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
     endif
     
      !write(*,*) 'DEBUG L waves ',i,j,L1,L2,L3,L4
      
     ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)  
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)
       
     ! Update ghost cells
     ! 2nd order
     q(i-1,j,QU)    = q(i+1,j,QU) - 2.0d0*dx*dudx
     q(i-1,j,QV)    = q(i+1,j,QV) - 2.0d0*dx*dvdx 
     q(i-1,j,QRHO)  = q(i+1,j,QRHO)  - 2.0d0*dx*drhodx
     q(i-1,j,QPRES) = q(i+1,j,QPRES) - 2.0d0*dx*dpdx
   
     !---------------- 
     q(i-2,j,QU)    = -2.0d0*q(i+1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i-1,j,QU) + 6.0d0*dx*dudx
     q(i-2,j,QV)    = -2.0d0*q(i+1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i-1,j,QV) + 6.0d0*dx*dvdx
     q(i-2,j,QRHO)    = -2.0d0*q(i+1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i-1,j,QRHO) + 6.0d0*dx*drhodx
     q(i-2,j,QPRES)    = -2.0d0*q(i+1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i-1,j,QPRES) + 6.0d0*dx*dpdx
 
     q(i-3,j,QU)    = 3.0d0*q(i+1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i-1,j,QU) + 6.0d0*q(i-2,j,QU) - 12.0d0*dx*dudx
     q(i-3,j,QV)    = 3.0d0*q(i+1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i-1,j,QV) + 6.0d0*q(i-2,j,QV) - 12.0d0*dx*dvdx
     q(i-3,j,QRHO)    = 3.0d0*q(i+1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i-1,j,QRHO) + 6.0d0*q(i-2,j,QRHO) - 12.0d0*dx*drhodx
     q(i-3,j,QPRES)    = 3.0d0*q(i+1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i-1,j,QPRES) + 6.0d0*q(i-2,j,QPRES) - 12.0d0*dx*dpdx
 
     q(i-4,j,QU)    = -2.0d0*q(i+1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i-1,j,QU) - 12.0d0*q(i-2,j,QU)  &
                     + 4.0d0*q(i-3,j,QU) + 12.0d0*dx*dudx
     q(i-4,j,QV)    = -2.0d0*q(i+1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i-1,j,QV) - 12.0d0*q(i-2,j,QV) &
                     + 4.0d0*q(i-3,j,QV) + 12.0d0*dx*dvdx
     q(i-4,j,QRHO)  = -2.0d0*q(i+1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i-1,j,QRHO) - 12.0d0*q(i-2,j,QRHO) &
                     + 4.0d0*q(i-3,j,QRHO) + 12.0d0*dx*drhodx
     q(i-4,j,QPRES) = -2.0d0*q(i+1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i-1,j,QPRES) - 12.0d0*q(i-2,j,QPRES) &
                     + 4.0d0*q(i-3,j,QPRES) + 12.0d0*dx*dpdx
 
     if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i-1,j,QU)    = -q(i,j,QU)
       q(i-2,j,QU)    = -q(i+1,j,QU)
       q(i-3,j,QU)    = -q(i+2,j,QU)
       q(i-4,j,QU)    = -q(i+3,j,QU)

       q(i-1,j,QV)    = wall_sign*q(i,j,QV)
       q(i-2,j,QV)    = wall_sign*q(i+1,j,QV)
       q(i-3,j,QV)    = wall_sign*q(i+2,j,QV)
       q(i-4,j,QV)    = wall_sign*q(i+3,j,QV)
       
       q(i-1,j,QRHO)  = q(i,j,QRHO)
       q(i-2,j,QRHO)  = q(i+1,j,QRHO)
       q(i-3,j,QRHO)  = q(i+2,j,QRHO)
       q(i-4,j,QRHO)  = q(i+3,j,QRHO)
       
       q(i-1,j,QPRES)  = q(i,j,QPRES)
       q(i-2,j,QPRES)  = q(i+1,j,QPRES)
       q(i-3,j,QPRES)  = q(i+2,j,QPRES)
       q(i-4,j,QPRES)  = q(i+3,j,QPRES)
     
     end if

 
     ! Recompute missing values thanks to EOS
     do hop=domlo(1)-1,-4,-1
     
       eos_state % p        = q(hop,j,QPRES )
       eos_state % rho      = q(hop,j,QRHO  )
       eos_state % massfrac = q(hop,j,QFS:QFS+nspec-1)
       eos_state % aux      = q(hop,j,QFX:QFX+naux-1)

       call eos_rp(eos_state)
       q(hop,j,QTEMP)  = eos_state % T
       q(hop,j,QREINT) = eos_state % e * q(hop,j,QRHO)
       q(hop,j,QGAME)  = q(hop,j,QPRES) / q(hop,j,QREINT) + ONE
       
       qaux(hop,j,QDPDR)  = eos_state % dpdr_e
       qaux(hop,j,QDPDE)  = eos_state % dpde
       qaux(hop,j,QGAMC)  = eos_state % gam1
       qaux(hop,j,QC   )  = eos_state % cs
       qaux(hop,j,QCSML)  = max(small, small * qaux(hop,j,QC))

       uin(hop,j,URHO )  = eos_state % rho 
       uin(hop,j,UMX  )  = q(hop,j,QU ) * eos_state % rho 
       uin(hop,j,UMY  )  = q(hop,j,QV ) * eos_state % rho 
       uin(hop,j,UMZ  ) = 0.0d0
       uin(hop,j,UEINT) = eos_state % rho   *  eos_state % e
       uin(hop,j,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (uin(hop,j,UMX)**2 + uin(hop,j,UMY)**2))
       uin(hop,j,UTEMP) = eos_state % T
       do n=1, nspec
          uin(hop,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do   
       
     enddo

   enddo
 end if
  
 !--------------------------------------------------------------------------   
 ! upper X
 !--------------------------------------------------------------------------
 
 if ((q_hi(1) > domhi(1)) .and. (physbc_hi(1) /= Interior)) then
   i = domhi(1)
   
   !write(*,*) 'DEBUG IN THE UPPER X '
   
   do j = q_lo(2)+1,q_hi(2)-1
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((j == domlo(2)) .or. (j == domhi(2))) cycle !Doing that to avoid ghost cells already filled by corners 
     endif
     
     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy
     
     !2nd order
     dpdx = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i-1,j,QPRES)+0.5d0*q(i-2,j,QPRES))/dx
     dudx = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i-1,j,QU)+0.5d0*q(i-2,j,QU))/dx
     dvdx = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i-1,j,QV)+0.5d0*q(i-2,j,QV))/dx
     drhodx = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i-1,j,QRHO)+0.5d0*q(i-2,j,QRHO))/dx
       
     ! Derivative along y
     ! 2nd order Central for interior point
     dpdy = (q(i,j+1,QPRES)-q(i,j-1,QPRES))/(2.0d0*dy)
     dudy = (q(i,j+1,QU)-q(i,j-1,QU))/(2.0d0*dy)
     dvdy = (q(i,j+1,QV)-q(i,j-1,QV))/(2.0d0*dy)
     drhody = (q(i,j+1,QRHO)-q(i,j-1,QRHO))/(2.0d0*dy)
  
     ! Compute transverse terms
     T1 = (q(i,j,QV)*(dpdy - q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T2 = (q(i,j,QV)*((qaux(i,j,QC)*qaux(i,j,QC)*drhody)-dpdy)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dvdy) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)
     T3 = ((q(i,j,QV)*dvdy))+(dpdy/q(i,j,QRHO))
     T4 = (q(i,j,QV)*(dpdy + q(i,j,QRHO)*qaux(i,j,QC)*dudy)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dvdy)

     ! Calling user target BC values 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,1,-1,.false.,bc_type,bc_params)
     if ((j < domlo(2)-1) .or. (j > domhi(2)+1)) then
       continue
     else
       x_bcMask(i+1,j) = bc_type
     endif

     eos_state %  T = U_ext(UTEMP)
     eos_state %  rho = U_ext(URHO)
     eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
     call eos_rt(eos_state)
     INLET_VX = U_ext(UMX)/U_ext(URHO)
     INLET_VY = U_ext(UMY)/U_ext(URHO)
     INLET_TEMPERATURE = U_ext(UTEMP)
     INLET_PRESSURE = eos_state%p

     mach_local_hi_x = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)

     ! Compute LODI equations                        
     if (bc_type .eq. Inflow) then
     
       relax_T = bc_params(1)      
       relax_U = bc_params(2)
       relax_V = bc_params(3)
       beta = bc_params(5)

       L1 = relax_U * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_hi_x*mach_local_hi_x)/probhi(1)) * &
            (q(i,j,QU) - INLET_VX)  - T1
       L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(1)) * (q(i,j,QTEMP) - INLET_TEMPERATURE) - T2
         
       L3 = relax_V * (qaux(i,j,QC)/probhi(1)) * (q(i,j,QV) - INLET_VY) ! - ((1.0d0 - beta)*T3)
       L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx)  
         
     elseif (bc_type .eq. Outflow) then
     
       ! We find that using a local Mach number gives better results for high Mach nb.
       ! This is in contradiction with Granet AIAA 2010
       ! However for low Mach number a surface averaged Mach number is much more better
       ! as reported in the paper of Granet
       sigma_out = bc_params(6)
       beta = mach_local_hi_x 

       Kout = sigma_out*(1.0d0 - (mach_local_hi_x**2.0d0))*qaux(i,j,QC)/(probhi(1))

       L1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE)) - ((1.0d0 - beta)*T1)
       L2 = q(i,j,QU) * ( ((qaux(i,j,QC)**2.0d0)*drhodx) - dpdx) 
       L3 = q(i,j,QU) * dvdx 
       L4 = (q(i,j,QU)+qaux(i,j,QC))* (dpdx + (q(i,j,QRHO)*qaux(i,j,QC))*dudx) 
   
     elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
       ! Values long Y will be computed by mirror functions below
       
     else
       call bl_error("Error:: This BC is not yet implemented for hi_x in characteristic form")
     endif
   
     if (q(i,j,QU) == 0.0d0) then
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
     else
       L1 = L1 / (q(i,j,QU)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QU)
       L3 = L3 / q(i,j,QU)
       L4 = L4 / (q(i,j,QU)+qaux(i,j,QC))
     endif
       
     ! Compute new spatial derivative
     drhodx = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)  
     dudx   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dvdx   = L3
     dpdx   = 0.5d0*(L1+L4)

     ! 2nd order
     q(i+1,j,QU)    = q(i-1,j,QU) + 2.0d0*dx*dudx
     q(i+1,j,QV)    = q(i-1,j,QV) + 2.0d0*dx*dvdx
     q(i+1,j,QRHO)  = q(i-1,j,QRHO)  + 2.0d0*dx*drhodx
     q(i+1,j,QPRES) = q(i-1,j,QPRES) + 2.0d0*dx*dpdx

     !----------------
     q(i+2,j,QU)    = -2.0d0*q(i-1,j,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i+1,j,QU) - 6.0d0*dx*dudx
     q(i+2,j,QV)    = -2.0d0*q(i-1,j,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i+1,j,QV) - 6.0d0*dx*dvdx
     q(i+2,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i+1,j,QRHO) - 6.0d0*dx*drhodx
     q(i+2,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i+1,j,QPRES) - 6.0d0*dx*dpdx

     q(i+3,j,QU)    = 3.0d0*q(i-1,j,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i+1,j,QU) + 6.0d0*q(i+2,j,QU) + 12.0d0*dx*dudx
     q(i+3,j,QV)    = 3.0d0*q(i-1,j,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i+1,j,QV) + 6.0d0*q(i+2,j,QV) + 12.0d0*dx*dvdx
     q(i+3,j,QRHO)  = 3.0d0*q(i-1,j,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i+1,j,QRHO) + 6.0d0*q(i+2,j,QRHO) + 12.0d0*dx*drhodx
     q(i+3,j,QPRES)  = 3.0d0*q(i-1,j,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i+1,j,QPRES) + 6.0d0*q(i+2,j,QPRES) + 12.0d0*dx*dpdx

     q(i+4,j,QU)    = -2.0d0*q(i-1,j,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i+1,j,QU) - 12.0d0*q(i+2,j,QU)  & 
                     + 4.0d0*q(i+3,j,QU) - 12.0d0*dx*dudx 
     q(i+4,j,QV)    = -2.0d0*q(i-1,j,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i+1,j,QV) - 12.0d0*q(i+2,j,QV) &
                       + 4.0d0*q(i+3,j,QV) - 12.0d0*dx*dvdx
     q(i+4,j,QRHO)  = -2.0d0*q(i-1,j,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i+1,j,QRHO) - 12.0d0*q(i+2,j,QRHO) &
                     + 4.0d0*q(i+3,j,QRHO) - 12.0d0*dx*drhodx
     q(i+4,j,QPRES) = -2.0d0*q(i-1,j,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i+1,j,QPRES) - 12.0d0*q(i+2,j,QPRES) &
                     + 4.0d0*q(i+3,j,QPRES) - 12.0d0*dx*dpdx

     if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if
       
       q(i+1,j,QU)    = -q(i,j,QU)
       q(i+2,j,QU)    = -q(i-1,j,QU)
       q(i+3,j,QU)    = -q(i-2,j,QU)
       q(i+4,j,QU)    = -q(i-3,j,QU)

       q(i+1,j,QV)    = wall_sign*q(i,j,QV)
       q(i+2,j,QV)    = wall_sign*q(i-1,j,QV)
       q(i+3,j,QV)    = wall_sign*q(i-2,j,QV)
       q(i+4,j,QV)    = wall_sign*q(i-3,j,QV)

       q(i+1,j,QRHO)  = q(i,j,QRHO)
       q(i+2,j,QRHO)  = q(i-1,j,QRHO)
       q(i+3,j,QRHO)  = q(i-2,j,QRHO)
       q(i+4,j,QRHO)  = q(i-3,j,QRHO)

       q(i+1,j,QPRES)  = q(i,j,QPRES)
       q(i+2,j,QPRES)  = q(i-1,j,QPRES)
       q(i+3,j,QPRES)  = q(i-2,j,QPRES)
       q(i+4,j,QPRES)  = q(i-3,j,QPRES)
     
     end if

     ! Recompute missing values thanks to EOS
     do hop= domhi(1)+1,domhi(1)+4,1
     
       eos_state % p        = q(hop,j,QPRES )
       eos_state % rho      = q(hop,j,QRHO  )
       eos_state % massfrac = q(hop,j,QFS:QFS+nspec-1)
       eos_state % aux      = q(hop,j,QFX:QFX+naux-1)

       call eos_rp(eos_state)
       q(hop,j,QTEMP)  = eos_state % T
       q(hop,j,QREINT) = eos_state % e * q(hop,j,QRHO)
       q(hop,j,QGAME)  = q(hop,j,QPRES) / q(hop,j,QREINT) + ONE
       
       qaux(hop,j,QDPDR)  = eos_state % dpdr_e
       qaux(hop,j,QDPDE)  = eos_state % dpde
       qaux(hop,j,QGAMC)  = eos_state % gam1
       qaux(hop,j,QC   )  = eos_state % cs
       qaux(hop,j,QCSML)  = max(small, small * qaux(hop,j,QC))

       uin(hop,j,URHO )  = eos_state % rho 
       uin(hop,j,UMX  )  = q(hop,j,QU ) * eos_state % rho 
       uin(hop,j,UMY  )  = q(hop,j,QV ) * eos_state % rho 
       uin(hop,j,UMZ  ) = 0.0d0
       uin(hop,j,UEINT) = eos_state % rho   *  eos_state % e
       uin(hop,j,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (uin(hop,j,UMX)**2 + uin(hop,j,UMY)**2))
       uin(hop,j,UTEMP) = eos_state % T
       do n=1, nspec
          uin(hop,j,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do 
     enddo

  enddo
endif



 !--------------------------------------------------------------------------   
 ! lower Y
 !--------------------------------------------------------------------------
 
 if ((q_lo(2) < domlo(2)) .and. (physbc_lo(2) /= Interior)) then
 
   j = domlo(2)
   
   !write(*,*) 'DEBUG IN THE LOWER Y '
   
   do i = q_lo(1)+1,q_hi(1)-1
   
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy

     !2nd order
     dpdy = ((-3.0d0/2.0d0)*q(i,j,QPRES)+2.0d0*q(i,j+1,QPRES)-0.5d0*q(i,j+2,QPRES))/dy
     dudy = ((-3.0d0/2.0d0)*q(i,j,QU)+2.0d0*q(i,j+1,QU)-0.5d0*q(i,j+2,QU))/dy
     dvdy = ((-3.0d0/2.0d0)*q(i,j,QV)+2.0d0*q(i,j+1,QV)-0.5d0*q(i,j+2,QV))/dy
     drhody = ((-3.0d0/2.0d0)*q(i,j,QRHO)+2.0d0*q(i,j+1,QRHO)-0.5d0*q(i,j+2,QRHO))/dy

     ! Derivative along x
     ! Central
     dpdx = (q(i+1,j,QPRES)-q(i-1,j,QPRES))/(2.0d0*dx)
     dudx = (q(i+1,j,QU)-q(i-1,j,QU))/(2.0d0*dx)
     dvdx = (q(i+1,j,QV)-q(i-1,j,QV))/(2.0d0*dx)
     drhodx = (q(i+1,j,QRHO)-q(i-1,j,QRHO))/(2.0d0*dx)
      
     ! Compute transverse terms
     T1 = (q(i,j,QU)*(dpdx - q(i,j,QRHO)*qaux(i,j,QC)*dvdx)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
     T2 = (q(i,j,QU)*((qaux(i,j,QC)*qaux(i,j,QC)*drhodx)-dpdx)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dudx) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
     T3 = ((q(i,j,QU)*dudx))+(dpdx/q(i,j,QRHO))
     T4 = (q(i,j,QU)*(dpdx + q(i,j,QRHO)*qaux(i,j,QC)*dvdx)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
                    
     ! Calling user target BC values 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,1,.false.,bc_type,bc_params)
     if ((i < domlo(1)-1) .or. (i > domhi(1)+1)) then
       continue
     else
       y_bcMask(i,j) = bc_type
     endif

     eos_state %  T = U_ext(UTEMP)
     eos_state %  rho = U_ext(URHO)
     eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
     call eos_rt(eos_state)
     INLET_VX = U_ext(UMX)/U_ext(URHO)
     INLET_VY = U_ext(UMY)/U_ext(URHO)
     INLET_TEMPERATURE = U_ext(UTEMP)
     INLET_PRESSURE = eos_state%p

     mach_local_lo_y = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)
     
     ! Compute LODI equations
     if (bc_type .eq. Inflow) then
     
       relax_T = bc_params(1)      
       relax_U = bc_params(2)
       relax_V = bc_params(3)
       beta =  bc_params(5)

       L1 = (q(i,j,QV)-qaux(i,j,QC))* (dpdy - (q(i,j,QRHO)*qaux(i,j,QC))*dvdy)  
       L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(2)) * (q(i,j,QTEMP) - INLET_TEMPERATURE) - T2
       
       L3 = relax_U * (qaux(i,j,QC)/probhi(2)) * (q(i,j,QU) - INLET_VX)  !- ((1.0d0 - beta)*T3)
       L4 = relax_V * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_lo_y*mach_local_lo_y)/probhi(2)) * &
            (q(i,j,QV) - INLET_VY) - ((1.0d0 - beta)*T4)
     
     elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
       ! Values long Y will be computed by mirror functions below
    
     elseif (bc_type .eq. Outflow) then
     
       sigma_out = bc_params(6)
       beta = mach_local_lo_y
       Kout = sigma_out*(1.0d0 - (mach_local_lo_y**2.0d0))*qaux(i,j,QC)/(probhi(2))
     
       L1 = (q(i,j,QV)-qaux(i,j,QC))* (dpdy - (q(i,j,QRHO)*qaux(i,j,QC))*dvdy)
       L2 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy)
       L3 = q(i,j,QV) * dudy
       L4 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE)) - ((1.0d0 - beta)*T4)
       
     else
       call bl_error("Error:: This BC is not yet implemented for lo_y in characteristic form")
     endif

     if (q(i,j,QV) == 0.0d0) then
       L1 = L1 / (q(i,j,QV)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QV)+qaux(i,j,QC))
     else
       L1 = L1 / (q(i,j,QV)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QV)
       L3 = L3 / q(i,j,QV)
       L4 = L4 / (q(i,j,QV)+qaux(i,j,QC))
     endif 
       
     ! Compute new spatial derivative
     drhody = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)  
     dvdy   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = L3
     dpdy   = 0.5d0*(L1+L4)
     
     ! Update ghost cells
       
     q(i,j-1,QU)    = q(i,j+1,QU) - 2.0d0*dy*dudy
     q(i,j-1,QV)    = q(i,j+1,QV) - 2.0d0*dy*dvdy 
     q(i,j-1,QRHO)  = q(i,j+1,QRHO)  - 2.0d0*dy*drhody
     q(i,j-1,QPRES) = q(i,j+1,QPRES) - 2.0d0*dy*dpdy

     !----------------
     q(i,j-2,QU)    = -2.0d0*q(i,j+1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j-1,QU) + 6.0d0*dy*dudy
     q(i,j-2,QV)    = -2.0d0*q(i,j+1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j-1,QV) + 6.0d0*dy*dvdy
     q(i,j-2,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j-1,QRHO) + 6.0d0*dy*drhody
     q(i,j-2,QPRES) = -2.0d0*q(i,j+1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j-1,QPRES) + 6.0d0*dy*dpdy
     
     q(i,j-3,QU)    = 3.0d0*q(i,j+1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j-1,QU) + 6.0d0*q(i,j-2,QU) - 12.0d0*dy*dudy
     q(i,j-3,QV)    = 3.0d0*q(i,j+1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j-1,QV) + 6.0d0*q(i,j-2,QV) - 12.0d0*dy*dvdy
     q(i,j-3,QRHO)  = 3.0d0*q(i,j+1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j-1,QRHO) + 6.0d0*q(i,j-2,QRHO) - 12.0d0*dy*drhody
     q(i,j-3,QPRES) = 3.0d0*q(i,j+1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j-1,QPRES) + 6.0d0*q(i,j-2,QPRES) - 12.0d0*dy*dpdy
     
     q(i,j-4,QU)    = -2.0d0*q(i,j+1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j-1,QU) - 12.0d0*q(i,j-2,QU)  &
                     + 4.0d0*q(i,j-3,QU) + 12.0d0*dy*dudy
     q(i,j-4,QV)    = -2.0d0*q(i,j+1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j-1,QV) - 12.0d0*q(i,j-2,QV) &
                     + 4.0d0*q(i,j-3,QV) + 12.0d0*dy*dvdy
     q(i,j-4,QRHO)  = -2.0d0*q(i,j+1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j-1,QRHO) - 12.0d0*q(i,j-2,QRHO) &
                     + 4.0d0*q(i,j-3,QRHO) + 12.0d0*dy*drhody
     q(i,j-4,QPRES) = -2.0d0*q(i,j+1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j-1,QPRES) - 12.0d0*q(i,j-2,QPRES) &
                     + 4.0d0*q(i,j-3,QPRES) + 12.0d0*dy*dpdy

     if ((bc_type .eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then

       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i,j-1,QU)    = wall_sign*q(i,j,QU)
       q(i,j-2,QU)    = wall_sign*q(i,j+1,QU)
       q(i,j-3,QU)    = wall_sign*q(i,j+2,QU)
       q(i,j-4,QU)    = wall_sign*q(i,j+3,QU)

       q(i,j-1,QV)    = -q(i,j,QV)
       q(i,j-2,QV)    = -q(i,j+1,QV)
       q(i,j-3,QV)    = -q(i,j+2,QV)
       q(i,j-4,QV)    = -q(i,j+3,QV)

       q(i,j-1,QRHO)  = q(i,j,QRHO)
       q(i,j-2,QRHO)  = q(i,j+1,QRHO)
       q(i,j-3,QRHO)  = q(i,j+2,QRHO)
       q(i,j-4,QRHO)  = q(i,j+3,QRHO)

       q(i,j-1,QPRES)  = q(i,j,QPRES)
       q(i,j-2,QPRES)  = q(i,j+1,QPRES)
       q(i,j-3,QPRES)  = q(i,j+2,QPRES)
       q(i,j-4,QPRES)  = q(i,j+3,QPRES)

     end if  

     ! Recompute missing values thanks to EOS
     do hop= domlo(2)-1,-4,-1
     
       eos_state % p        = q(i,hop,QPRES )
       eos_state % rho      = q(i,hop,QRHO  )
       eos_state % massfrac = q(i,hop,QFS:QFS+nspec-1)
       eos_state % aux      = q(i,hop,QFX:QFX+naux-1)

       call eos_rp(eos_state)
       q(i,hop,QTEMP)  = eos_state % T
       q(i,hop,QREINT) = eos_state % e * q(i,hop,QRHO)
       q(i,hop,QGAME)  = q(i,hop,QPRES) / q(i,hop,QREINT) + ONE
       
       qaux(i,hop,QDPDR)  = eos_state % dpdr_e
       qaux(i,hop,QDPDE)  = eos_state % dpde
       qaux(i,hop,QGAMC)  = eos_state % gam1
       qaux(i,hop,QC   )  = eos_state % cs
       qaux(i,hop,QCSML)  = max(small, small * qaux(i,hop,QC))

       uin(i,hop,URHO )  = eos_state % rho 
       uin(i,hop,UMX  )  = q(i,hop,QU ) * eos_state % rho 
       uin(i,hop,UMY  )  = q(i,hop,QV ) * eos_state % rho 
       uin(i,hop,UMZ  ) = 0.0d0
       uin(i,hop,UEINT) = eos_state % rho   *  eos_state % e
       uin(i,hop,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (uin(i,hop,UMX)**2 + uin(i,hop,UMY)**2))
       uin(i,hop,UTEMP) = eos_state % T
       do n=1, nspec
          uin(i,hop,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do   

     enddo
   
   enddo
 end if
     
!--------------------------------------------------------------------------   
! upper Y
!--------------------------------------------------------------------------

 if ((q_hi(2) > domhi(2)) .and. (physbc_hi(2) /= Interior)) then
 
   j = domhi(2)
   
   !write(*,*) 'DEBUG IN THE UPPER Y '
   
   do i = q_lo(1)+1,q_hi(1)-1
     
     if ( flag_nscbc_isAnyPerio == 0) then
       if ((i == domlo(1)) .or. (i == domhi(1))) cycle !Doing that to avoid ghost cells already filled by corners
     endif
     
     x   = (dble(i)+HALF)*dx
     y   = (dble(j)+HALF)*dy 
   
     !2nd order
     dpdy = ((3.0d0/2.0d0)*q(i,j,QPRES)-2.0d0*q(i,j-1,QPRES)+0.5d0*q(i,j-2,QPRES))/dy
     dudy = ((3.0d0/2.0d0)*q(i,j,QU)-2.0d0*q(i,j-1,QU)+0.5d0*q(i,j-2,QU))/dy
     dvdy = ((3.0d0/2.0d0)*q(i,j,QV)-2.0d0*q(i,j-1,QV)+0.5d0*q(i,j-2,QV))/dy
     drhody = ((3.0d0/2.0d0)*q(i,j,QRHO)-2.0d0*q(i,j-1,QRHO)+0.5d0*q(i,j-2,QRHO))/dy

     ! Derivative along x
     ! Central
     dpdx = (q(i+1,j,QPRES)-q(i-1,j,QPRES))/(2.0d0*dx)
     dudx = (q(i+1,j,QU)-q(i-1,j,QU))/(2.0d0*dx)
     dvdx = (q(i+1,j,QV)-q(i-1,j,QV))/(2.0d0*dx)
     drhodx = (q(i+1,j,QRHO)-q(i-1,j,QRHO))/(2.0d0*dx)
      
     ! Compute transverse terms
     T1 = (q(i,j,QU)*(dpdx - q(i,j,QRHO)*qaux(i,j,QC)*dvdx)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
     T2 = (q(i,j,QU)*((qaux(i,j,QC)*qaux(i,j,QC)*drhodx)-dpdx)) + &
          (qaux(i,j,QC)*qaux(i,j,QC)*q(i,j,QRHO)*dudx) - (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
     T3 = ((q(i,j,QU)*dudx))+(dpdx/q(i,j,QRHO))
     T4 = (q(i,j,QU)*(dpdx + q(i,j,QRHO)*qaux(i,j,QC)*dvdx)) + (qaux(i,j,QGAMC) * q(i,j,QPRES)*dudx)
     
     ! Calling user target BC values 
     call bcnormal([x,y,0.0d0],U_dummy,U_ext,2,-1,.false.,bc_type,bc_params)
     if ((i < domlo(1)-1) .or. (i > domhi(1)+1)) then
       continue
     else
       y_bcMask(i,j+1) = bc_type
     endif
     
     eos_state %  T = U_ext(UTEMP)
     eos_state %  rho = U_ext(URHO)
     eos_state % massfrac(1:nspec) = u_ext(UFS:UFS+nspec-1) / U_ext(URHO)
     call eos_rt(eos_state)
     INLET_VX = U_ext(UMX)/U_ext(URHO)
     INLET_VY = U_ext(UMY)/U_ext(URHO)
     INLET_TEMPERATURE = U_ext(UTEMP)
     INLET_PRESSURE = eos_state%p
     
     mach_local_hi_y = dsqrt(q(i,j,QU)**2.0d0 + q(i,j,QV)**2.0d0)/qaux(i,j,QC)
       
     ! Compute LODI equations
     if (bc_type .eq. Inflow) then
     
       relax_T = bc_params(1)      
       relax_U = bc_params(2)
       relax_V = bc_params(3)
       beta =  bc_params(5)

       ! Compute relax parameter
       L1 = relax_V * ((q(i,j,QRHO)*qaux(i,j,QC)**2.0d0)*(1.0d0-mach_local_hi_y*mach_local_hi_y)/probhi(2)) * &
            (q(i,j,QV) - INLET_VY)  - T1
       L2 = relax_T * (q(i,j,QRHO)*qaux(i,j,QC)*qaux(i,j,QRSPEC)/probhi(2)) * (q(i,j,QTEMP) - INLET_TEMPERATURE) - T2
       
       L3 = relax_U * (qaux(i,j,QC)/probhi(2)) * (q(i,j,QU) - INLET_VX) - ((1.0d0 - beta)*T3)
       L4 = (q(i,j,QV)+qaux(i,j,QC))* (dpdy + (q(i,j,QRHO)*qaux(i,j,QC))*dvdy)
             
     elseif ((bc_type == SlipWall).or.(bc_type == NoSlipWall)) then
     
       ! Values long Y will be computed by mirror functions below
       
     elseif (bc_type .eq. Outflow) then
     
       sigma_out = bc_params(6)
       beta = mach_local_hi_y 

       Kout = sigma_out*(1.0d0 - (mach_local_hi_y**2.0d0))*qaux(i,j,QC)/(probhi(2))
   
       L1 = (Kout*(q(i,j,QPRES) - INLET_PRESSURE)) - ((1.0d0 - beta )*T1)
       L2 = q(i,j,QV) * ( ((qaux(i,j,QC)**2.0d0)*drhody) - dpdy)
       L3 = q(i,j,QV) * dudy
       L4 = (q(i,j,QV)+qaux(i,j,QC))* (dpdy + (q(i,j,QRHO)*qaux(i,j,QC))*dvdy)
        
     else
       call bl_error("Error:: This BC is not yet implemented for hi_y in characteristic form")
     endif
    
     if (q(i,j,QV) == 0.0d0) then
       L1 = L1 / (q(i,j,QV)-qaux(i,j,QC))
       L2 = 0.0d0
       L3 = 0.0d0
       L4 = L4 / (q(i,j,QV)+qaux(i,j,QC))
     else
       L1 = L1 / (q(i,j,QV)-qaux(i,j,QC))
       L2 = L2 / q(i,j,QV)
       L3 = L3 / q(i,j,QV)
       L4 = L4 / (q(i,j,QV)+qaux(i,j,QC))
     endif
      
     ! Compute new spatial derivative
     drhody = (L2 + 0.5d0*(L1 + L4))/(qaux(i,j,QC)**2.0d0)  
     dvdy   = (L4-L1)/(2.0d0*qaux(i,j,QC)*q(i,j,QRHO))
     dudy   = L3
     dpdy   = 0.5d0*(L1+L4)
      
     ! Update ghost cells

     q(i,j+1,QU)    = q(i,j-1,QU) + 2.0d0*dy*dudy
     q(i,j+1,QV)    = q(i,j-1,QV) + 2.0d0*dy*dvdy 
     q(i,j+1,QRHO)  = q(i,j-1,QRHO)  + 2.0d0*dy*drhody
     q(i,j+1,QPRES) = q(i,j-1,QPRES) + 2.0d0*dy*dpdy

     !----------------
     q(i,j+2,QU)    = -2.0d0*q(i,j-1,QU) - 3.0d0*q(i,j,QU) + 6.0d0*q(i,j+1,QU) - 6.0d0*dy*dudy
     q(i,j+2,QV)    = -2.0d0*q(i,j-1,QV) - 3.0d0*q(i,j,QV) + 6.0d0*q(i,j+1,QV) - 6.0d0*dy*dvdy
     q(i,j+2,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 3.0d0*q(i,j,QRHO) + 6.0d0*q(i,j+1,QRHO) - 6.0d0*dy*drhody
     q(i,j+2,QPRES) = -2.0d0*q(i,j-1,QPRES) - 3.0d0*q(i,j,QPRES) + 6.0d0*q(i,j+1,QPRES) - 6.0d0*dy*dpdy

     q(i,j+3,QU)    = 3.0d0*q(i,j-1,QU) +10.0d0*q(i,j,QU) - 18.0d0*q(i,j+1,QU) + 6.0d0*q(i,j+2,QU) + 12.0d0*dy*dudy
     q(i,j+3,QV)    = 3.0d0*q(i,j-1,QV) +10.0d0*q(i,j,QV) - 18.0d0*q(i,j+1,QV) + 6.0d0*q(i,j+2,QV) + 12.0d0*dy*dvdy
     q(i,j+3,QRHO)  = 3.0d0*q(i,j-1,QRHO) +10.0d0*q(i,j,QRHO) - 18.0d0*q(i,j+1,QRHO) + 6.0d0*q(i,j+2,QRHO) + 12.0d0*dy*drhody
     q(i,j+3,QPRES) = 3.0d0*q(i,j-1,QPRES) +10.0d0*q(i,j,QPRES) - 18.0d0*q(i,j+1,QPRES) + 6.0d0*q(i,j+2,QPRES) + 12.0d0*dy*dpdy

     q(i,j+4,QU)    = -2.0d0*q(i,j-1,QU) - 13.0d0*q(i,j,QU) + 24.0d0*q(i,j+1,QU) - 12.0d0*q(i,j+2,QU)  & 
                     + 4.0d0*q(i,j+3,QU) - 12.0d0*dy*dudy 
     q(i,j+4,QV)    = -2.0d0*q(i,j-1,QV) - 13.0d0*q(i,j,QV) + 24.0d0*q(i,j+1,QV) - 12.0d0*q(i,j+2,QV) &
                     + 4.0d0*q(i,j+3,QV) - 12.0d0*dy*dvdy
     q(i,j+4,QRHO)  = -2.0d0*q(i,j-1,QRHO) - 13.0d0*q(i,j,QRHO) + 24.0d0*q(i,j+1,QRHO) - 12.0d0*q(i,j+2,QRHO) &
                     + 4.0d0*q(i,j+3,QRHO) - 12.0d0*dy*drhody
     q(i,j+4,QPRES) = -2.0d0*q(i,j-1,QPRES) - 13.0d0*q(i,j,QPRES) + 24.0d0*q(i,j+1,QPRES) - 12.0d0*q(i,j+2,QPRES) &
                     + 4.0d0*q(i,j+3,QPRES) - 12.0d0*dy*dpdy
                     
                     
     if ((bc_type.eq. NoSlipWall).or.(bc_type .eq. SlipWall)) then
     
       if (bc_type .eq. NoSlipWall) then
         wall_sign = -1.0d0
       else if (bc_type .eq. SlipWall)  then
         wall_sign = 1.0d0
       end if

       q(i,j+1,QU)    = wall_sign*q(i,j,QU)
       q(i,j+2,QU)    = wall_sign*q(i,j-1,QU)
       q(i,j+3,QU)    = wall_sign*q(i,j-2,QU)
       q(i,j+4,QU)    = wall_sign*q(i,j-3,QU)

       q(i,j+1,QV)    = -q(i,j,QV)
       q(i,j+2,QV)    = -q(i,j-1,QV)
       q(i,j+3,QV)    = -q(i,j-2,QV)
       q(i,j+4,QV)    = -q(i,j-3,QV)

       q(i,j+1,QRHO)  = q(i,j,QRHO)
       q(i,j+2,QRHO)  = q(i,j-1,QRHO)
       q(i,j+3,QRHO)  = q(i,j-2,QRHO)
       q(i,j+4,QRHO)  = q(i,j-3,QRHO)

       q(i,j+1,QPRES)  = q(i,j,QPRES)
       q(i,j+2,QPRES)  = q(i,j-1,QPRES)
       q(i,j+3,QPRES)  = q(i,j-2,QPRES)
       q(i,j+4,QPRES)  = q(i,j-3,QPRES)
     
     end if
            
     ! Recompute missing values thanks to EOS
     do hop = domhi(2)+1,domhi(2)+4,1
     
       eos_state % p        = q(i,hop,QPRES )
       eos_state % rho      = q(i,hop,QRHO  )
       eos_state % massfrac = q(i,hop,QFS:QFS+nspec-1)
       eos_state % aux      = q(i,hop,QFX:QFX+naux-1)

       call eos_rp(eos_state)
       q(i,hop,QTEMP)  = eos_state % T
       q(i,hop,QREINT) = eos_state % e * q(i,hop,QRHO)
       q(i,hop,QGAME)  = q(i,hop,QPRES) / q(i,hop,QREINT) + ONE
       
       qaux(i,hop,QDPDR)  = eos_state % dpdr_e
       qaux(i,hop,QDPDE)  = eos_state % dpde
       qaux(i,hop,QGAMC)  = eos_state % gam1
       qaux(i,hop,QC   )  = eos_state % cs
       qaux(i,hop,QCSML)  = max(small, small * qaux(i,hop,QC))

       uin(i,hop,URHO )  = eos_state % rho 
       uin(i,hop,UMX  )  = q(i,hop,QU ) * eos_state % rho 
       uin(i,hop,UMY  )  = q(i,hop,QV ) * eos_state % rho 
       uin(i,hop,UMZ  ) = 0.0d0
       uin(i,hop,UEINT) = eos_state % rho   *  eos_state % e
       uin(i,hop,UEDEN) = eos_state % rho  &
          * (eos_state % e + 0.5d0 * (uin(i,hop,UMX)**2 + uin(i,hop,UMY)**2))
       uin(i,hop,UTEMP) = eos_state % T
       do n=1, nspec
          uin(i,hop,UFS+n-1) = eos_state % rho  *  eos_state % massfrac(n)
       end do 
       
     enddo
   
   enddo
end if



call destroy(eos_state)

end subroutine impose_NSCBC

end module gc_nscbc_mod
