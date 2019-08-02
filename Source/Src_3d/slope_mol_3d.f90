module slope_module
  
  implicit none

  private

  public slopex, slopey, slopez

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine slopex(q,flatn,qd_lo,qd_hi, &
                        dqxal,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv, &
                        qaux, qa_lo, qa_hi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqxal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)

      integer i, j, k, n,nsp

      real(amrex_real), contiguous, pointer :: dlft(:,:), drgt(:,:)
      double precision :: slop, dq1
  
      
      real(amrex_real), contiguous, pointer :: dsgn(:)
      real(amrex_real), contiguous, pointer :: dlim(:,:,:,:),df(:,:,:,:),dcen(:,:,:,:)


      call bl_allocate (dlft, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (drgt, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (dsgn, ilo1-2,ihi1+2)
      call bl_allocate (dcen, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
!     call bl_allocate (df, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
      call bl_allocate (dlim, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqxal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else

!ONLY X SLOPE CODED HERE

            ! Compute slopes in first coordinate direction
         do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1

!              do i = ilo1-2, ihi1+2
               do i = ilo1-1, ihi1+1

!  du, dp, dv, dw in first four slots
                  dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                  dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                  dlft(i,3) = q(i,j,k,QV) - q(i-1,j,k,QV) 
                  dlft(i,4) = q(i,j,k,QW) - q(i-1,j,k,QW) 

                  drgt(i,1) = 0.5d0*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                  drgt(i,2) = 0.5d0*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                  drgt(i,3) = q(i+1,j,k,QV) - q(i,j,k,QV) 
                  drgt(i,4) = q(i+1,j,k,QW) - q(i,j,k,QW) 

                  do nsp = 1, nspecies
                      dlft(i,4+nsp) = q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1)-  q(i-1,j,k,QRHO)*q(i-1,j,k,QFS+nsp-1) -  &
                                q(i ,j,k,QFS+nsp-1)* (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/qaux(i,j,k,QC)**2
                      drgt(i,4+nsp) = q(i+1,j,k,QRHO)*q(i+1,j,k,QFS+nsp-1)-  q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1) -  &
                                q(i ,j,k,QFS+nsp-1)* (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/qaux(i,j,k,QC)**2
                  enddo
               enddo
 



               do n=1,4+nspecies

               ! First compute Fromm slopes
!                 do i = ilo1-2, ihi1+2
                  do i = ilo1-1, ihi1+1
                     dcen(i,j,k,n) = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn(i) = sign(ONE, dcen(i,j,k,n))
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim(i,j,k,n) = slop
                     else
                        dlim(i,j,k,n) = ZERO
                     endif
!                    df(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                     dqxal(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                  enddo
   
                  ! Now compute limited fourth order slopes
!                 do i = ilo1-1, ihi1+1
!                    dsgn(i) = sign(ONE, dcen(i,j,k,n))
!                    dq1       = FOUR3RD*dcen(i,j,k,n) - SIXTH*(df(i+1,j,k,n) + df(i-1,j,k,n))
!                    dqxal(i,j,k,n) = flatn(i,j,k)*dsgn(i)*min(dlim(i,j,k,n),abs(dq1))
!                 enddo
               enddo

            enddo
         enddo
      end if ! If plm_iorder .eq. 1

      call bl_deallocate (dlft)
      call bl_deallocate (drgt)
      call bl_deallocate (dsgn)
      call bl_deallocate (dlim)
 !    call bl_deallocate (  df)
      call bl_deallocate (dcen)

      end subroutine slopex

      subroutine slopey(q,flatn,qd_lo,qd_hi, &
                        dqyal,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv, &
                        qaux, qa_lo, qa_hi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)      
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqyal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nv)

      integer i, j, k, n

      double precision, pointer :: dlft(:,:), drgt(:,:)
      double precision :: slop, dq1
  
      
      real(amrex_real), contiguous, pointer::dsgn(:),dlim(:,:,:,:),df(:,:,:,:),dcen(:,:,:,:)

      call bl_allocate (dlft, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (drgt, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (dsgn, ilo1-2,ihi1+2)
!     call bl_allocate (df, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
      call bl_allocate (dcen, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
      call bl_allocate (dlim, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)


      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqyal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else


            ! Compute slopes in first coordinate direction
         do k = ilo3-1, ihi3+1
!           do j = ilo2-2, ihi2+2
            do j = ilo2-1, ihi2+1

               do i = ilo1-1, ihi1+1

                     dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,3) = q(i,j,k,QU) - q(i,j-1,k,QU) 
                     dlft(i,4) = q(i,j,k,QW) - q(i,j-1,k,QW) 
                    
                     do nsp = 1, nspecies
                        dlft(i,4+nsp) = q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1)-  q(i-1,j,k,QRHO)*q(i,j-1,k,QFS+nsp-1) -  &
                                q(i ,j,k,QFS+nsp-1)* (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/qaux(i,j,k,QC)**2
                     enddo

                     drgt(i,1) = 0.5d0*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,2) = 0.5d0*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,3) = q(i,j+1,k,QU) - q(i,j,k,QU) 
                     drgt(i,4) = q(i,j+1,k,QW) - q(i,j,k,QW) 

                     do nsp = 1, nspecies
                        drgt(i,4+nsp) = q(i,j+1,k,QRHO)*q(i,j+1,k,QFS+nsp-1) - q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1) -  &
                           q(i ,j,k,QFS+nsp-1)* (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/qaux(i,j,k,QC)**2
                     enddo

               do n=1,4+nspecies

                  ! First compute Fromm slopes
                  do i = ilo1-1, ihi1+1
                     dcen(i,j,k,n) = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn(i) = sign(ONE, dcen(i,j,k,n))
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim(i,j,k,n) = slop
                     else
                        dlim(i,j,k,n) = ZERO
                     endif
!                    df(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                     dqzal(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                  enddo

               enddo
            enddo
         enddo



               ! Now compute limited fourth order slopes
!         do n = 1,nv
!           do k = ilo3-1, ihi3+1
!              do j = ilo2-1, ihi2+1
!                 do i = ilo1-1, ihi1+1
!                    dsgn(i) = sign(ONE, dcen(i,j,k,n))
!                    dq1       = FOUR3RD*dcen(i,j,k,n) - SIXTH*(df(i,j+1,k,n) + df(i,j-1,k,n))
!                    dqyal(i,j,k,n) = flatn(i,j,k)*dsgn(i)*min(dlim(i,j,k,n),abs(dq1))
!                 enddo
!              enddo
!           enddo
!        enddo
      endif

      call bl_deallocate (dlft)
      call bl_deallocate (drgt)
      call bl_deallocate (dsgn)
      call bl_deallocate (dlim)
!     call bl_deallocate (  df)
      call bl_deallocate (dcen)

      end subroutine slopey

      subroutine slopez(q,flatn,qd_lo,qd_hi, &
                        dqzal,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv, &
                        qaux, qa_lo, qa_hi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv

      double precision, intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision, intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqzal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)

      integer i, j, k, n

      double precision, pointer :: dlft(:,:), drgt(:,:)
      double precision :: slop, dq1   
      
      double precision, pointer:: dsgn(:),dlim(:,:,:,:),df(:,:,:,:),dcen(:,:,:,:)

      call bl_allocate (dlft, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (drgt, ilo1-2,ihi1+2,1,nv)
      call bl_allocate (dsgn, ilo1-2,ihi1+2)
!     call bl_allocate (df, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
      call bl_allocate (dcen, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)
      call bl_allocate (dlim, ilo1-2,ihi1+2,ilo2-2,ihi2+2,ilo3-2,ihi3+2,1,nv)


      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqzal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else


            ! Compute slopes in first coordinate direction
!        do k = ilo3-2, ihi3+2
         do k = ilo3-1, ihi3+1
            do j = ilo2-1, ihi2+1

               do i = ilo1-1, ihi1+1

                  dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                  dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                  dlft(i,3) = q(i,j,k,QU) - q(i,j,k-1,QU) 
                  dlft(i,4) = q(i,j,k,QV) - q(i,j,k-1,QV) 

                  do nsp = 1, nspecies
                     dlft(i,4+nsp) = q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1)-  q(i,j,k-1,QRHO)*q(i,j,k-1,QFS+nsp-1) -  &
                             q(i ,j,k,QFS+nsp-1)* (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/qaux(i,j,k,QC)**2
                  enddo

                  drgt(i,1) = 0.5d0*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                  drgt(i,2) = 0.5d0*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                  drgt(i,3) = q(i,j,k+1,QU) - q(i,j,k,QU) 
                  drgt(i,3) = q(i,j,k+1,QV) - q(i,j,k,QV)

                  do nsp = 1, nspecies
                     drgt(i,4+nsp) = q(i,j,k+1,QRHO)*q(i,j,k+1,QFS+nsp-1) - q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1) -  &
                        q(i ,j,k,QFS+nsp-1)* (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/qaux(i,j,k,QC)**2
                  enddo

               enddo

               do n=1,4+nspecies

                  ! First compute Fromm slopes
                  do i = ilo1-1, ihi1+1
                     dcen(i,j,k,n) = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn(i) = sign(ONE, dcen(i,j,k,n))
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim(i,j,k,n) = slop
                     else
                        dlim(i,j,k,n) = ZERO
                     endif
!                    df(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                     dqzal(i,j,k,n) = dsgn(i)*min( dlim(i,j,k,n), abs(dcen(i,j,k,n)) )
                  enddo

               enddo
            enddo
         enddo

               ! Now compute limited fourth order slopes
!        do n = 1,nv
!           do k = ilo3-1, ihi3+1
!              do j = ilo2-1, ihi2+1
!                 do i = ilo1-1, ihi1+1
!                    dsgn(i) = sign(ONE, dcen(i,j,k,n))
!                    dq1       = FOUR3RD*dcen(i,j,k,n) - SIXTH*(df(i,j,k+1,n) + df(i,j,k-1,n))
!                    dqzal(i,j,k,n) = flatn(i,j,k)*dsgn(i)*min(dlim(i,j,k,n),abs(dq1))
!                 enddo
!              enddo
!           enddo
!        enddo

      endif ! End if plm_iorder .eq. 1

      call bl_deallocate (dlft)
      call bl_deallocate (drgt)
      call bl_deallocate (dsgn)
      call bl_deallocate (dlim)
!     call bl_deallocate (  df)
      call bl_deallocate (dcen)

      end subroutine slopez

end module slope_module
