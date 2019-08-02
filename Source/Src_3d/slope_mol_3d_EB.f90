module slope_module
  
  use amrex_ebcellflag_module, only : get_neighbor_cells, is_covered_cell

  implicit none

  private

  public slopex, slopey, slopez

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine slopex(q,flatn,qd_lo,qd_hi, &
                        dqxal,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv,nva,&
                        domlo,domhi,&
                        qaux, qa_lo, qa_hi, &
                        flag, fglo, fghi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies
      use prob_params_module, only : physbc_lo, physbc_hi, Inflow

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv, nva
      integer, intent(in) :: domlo(3),domhi(3)

      integer, intent(in) :: fglo(3),fghi(3)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nva)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqxal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)

      integer i, j, k, n, nsp

      real(amrex_real) :: dlft(ilo1:ihi1,1:nv), drgt(ilo1:ihi1,1:nv)
      real(amrex_real) :: slop, dsgn, dlim, dcen

      integer :: nbr(-1:1,-1:1,-1:1)

      logical :: flagArrayL(ilo1:ihi1)
      logical :: flagArrayR(ilo1:ihi1)
      
      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3, ihi3
               do j = ilo2, ihi2
                  do i = ilo1, ihi1
                     dqxal(i,j,k,n) = ZERO
                  enddo
               enddo
            enddo
         enddo

      else

         do k = ilo3, ihi3
            do j = ilo2, ihi2

               do n=1,nv
                  do i = ilo1, ihi1
                     dlft(i,n) = 0.d0
                     drgt(i,n) = 0.d0
                  enddo
               enddo

               do i = ilo1, ihi1
                  call get_neighbor_cells( flag(i,j,k), nbr )
                  flagArrayL(i) = nbr(-1,0,0).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
                  flagArrayR(i) = nbr(+1,0,0).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
               enddo

               do i = ilo1, ihi1
                  if (flagArrayL(i)) then
                     dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                     dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU))
                     dlft(i,3) = q(i,j,k,QV) - q(i-1,j,k,QV)
                     dlft(i,4) = q(i,j,k,QW) - q(i-1,j,k,QW)
                  endif
                  if (flagArrayR(i)) then
                     drgt(i,1) = 0.5d0*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                     drgt(i,2) = 0.5d0*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU))
                     drgt(i,3) = q(i+1,j,k,QV) - q(i,j,k,QV)
                     drgt(i,4) = q(i+1,j,k,QW) - q(i,j,k,QW)
                  endif
               enddo

               do nsp = 1, nspecies
                  do i = ilo1, ihi1
                     if (flagArrayL(i)) then
                        dlft(i,4+nsp) = q(i  ,j,k,QRHO)*q(i,j,k,QFS+nsp-1)-q(i-1,j,k,QRHO)*q(i-1,j,k,QFS+nsp-1) &
                             -    q(i,j,k,QFS+nsp-1)*(q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/qaux(i,j,k,QC)**2
                     endif
                     if (flagArrayR(i)) then
                        drgt(i,4+nsp) = q(i+1,j,k,QRHO)*q(i+1,j,k,QFS+nsp-1)-q(i,j,k,QRHO)*q(i ,j,k,QFS+nsp-1) -  &
                           q(i,j,k,QFS+nsp-1)* (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/qaux(i,j,k,QC)**2
                     endif
                  enddo
               enddo

               do n=1, nv
                  do i = ilo1, ihi1
                     dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn = sign(ONE, dcen)
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim = slop
                     else
                        dlim = ZERO
                     endif
                     dqxal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                  enddo
               enddo
            enddo
         enddo
      endif

    end subroutine slopex

    subroutine slopey(q,flatn,qd_lo,qd_hi, &
         dqyal,qpd_lo,qpd_hi, &
         ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv,nva,&
         domlo,domhi, &
         qaux, qa_lo, qa_hi, &
         flag, fglo, fghi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies
      use prob_params_module, only : physbc_lo, physbc_hi, Inflow

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)      
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv, nva
      integer, intent(in) :: domlo(3),domhi(3)

      integer, intent(in) :: fglo(3),fghi(3)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqyal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nva)

      integer i, j, k, n, nsp

      real(amrex_real) :: dlft(ilo1:ihi1,1:nv), drgt(ilo1:ihi1,1:nv)
      double precision :: slop, dsgn, dlim, dcen

      integer :: nbr(-1:1,-1:1,-1:1)

      logical :: flagArrayL(ilo1:ihi1)
      logical :: flagArrayR(ilo1:ihi1)

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3, ihi3
               do j = ilo2, ihi2
                  do i = ilo1, ihi1
                     dqyal(i,j,k,n) = ZERO
                  enddo
               enddo
            enddo
         enddo

      else
         
         do k = ilo3, ihi3
            do j = ilo2, ihi2

               do n=1,nv
                  do i = ilo1, ihi1
                     dlft(i,n) = 0.d0
                     drgt(i,n) = 0.d0
                  enddo
               enddo
               
               do i = ilo1, ihi1
                  call get_neighbor_cells( flag(i,j,k), nbr )
                  flagArrayL(i) = nbr(0,-1,0).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
                  flagArrayR(i) = nbr(0,+1,0).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
               enddo

               do i = ilo1, ihi1
                  if (flagArrayL(i)) then
                     dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV))
                     dlft(i,3) = q(i,j,k,QU) - q(i,j-1,k,QU)
                     dlft(i,4) = q(i,j,k,QW) - q(i,j-1,k,QW)
                  endif
                  if (flagArrayR(i)) then
                     drgt(i,1) = 0.5d0*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,2) = 0.5d0*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV))
                     drgt(i,3) = q(i,j+1,k,QU) - q(i,j,k,QU)
                     drgt(i,4) = q(i,j+1,k,QW) - q(i,j,k,QW)
                  endif

               enddo

               do nsp = 1, nspecies
                  do i = ilo1, ihi1
                     if (flagArrayL(i)) then
                        dlft(i,4+nsp) = q(i,j  ,k,QRHO)*q(i,j  ,k,QFS+nsp-1)-q(i,j-1,k,QRHO)*q(i,j-1,k,QFS+nsp-1) &
                             -    q(i,j,k,QFS+nsp-1)*(q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/qaux(i,j,k,QC)**2
                     endif
                     if (flagArrayR(i)) then
                        drgt(i,4+nsp) = q(i,j+1,k,QRHO)*q(i,j+1,k,QFS+nsp-1)-q(i,j,k,QRHO)*q(i,j,k,QFS+nsp-1) &
                             -    q(i,j,k,QFS+nsp-1)*(q(i,j+1,k,QPRES) - q(i,j  ,k,QPRES))/qaux(i,j,k,QC)**2
                     endif
                  enddo
               enddo

               do n=1, nv
                  do i = ilo1, ihi1
                     dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn = sign(ONE, dcen)
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim = slop
                     else
                        dlim = ZERO
                     endif
                     dqyal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                  enddo
               enddo
            enddo
          enddo
      endif

    end subroutine slopey

    subroutine slopez(q,flatn,qd_lo,qd_hi, &
         dqzal,qpd_lo,qpd_hi, &
         ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,nv,nva, & 
         domlo,domhi, &
         qaux, qa_lo, qa_hi, &
         flag, fglo, fghi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies
      use prob_params_module, only : physbc_lo, physbc_hi, Inflow

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: qa_lo(3), qa_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, ilo3, ihi3, nv, nva
      integer, intent(in) :: domlo(3),domhi(3)

      integer, intent(in) :: fglo(3),fghi(3)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      double precision :: dqzal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),nva)

      integer i, j, k, n, nsp

      real(amrex_real) :: dlft(ilo1:ihi1,1:nv), drgt(ilo1:ihi1,1:nv)
      double precision :: slop, dsgn, dlim, dcen
      
      integer :: nbr(-1:1,-1:1,-1:1)

      logical :: flagArrayL(ilo1:ihi1)
      logical :: flagArrayR(ilo1:ihi1)


      if(plm_iorder.eq.1) then

         do n = 1, nv
            do k = ilo3, ihi3
               do j = ilo2, ihi2
                  do i = ilo1, ihi1
                     dqzal(i,j,k,n) = ZERO
               enddo
            enddo
            enddo
         enddo

      else

         do k = ilo3, ihi3
            do j = ilo2, ihi2

               do n=1,nv
                  do i = ilo1, ihi1
                     dlft(i,n) = 0.d0
                     drgt(i,n) = 0.d0
                  enddo
               enddo
               
               do i = ilo1, ihi1
                  call get_neighbor_cells( flag(i,j,k), nbr )
                  flagArrayL(i) = nbr(0,0,-1).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
                  flagArrayR(i) = nbr(0,0,+1).eq.1 .and. .not. is_covered_cell(flag(i,j,k))
               enddo

               do i = ilo1, ihi1
                  if (flagArrayL(i)) then
                     dlft(i,1) = 0.5d0*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,2) = 0.5d0*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW))
                     dlft(i,3) = q(i,j,k,QU) - q(i,j,k-1,QU)
                     dlft(i,4) = q(i,j,k,QV) - q(i,j,k-1,QV)
                  endif
                  if (flagArrayR(i)) then
                     drgt(i,1) = 0.5d0*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          - 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,2) = 0.5d0*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/qaux(i,j,k,QC) &
                          + 0.5d0*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW))
                     drgt(i,3) = q(i,j,k+1,QU) - q(i,j,k,QU)
                     drgt(i,4) = q(i,j,k+1,QV) - q(i,j,k,QV)
                  endif
               enddo

               do nsp = 1, nspecies
                  do i = ilo1, ihi1
                     if (flagArrayL(i)) then
                        dlft(i,4+nsp) = q(i,j,k  ,QRHO)*q(i,j,k,  QFS+nsp-1)-q(i,j,k-1,QRHO)*q(i,j,k-1,QFS+nsp-1) &
                             -    q(i,j,k,QFS+nsp-1)*(q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/qaux(i,j,k,QC)**2
                     endif
                     if (flagArrayR(i)) then
                        drgt(i,4+nsp) = q(i,j,k+1,QRHO)*q(i,j,k+1,QFS+nsp-1)-q(i,j,k  ,QRHO)*q(i,j,k  ,QFS+nsp-1) &
                             -    q(i,j,k,QFS+nsp-1)*(q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/qaux(i,j,k,QC)**2
                     endif
                  enddo
               enddo

               do n=1, nv
                  do i = ilo1, ihi1
                     dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                     dsgn = sign(ONE, dcen)
                     slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                     if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                        dlim = slop
                     else
                        dlim = ZERO
                     endif
                     dqzal(i,j,k,n) = dsgn*min( dlim, abs(dcen) )
                  enddo
               enddo
            enddo
         enddo
      endif

    end subroutine slopez

end module slope_module
