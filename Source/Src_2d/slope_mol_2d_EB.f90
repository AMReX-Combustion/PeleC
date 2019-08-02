module slope_module
  
  use amrex_ebcellflag_module, only : get_neighbor_cells, is_covered_cell

  implicit none

  private

  public slopex, slopey

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine slopex(q,flatn,qd_lo,qd_hi, &
                        dqxal,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ihi1,ihi2,nv,nva,&
                        domlo,domhi,&
                        qaux, qa_lo, qa_hi, &
                        flag, fglo, fghi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies

      implicit none

      integer          :: qd_lo(2), qd_hi(2)
      integer          :: qpd_lo(2),qpd_hi(2)
      integer          :: qa_lo(2), qa_hi(2)
      integer          :: ilo1, ilo2, ihi1, ihi2, nv, nva
      integer, intent(in) :: domlo(2), domhi(2)

      integer, intent(in) :: fglo(2),fghi(2)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2))

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),nva)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2))
      double precision :: dqxal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),nv)

      integer i, j, n, nsp

      real(amrex_real) :: dlft(ilo1:ihi1,1:nv), drgt(ilo1:ihi1,1:nv)
      double precision :: slop, dsgn, dlim, dcen

      integer :: nbr(-1:1,-1:1)

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do j = ilo2, ihi2
               do i = ilo1, ihi1
                  dqxal(i,j,n) = ZERO
               enddo
            enddo
         enddo

      else

         ! Compute slopes in first coordinate direction
         do j = ilo2, ihi2
            do i = ilo1, ihi1

               call get_neighbor_cells(flag(i,j), nbr)

               if (nbr(-1,0).eq.1 .and. .not. is_covered_cell(flag(i,j)) ) then
                  dlft(i,1) = 0.5d0*(q(i,j,QPRES)-q(i-1,j,QPRES))/qaux(i,j,QC) - 0.5d0*q(i,j,QRHO)*(q(i,j,QU) - q(i-1,j,QU))
                  dlft(i,2) = 0.5d0*(q(i,j,QPRES)-q(i-1,j,QPRES))/qaux(i,j,QC) + 0.5d0*q(i,j,QRHO)*(q(i,j,QU) - q(i-1,j,QU))
                  dlft(i,3) = q(i,j,QV) - q(i-1,j,QV)
                  dlft(i,4) = 0.d0

                  do nsp = 1, nspecies
                     dlft(i,4+nsp) = q(i  ,j,QRHO)*q(i,j,QFS+nsp-1)-q(i-1,j,QRHO)*q(i-1,j,QFS+nsp-1) &
                          -    q(i,j,QFS+nsp-1)*(q(i  ,j,QPRES) - q(i-1,j,QPRES))/qaux(i,j,QC)**2
                  enddo

                  dlft(i,4+nspecies+1:nv) = 0.0
               else
                  dlft(i,:) = 0.d0
               endif

               if (nbr(1,0).eq.1 .and. .not. is_covered_cell(flag(i,j)) ) then
                  drgt(i,1) = 0.5d0*(q(i+1,j,QPRES)-q(i,j,QPRES))/qaux(i,j,QC) - 0.5d0*q(i,j,QRHO)*(q(i+1,j,QU) - q(i,j,QU))
                  drgt(i,2) = 0.5d0*(q(i+1,j,QPRES)-q(i,j,QPRES))/qaux(i,j,QC) + 0.5d0*q(i,j,QRHO)*(q(i+1,j,QU) - q(i,j,QU))
                  drgt(i,3) = q(i+1,j,QV) - q(i,j,QV)
                  drgt(i,4) = 0.d0

                  do nsp = 1, nspecies
                     drgt(i,4+nsp) = q(i+1,j,QRHO)*q(i+1,j,QFS+nsp-1)-q(i,j,QRHO)*q(i  ,j,QFS+nsp-1) &
                          -    q(i,j,QFS+nsp-1)*(q(i+1,j,QPRES) - q(i  ,j,QPRES))/qaux(i,j,QC)**2
                  enddo

                  drgt(i,4+nspecies+1:nv) = 0.0
               else
                  drgt(i,:) = 0.d0
               endif

            enddo

            do n=1, nv

               ! Compute Fromm slopes
               do i = ilo1, ihi1
                  dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                  dsgn = sign(ONE, dcen)
                  slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                  if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                     dlim = slop
                  else
                     dlim = ZERO
                  endif
                  dqxal(i,j,n) = dsgn*min( dlim, abs(dcen) )
               enddo
   
            enddo
         enddo
      endif ! If plm_iorder .eq. 1

    end subroutine slopex

    subroutine slopey(q,flatn,qd_lo,qd_hi, &
                      dqyal,qpd_lo,qpd_hi, &
                      ilo1,ilo2,ihi1,ihi2,nv,nva,&
                      domlo,domhi, &
                      qaux, qa_lo, qa_hi, &
                      flag, fglo, fghi)

      use amrex_fort_module, only : amrex_real
      use amrex_mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use amrex_constants_module
      use network, only : nspecies

      implicit none

      integer          :: qd_lo(2), qd_hi(2)
      integer          :: qpd_lo(2),qpd_hi(2)
      integer          :: qa_lo(2), qa_hi(2)      
      integer          :: ilo1, ilo2, ihi1, ihi2, nv, nva
      integer, intent(in) :: domlo(2),domhi(2)

      integer, intent(in) :: fglo(2),fghi(2)
      integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2))

      double precision :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),nv)
      double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2))
      double precision :: dqyal(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),nv)
      double precision :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),nva)

      integer i, j, n, nsp

      real(amrex_real) :: dlft(ilo1:ihi1,1:nv), drgt(ilo1:ihi1,1:nv)
      double precision :: slop, dsgn, dlim, dcen

      integer ::   nbr(-1:1,-1:1)

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do j = ilo2, ihi2
               do i = ilo1, ihi1
                  dqyal(i,j,n) = ZERO
               enddo
            enddo
         enddo

      else

         ! Compute slopes in second coordinate direction
         do j = ilo2, ihi2
            do i = ilo1, ihi1

               call get_neighbor_cells(flag(i,j), nbr)

               if (nbr(0,-1).eq.1 .and. .not. is_covered_cell(flag(i,j)) ) then
                  dlft(i,1) = 0.5d0*(q(i,j,QPRES)-q(i,j-1,QPRES))/qaux(i,j,QC) - 0.5d0*q(i,j,QRHO)*(q(i,j,QV) - q(i,j-1,QV))
                  dlft(i,2) = 0.5d0*(q(i,j,QPRES)-q(i,j-1,QPRES))/qaux(i,j,QC) + 0.5d0*q(i,j,QRHO)*(q(i,j,QV) - q(i,j-1,QV))
                  dlft(i,3) = q(i,j,QU) - q(i,j-1,QU)
                  dlft(i,4) = 0.d0

                  do nsp = 1, nspecies
                     dlft(i,4+nsp) = q(i,j  ,QRHO)*q(i,j  ,QFS+nsp-1)-q(i,j-1,QRHO)*q(i,j-1,QFS+nsp-1) &
                          -    q(i,j,QFS+nsp-1)*(q(i,j  ,QPRES) - q(i,j-1,QPRES))/qaux(i,j,QC)**2
                  enddo

                  dlft(i,4+nspecies+1:nv) = 0.0
               else
                  dlft(i,:) = 0.d0
               endif

               if (nbr(0,1).eq.1 .and. .not. is_covered_cell(flag(i,j)) ) then
                  drgt(i,1) = 0.5d0*(q(i,j+1,QPRES)-q(i,j,QPRES))/qaux(i,j,QC) - 0.5d0*q(i,j,QRHO)*(q(i,j+1,QV) - q(i,j,QV))
                  drgt(i,2) = 0.5d0*(q(i,j+1,QPRES)-q(i,j,QPRES))/qaux(i,j,QC) + 0.5d0*q(i,j,QRHO)*(q(i,j+1,QV) - q(i,j,QV))
                  drgt(i,3) = q(i,j+1,QU) - q(i,j,QU)
                  drgt(i,4) = 0.d0

                  do nsp = 1, nspecies
                     drgt(i,4+nsp) = q(i,j+1,QRHO)*q(i,j+1,QFS+nsp-1)-q(i,j,QRHO)*q(i,j,QFS+nsp-1) &
                          -    q(i,j,QFS+nsp-1)*(q(i,j+1,QPRES) - q(i,j  ,QPRES))/qaux(i,j,QC)**2
                  enddo

                  drgt(i,4+nspecies+1:nv) = 0.0
               else
                  drgt(i,:) = 0.d0
               endif
            enddo

            do n=1, nv
               ! Compute Fromm slopes
               do i = ilo1, ihi1
                  dcen = 0.5d0 * (dlft(i,n)+drgt(i,n))
                  dsgn = sign(ONE, dcen)
                  slop = 2.d0* min( abs(dlft(i,n)), abs(drgt(i,n)) )
                  if (dlft(i,n)*drgt(i,n) .ge. ZERO) then
                     dlim = slop
                  else
                     dlim = ZERO
                  endif
                  dqyal(i,j,n) = dsgn*min( dlim, abs(dcen) )
               enddo
            enddo
         enddo
      endif  ! If plm_iorder .eq. 1

    end subroutine slopey

end module slope_module
