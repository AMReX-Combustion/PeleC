module nbrsTest_nd_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  use amrex_ebcellflag_module, only : get_neighbor_cells
  use pelec_eb_stencil_types_module, only : eb_bndry_geom, eb_bndry_sten, face_sten
  use amrex_constants_module, only: ONE, HALF, TWO, M_PI, FOUR3RD

  implicit none

contains
  pure logical function is_inside (i,j,lo,hi)
    implicit none
    integer, intent(in) :: i,j,lo(2),hi(2)
    is_inside = i.ge.lo(1) .and. i.le.hi(1) &
         .and.  j.ge.lo(2) .and. j.le.hi(2)
  end function is_inside

  subroutine pc_fill_bndry_grad_stencil_amrex(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil_amrex")

      ! Work in process - least squares boundary stencil capability. Currently doesn't work.
      implicit none
      integer,            intent(in   ) :: lo(0:2),hi(0:2)
      integer,            intent(in   ) :: Nebg, Nsten
      type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
      type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)
      real(amrex_real),   intent(in   ) :: dx
      integer :: i, j, k, m, c(dim), s(dim), iv(dim), ivs(dim), sh(dim),  baseiv(dim)
      real(amrex_real) :: cy(-1:1),cz(-1:1), bcs, tsten(-1:1,-1:1,-1:1)
      real(amrex_real), dimension(0:2,0:2,0:2,0:2) :: psten, rsten, sten

      call bl_error('pc_fill_bndry_grad_stencil_amrex has not been implemented in 2D')

    end subroutine pc_fill_bndry_grad_stencil_amrex

    subroutine pc_fill_bndry_grad_stencil_ls(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil_ls")

      ! Work in process - least squares boundary stencil capability. Currently doesn't work.
      implicit none
      integer,            intent(in   ) :: lo(0:2),hi(0:2)
      integer,            intent(in   ) :: Nebg, Nsten
      type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
      type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)
      real(amrex_real),   intent(in   ) :: dx
      integer :: i, j, k, m, c(dim), s(dim), iv(dim), ivs(dim), sh(dim),  baseiv(dim)
      real(amrex_real) :: cy(-1:1),cz(-1:1), bcs, tsten(-1:1,-1:1,-1:1)
      real(amrex_real), dimension(0:2,0:2,0:2,0:2) :: psten, rsten, sten

      call bl_error('pc_fill_bndry_grad_stencil_ls has not been implemented in 2D')

  end subroutine pc_fill_bndry_grad_stencil_ls


  subroutine pc_fill_bndry_grad_stencil(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil")

    implicit none
    integer,            intent(in   ) :: lo(0:1),hi(0:1)
    integer,            intent(in   ) :: Nebg, Nsten
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
    type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)
    real(amrex_real),   intent(in   ) :: dx
    !integer :: i, j, c(dim), s(dim), iv(dim), ivs(dim), sh(dim), L, baseiv(dim)
    !real(amrex_real) :: n(dim), b(dim), x(2), y(2), d(2), fac, sten(0:2,0:2), AREA

    integer :: i, j, m, c(dim), s(dim), iv(dim), ivs(dim), sh(dim), L, baseiv(dim)
    real(amrex_real) :: n(dim), b(dim), x(2), y(2), d(2), fac, sten(-1:1,-1:1), AREA
    real(amrex_real) :: cy(-1:1), bcs, tsten(-1:1,-1:1)
    integer :: ii,jj
    
    AREA = dx**(dim-1)
    fac = AREA / dx

    do L = 0, Nsten-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) ) then

          n = ebg(L)%eb_normal
          ivs = ebg(L)%iv

          c(1) = 1 + INT( MIN(1.d0, ABS(n(2)/n(1)) ) ) ! Either 1 or 2 only
          c(2) = 3 - c(1) ! Either 1 or 2 only

          s = SIGN(1.d0,n(c))
          b = ebg(L)%eb_centroid(c) * s

          baseiv = ivs + SIGN(1.d0,n) - 1 ! From ivs, move to center of stencil, then move to lower-left of that

          x(1) = 1
          x(2) = 2
          y = b(2) + (x - b(1))*ABS(n(c(2))/n(c(1)))

          sh(:) = 0
          if (y(1)<0 .or. y(2)<0) then
             sh(c(2)) = -s(2) ! Slide stencil down to avoid extrapolating, push up eb, shift down base later
             b(2) = b(2) + 1
             y(:) = b(2) + (x(:) - b(1))*ABS(n(c(2))/n(c(1)))
          endif

          d = SQRT( (x-b(1))**2 + (y-b(2))**2)


          sten = 0.d0

          do m=1,2
             cy(-1) = HALF*(y(m)-ONE)*(y(m)-TWO)
             cy( 0) =  -    y(m)     *(y(m)-TWO)
             cy( 1) = HALF* y(m)     *(y(m)-ONE)

             do jj=-1,1
                sten(m-1,jj) = cy(jj)
             enddo
          enddo

          sten(0,-1:1) = sten(0,-1:1) * d(2)/(d(1)*(d(2)-d(1)))
          sten(1,-1:1) = sten(1,-1:1) * d(1)/(d(2)*(d(1)-d(2)))
          bcs = - (d(1)+d(2))/(d(1)*d(2))


          ! Transform stencil into regular stencil structure
          do jj=-1,1
             do ii=-1,1
                iv(c(1)) = (ii+1) * s(1) + ivs(c(1)) - baseiv(c(1)) - 1
                iv(c(2)) = (jj+1) * s(2) + ivs(c(2)) - baseiv(c(2)) - 1
                tsten(iv(1),iv(2)) = sten(ii,jj)
             enddo
          enddo

          grad_stencil(L)%iv = ebg(L)%iv
          grad_stencil(L)%iv_base = baseiv + sh ! Shift base down, if required
          grad_stencil(L)%bcval = fac * ebg(L)%eb_area * bcs
          grad_stencil(L)%val   = fac * ebg(L)%eb_area * tsten 

       endif
    enddo

  end subroutine pc_fill_bndry_grad_stencil

  subroutine pc_apply_eb_boundry_flux_stencil(lo, hi, sten, Nsten, s, slo, shi, D, Dlo, Dhi, &
       bcval, Nvals, bcflux, Nflux, nc) &
       bind(C,name="pc_apply_eb_boundry_flux_stencil")

    implicit none
    integer,          intent(in   ) ::  lo(0:1),  hi(0:1)
    integer,          intent(in   ) :: Nsten, Nvals, Nflux, nc
    type(eb_bndry_sten), intent(in) :: sten(0:Nsten-1)
    real(amrex_real), intent(in   ) :: bcval(0:Nvals-1,1:nc)
    real(amrex_real), intent(inout) :: bcflux(0:Nflux-1,1:nc)
    integer,          intent(in)  :: slo(0:1), shi(0:1)
    integer,          intent(in)  :: Dlo(0:1), Dhi(0:1)
    real(amrex_real), intent(in)  :: s(slo(0):shi(0),slo(1):shi(1),1:nc)
    real(amrex_real), intent(in)  :: D(Dlo(0):Dhi(0),Dlo(1):Dhi(1),1:nc)
    integer :: i,j,L,n,ii,jj

    do L = 0, Nsten-1
       i = sten(L) % iv(0)
       j = sten(L) % iv(1)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) ) then

          ii = sten(L)%iv_base(0)
          jj = sten(L)%iv_base(1)

          do n=1,nc
             bcflux(L,1:nc) = D(i,j,n) * (bcval(L,n) * sten(L)%bcval + &
                  sum(sten(L)%val(-1:1,-1:1) * s(ii:ii+2,jj:jj+2,n)) )
          enddo

       endif
    enddo

  end subroutine pc_apply_eb_boundry_flux_stencil
       
  subroutine pc_apply_eb_boundry_visc_flux_stencil( &
       lo, hi,         &
       sten, Nsten,    &
       ebg,  Nebg,     &
       s,  slo,  shi,  &
       mu, mulo, muhi, &
       xi, xilo, xihi, &
       bcval, Nvals, bcflux, Nflux, nc) &
       bind(C,name="pc_apply_eb_boundry_visc_flux_stencil")

    implicit none
    integer,          intent(in   ) ::  lo(0:1),  hi(0:1)
    integer,          intent(in   ) :: Nsten, Nebg, Nvals, Nflux, nc
    type(eb_bndry_sten), intent(in) :: sten(0:Nsten-1)
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
    real(amrex_real), intent(in   ) :: bcval(0:Nvals-1,1:nc)
    real(amrex_real), intent(inout) :: bcflux(0:Nflux-1,1:nc)
    integer,          intent(in)  :: slo(0:1), shi(0:1)
    integer,          intent(in)  :: mulo(0:1), muhi(0:1)
    integer,          intent(in)  :: xilo(0:1), xihi(0:1)
    real(amrex_real), intent(in)  ::  s( slo(0):shi(0),  slo(1):shi(1) ,1:nc)
    real(amrex_real), intent(in)  :: mu(mulo(0):muhi(0),mulo(1):muhi(1))
    real(amrex_real), intent(in)  :: xi(xilo(0):xihi(0),xilo(1):xihi(1))
    integer :: i,j,L,M,ii,jj,iii,jjj
    real(amrex_real) :: Nmag, n(dim), t(dim)
    real(amrex_real) :: Uo(-1:1,-1:1,dim), Ut(-1:1,-1:1,dim)
    real(amrex_real) :: Qt(dim,dim), dUtdn(dim), tauDotN(dim), bco(dim), bct(dim)

    do L = 0, Nsten-1
       i = sten(L) % iv(0)
       j = sten(L) % iv(1)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) ) then

          ii = sten(L)%iv_base(0)
          jj = sten(L)%iv_base(1)

          Nmag = SQRT(ebg(L)%eb_normal(1)**2 + ebg(L)%eb_normal(2)**2)
          n(1) = ebg(L)%eb_normal(1) / Nmag
          n(2) = ebg(L)%eb_normal(2) / Nmag
          t(1) = -n(2)
          t(2) =  n(1)

          Qt(1,1) = n(1)
          Qt(1,2) = n(2)
          Qt(2,1) = t(1)
          Qt(2,2) = t(2)

          ! Transform velocities at all stencil points to coordinates aligned with EB
          Uo(-1:1,-1:1,1:dim) = s(ii:ii+2,jj:jj+2,1:dim)
          do jjj=-1,1
             do iii=-1,1
                do M=1,dim
                   Ut(iii,jjj,M) = Qt(M,1) * Uo(iii,jjj,1) + Qt(M,2) * Uo(iii,jjj,2)
                enddo
             enddo
          enddo

          ! Transform eb boundary velocities to coordinates aligned with EB
          bco(1:dim) = bcval(L,1:dim)
          bct(1) = Qt(1,1) * bco(1) + Qt(1,2)*bco(2)
          bct(2) = Qt(2,1) * bco(1) + Qt(2,2)*bco(2)

          ! Compute normal derivative (times eb area) using precomputed stencil
          dUtdn(1) = sum(sten(L)%val(-1:1,-1:1) * Ut(-1:1,-1:1,1)) + bct(1) * sten(L)%bcval
          dUtdn(2) = sum(sten(L)%val(-1:1,-1:1) * Ut(-1:1,-1:1,2)) + bct(2) * sten(L)%bcval

          tauDotN(1) = (FOUR3RD*mu(i,j) + xi(i,j)) * dUtdn(1)
          tauDotN(2) =          mu(i,j)            * dUtdn(2)

          bcflux(L,1) = Qt(1,1) * tauDotN(1) + Qt(2,1) * tauDotN(2)
          bcflux(L,2) = Qt(1,2) * tauDotN(1) + Qt(2,2) * tauDotN(2)

       endif
    enddo

  end subroutine pc_apply_eb_boundry_visc_flux_stencil

  subroutine pc_fill_flux_interp_stencil(lo, hi, slo, shi, sten, Nsten, idir, &
       fc, fclo, fchi, fa, falo, fahi) bind(C,name="pc_fill_flux_interp_stencil")

    implicit none
    integer,          intent(in)  ::  lo(0:1),  hi(0:1)
    integer,          intent(in)  :: slo(0:1), shi(0:1)
    integer,          intent(in)  :: Nsten, idir
    type(face_sten),intent(inout) :: sten(0:Nsten-1)
    integer,          intent(in)  :: fclo(0:1), fchi(0:1)
    integer,          intent(in)  :: falo(0:1), fahi(0:1)
    real(amrex_real), intent(in)  :: fc( fclo(0):fchi(0),fclo(1):fchi(1),0:dim-1)
    real(amrex_real), intent(in)  :: fa( falo(0):fahi(0),falo(1):fahi(1))
    integer :: i,j,n,in,jn
    real(amrex_real) :: cx,cy,acx,acy

    if (idir .eq. 0) then
       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)

          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) ) then
             sten(n)%val(:) = 0.d0
             cy = fc(i,j,1)
             jn = SIGN(1.d0, cy)
             acy = ABS(cy)
             sten(n)%val(0) = fa(i,j) * (1.d0 - acy)
             sten(n)%val(jn) = fa(i,j) * acy
          endif
       enddo
    else
       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1)) then
             sten(n)%val(:) = 0.d0
             cx = fc(i,j,0)
             in = SIGN(1.d0, cx)
             acx = ABS(cx)
             sten(n)%val(0) = fa(i,j) * (1.d0 - acx)
             sten(n)%val(in) = fa(i,j) * acx
          endif
       enddo
    endif
  end subroutine pc_fill_flux_interp_stencil

  subroutine pc_apply_face_stencil(lo, hi, slo, shi, sten, Nsten, idir, vin, vin_lo, vin_hi, &
    vout, vout_lo, vout_hi, nc, in_place) bind(C,name="pc_apply_face_stencil")

    implicit none
    integer,          intent(in   ) ::  lo(0:1),  hi(0:1)
    integer,          intent(in   ) :: slo(0:1), shi(0:1)
    integer,          intent(in   ) :: Nsten, idir, nc, in_place
    type(face_sten),  intent(in   ) :: sten(0:Nsten-1)
    integer,          intent(in   ) ::  vin_lo(0:1),  vin_hi(0:1)
    integer,          intent(in   ) :: vout_lo(0:1), vout_hi(0:1)
    real(amrex_real), intent(in   ) ::  vin( vin_lo(0):vin_hi(0),  vin_lo(1):vin_hi(1),  1:nc)
    real(amrex_real), intent(inout) :: vout(vout_lo(0):vout_hi(0),vout_lo(1):vout_hi(1), 1:nc)
    integer :: i,j,L,n
    real(amrex_real), allocatable :: newval(:)

    if (in_place .eq. 1) then

       allocate(newval(0:Nsten-1))

       if (idir.eq.0) then
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   newval(L) = sum(sten(L)%val(-1:1) * vin(i,j-1:j+1,n))
                endif
             enddo

             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   vout(i,j,n) = newval(L)
                endif
             enddo

          enddo
       else
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   newval(L) = sum(sten(L)%val(-1:1) * vin(i-1:i+1,j,n))
                endif
             enddo

             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   vout(i,j,n) = newval(L)
                endif
             enddo

          enddo
       endif
       
       deallocate(newval)

    else

       if (idir.eq.0) then
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   vout(i,j,n) = sum(sten(L)%val(-1:1) * vin(i,j-1:j+1,n))
                endif
             enddo
          enddo
       else
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) ) then
                   vout(i,j,n) = sum(sten(L)%val(-1:1) * vin(i-1:i+1,j,n))
                endif
             enddo
          enddo
       endif

    endif

  end subroutine pc_apply_face_stencil

  subroutine pc_fix_div_and_redistribute( &
       lo, hi,             &
       sv_ebg, Ncut,       &
       flag, fglo, fghi,   &
       f0,   f0lo,   f0hi, &
       f1,   f1lo,   f1hi, &
       ebflux, nebflux,    &
       DC,   DClo,   DChi, &
       W,    Wlo,    Whi,  & 
       vf,   vflo,   vfhi, VOL, nc, &
       as_crse, rr_drho_crse, rdclo, rdchi, rr_flag_crse, rfclo, rfchi, &
       as_fine, dm_as_fine, dflo, dfhi, &
       levmsk, lmlo, lmhi,dt) bind(C,name="pc_fix_div_and_redistribute")
    
    use amrex_eb_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, &
         covered_by_fine=>fine_cell, reredistribution_threshold
    
    use meth_params_module, only: levmsk_notcovered 

    implicit none
    integer,          intent(in   ) ::  lo(0:1),  hi(0:1)
    integer,          intent(in   ) :: nc, Ncut, nebflux
    type(eb_bndry_geom), intent(in   ) :: sv_ebg(0:Ncut-1)
    integer,          intent(in   ) :: fglo(0:1),fghi(0:1)
    integer,          intent(in   ) :: f0lo(0:1),f0hi(0:1)
    integer,          intent(in   ) :: f1lo(0:1),f1hi(0:1)
    integer,          intent(in   ) :: DClo(0:1),DChi(0:1)
    integer,          intent(in   ) :: vflo(0:1),vfhi(0:1)
    integer,          intent(in   ) ::  Wlo(0:1), Whi(0:1)
    integer,          intent(in   ) :: flag(fglo(0):fghi(0),fglo(1):fghi(1))
    real(amrex_real), intent(in   ) ::   f0(f0lo(0):f0hi(0),f0lo(1):f0hi(1),1:nc)
    real(amrex_real), intent(in   ) ::   f1(f1lo(0):f1hi(0),f1lo(1):f1hi(1),1:nc)
    real(amrex_real), intent(inout) ::   ebflux(0:nebflux-1,1:nc)
    real(amrex_real), intent(inout) ::   DC(DClo(0):DChi(0),DClo(1):DChi(1),1:nc)
    real(amrex_real), intent(in   ) ::   vf(vflo(0):vfhi(0),vflo(1):vfhi(1))
    real(amrex_real), intent(in   ) ::    W( Wlo(0):Whi(0) , Wlo(1):Whi(1))
    real(amrex_real), intent(in   ) :: VOL,dt
    real(amrex_real) :: VOLINV, sum_kappa, sum_div, kappa_inv, DNC, sum_kappa_inv
    integer :: i,j,L,n
    integer :: nbr(-1:1,-1:1)
    
    logical, intent(in) :: as_crse, as_fine
    integer, intent(in), dimension(2) :: rdclo,rdchi,rfclo,rfchi,dflo,dfhi,lmlo,lmhi
    real(amrex_real), intent(inout) :: rr_drho_crse(rdclo(1):rdchi(1),rdclo(2):rdchi(2),nc)
    real(amrex_real), intent(out) :: dm_as_fine(dflo(1):dfhi(1),dflo(2):dfhi(2),nc)
    integer,  intent(in) ::  levmsk (lmlo(1):lmhi(1),lmlo(2):lmhi(2))
    integer,  intent(in) ::  rr_flag_crse(rfclo(1):rfchi(1),rfclo(2):rfchi(2))
    real(amrex_real) :: drho, tmp

    integer :: ii,jj,iii,jjj
    logical :: valid_dst_cell
    logical :: as_crse_crse_cell, as_crse_covered_cell, as_fine_valid_cell, as_fine_ghost_cell

    real(amrex_real), allocatable :: dM(:), HD(:)

    allocate(dM(0:Ncut-1))
    allocate(HD(0:Ncut-1))

    dm_as_fine=0.d0
    VOLINV = 1.d0 / VOL

    do n=1,nc

       ! Recompute conservative divergence, DC, on cut cells...need DC in 2 grow cells for final result
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0)-2 .and. i.le.hi(0)+2 &
               .and. j.ge.lo(1)-2 .and. j.le.hi(1)+2 ) then
             kappa_inv = 1.d0 / MAX(vf(i,j),1.d-12)
#ifdef _OPENMP
!$omp atomic read
#endif
             tmp = ebflux(L,n)
#ifdef _OPENMP
!$omp end atomic
#endif
             DC(i,j,n) = - (f0(i+1,j,n) - f0(i,j,n) &
                          + f1(i,j+1,n) - f1(i,j,n) &
                          + tmp) * VOLINV * kappa_inv
          endif
       enddo

       ! Compute non-conservative and hybrid divergence, DNC and HD, and redistribution mass dM in cut cells
       ! We will need this in 1 grow cells (see below), so it depends on having a conservative div in 2 grow cells
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 ) then
             
             call get_neighbor_cells(flag(i,j),nbr)
             sum_kappa = sum(nbr(-1:1,-1:1) * vf(i-1:i+1,j-1:j+1))
             sum_div =   sum(nbr(-1:1,-1:1) * vf(i-1:i+1,j-1:j+1) * DC(i-1:i+1,j-1:j+1,n))
             DNC = sum_div / sum_kappa
             dM(L) = vf(i,j)*(1.d0 - vf(i,j))*(DC(i,j,n) - DNC)
             HD(L) = vf(i,j)*DC(i,j,n) + (1.d0 - vf(i,j))*DNC
          endif
       enddo

       ! Now that we finished computing HD and dM everywhere, it is safe to increment DC to hold HD
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 ) then
             DC(i,j,n) = HD(L)
          endif
       enddo

       ! Redistribute dM - THIS REQUIRES THAT DC BE GOOD IN 1 GROW CELL
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 ) then

             call get_neighbor_cells(flag(i,j),nbr)
             nbr(0,0) = 0.d0 ! redistribute to all neighbors but me
             sum_kappa = sum(nbr(-1:1,-1:1) * vf(i-1:i+1,j-1:j+1) * W(i-1:i+1,j-1:j+1))
             sum_kappa_inv = 1.d0/sum_kappa
             DC(i-1:i+1,j-1:j+1,n) = DC(i-1:i+1,j-1:j+1,n) + dM(L) * nbr(-1:1,-1:1) * W(i-1:i+1,j-1:j+1) * sum_kappa_inv

             !re redistribution book keeping
             as_crse_crse_cell = .false.
             as_crse_covered_cell = .false.
             if (as_crse) then
                  as_crse_crse_cell = is_inside(i,j,lo,hi) .and. &
                           rr_flag_crse(i,j) .eq. crse_fine_boundary_cell
                  as_crse_covered_cell = rr_flag_crse(i,j) .eq. covered_by_fine
             end if

             as_fine_valid_cell = .false.  ! valid cells near box boundary
             as_fine_ghost_cell = .false.  ! ghost cells just outside valid region
             if (as_fine) then
                as_fine_valid_cell = is_inside(i,j,lo,hi)
                as_fine_ghost_cell = levmsk(i,j) .eq. levmsk_notcovered ! not covered by other grids
             end if

             do jj = -1,1
               do ii = -1,1
                  if((ii.ne. 0 .or. jj.ne.0) .and. nbr(ii,jj).eq.1) then

                      iii = i + ii
                      jjj = j + jj

                      drho = dM(L)*sum_kappa_inv*W(iii,jjj)

                      valid_dst_cell = is_inside(iii,jjj,lo,hi)

                      if (as_crse_crse_cell) then
                            if (rr_flag_crse(iii,jjj).eq.covered_by_fine &
                                       .and. vf(i,j).gt.reredistribution_threshold) then
                                rr_drho_crse(i,j,n) = rr_drho_crse(i,j,n) &
                                          + dt*drho*(vf(iii,jjj)/vf(i,j))
                            end if
                      end if

                      if (as_crse_covered_cell) then
                            if (valid_dst_cell) then
                               if (rr_flag_crse(iii,jjj).eq.crse_fine_boundary_cell &
                                          .and. vf(iii,jjj).gt.reredistribution_threshold) then
                                        ! the recipient is a crse/fine boundary cell
                                    rr_drho_crse(iii,jjj,n) = rr_drho_crse(iii,jjj,n) &
                                             - dt*drho
                                end if
                            end if
                     end if

                    if (as_fine_valid_cell) then
                            if (.not.valid_dst_cell) then
                                dm_as_fine(iii,jjj,n) = dm_as_fine(iii,jjj,n) &
                                          + dt*drho*vf(iii,jjj)
                            end if
                    end if

                    if (as_fine_ghost_cell) then
                            if (valid_dst_cell) then
                                     dm_as_fine(i,j,n) = dm_as_fine(i,j,n) &
                                         - dt*drho*vf(iii,jjj)
                            end if
                        end if
                    endif
                enddo
              enddo !end neighbor loop (re redistribution)
          endif
       enddo

    enddo  !component loop

    deallocate(dM)
    deallocate(HD)

  end subroutine pc_fix_div_and_redistribute

  subroutine pc_set_body_state(lo, hi, S, Slo, Shi, mask, mlo, mhi, b, nc, bval) bind(C,name="pc_set_body_state")

    implicit none
    integer,          intent(in   ) :: nc, bval
    integer,          intent(in   ) :: lo(1:2),hi(1:2)
    integer,          intent(in   ) :: Slo(1:2),Shi(1:2)
    integer,          intent(in   ) :: mlo(1:2),mhi(1:2)
    integer,          intent(in   ) :: mask(mlo(1):mhi(1),mlo(2):mhi(2))
    real(amrex_real), intent(inout) :: S(Slo(1):Shi(1),Slo(2):Shi(2),1:nc)
    real(amrex_real), intent(in   ) :: b(1:nc)
    integer :: i,j,n

    do n=1,nc
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             if (mask(i,j).eq.bval) then

#ifdef _OPENMP
!$omp atomic write
#endif
                S(i,j,n)=b(n)
#ifdef _OPENMP
!$omp end atomic
#endif
             endif
          enddo
       enddo
    enddo

  end subroutine pc_set_body_state

  subroutine pc_fill_sv_ebg(lo, hi, ebg, Nebg, vfrac, vflo, vfhi, bcent, blo, bhi, &
       apx, axlo, axhi, apy, aylo, ayhi) &
       bind(C,name="pc_fill_sv_ebg")

    implicit none
    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer, dimension(3), intent(in) :: vflo, vfhi, blo, bhi, axlo, axhi, aylo, ayhi
    integer,          intent(in   ) :: Nebg
    type(eb_bndry_geom), intent(inout) :: ebg(0:Nebg-1)
    integer :: i,j,L,n,ii,jj


    real(amrex_real), intent(in) :: vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    real(amrex_real), intent(in) :: bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(amrex_real), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))

    real(amrex_real) :: axm, axp, aym, ayp, apnorm, apnorminv
    do L = 0, Nebg-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) ) then

          axm = apx(i,j,axlo(3))
          axp = apx(i+1,j,axlo(3))
          aym = apy(i,j,aylo(3))
          ayp = apy(i,j+1,aylo(3))

          apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
          apnorminv = -1.d0 / apnorm
          ebg(L) % eb_normal(1) = (axm-axp) * apnorminv  ! pointing to the wall
          ebg(L) % eb_normal(2) = (aym-ayp) * apnorminv

          ebg(L) % eb_area = apnorm

          ebg(L) % eb_centroid(1) = bcent(i,j,blo(3),1)
          ebg(L) % eb_centroid(2) = bcent(i,j,blo(3),2)
       endif
    enddo

  end subroutine pc_fill_sv_ebg

end module nbrsTest_nd_module
