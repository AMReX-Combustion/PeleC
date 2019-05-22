module nbrsTest_nd_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_ebcellflag_module, only : get_neighbor_cells
  use pelec_eb_stencil_types_module, only : eb_bndry_geom, eb_bndry_sten, face_sten
  use amrex_constants_module, only: ONE, HALF, TWO, FOUR3RD

  implicit none

contains

  pure logical function is_inside (i,j,k,lo,hi)
    integer, intent(in) :: i,j,k,lo(3),hi(3)
    is_inside = i.ge.lo(1) .and. i.le.hi(1) &
         .and.  j.ge.lo(2) .and. j.le.hi(2) &
         .and.  k.ge.lo(3) .and. k.le.hi(3)
  end function is_inside

  pure subroutine mysort(c, n)
    real(amrex_real), intent(in) :: n(dim)
    integer, intent(out) :: c(dim)
    logical          :: mask(dim)
    integer :: m,imax(1)
    mask = .true.
    do m=1,dim
       imax = MAXLOC(ABS(n),mask)
       c(m) = imax(1)
       mask(imax(1)) = .false.
    enddo
  end subroutine mysort

  subroutine pc_fill_bndry_grad_stencil_amrex(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil_amrex")

    ! This one fills the stencil using the strategy in amrex_mlebabeclap_grad routine
    ! compute_dphidn_3d
    ! Currently work in progress (doesn't work)

    ! Notes from WZ 4/24/19 follow:

    ! What it does is
    !(1) find a point in the boundary normal direction that is `dx_eb` away from the boundary centroid.
    !(2) find out 8 neighboring cell centers and use them to do linear interpolation.
    !(3) compute `dphidn = (phi_interp - phi_b)/dx_eb`.
    ! Here `dx_eb =  max(0.3d0, (kappa*kappa-0.25d0)/(2.d0*kappa))`, where kappa is volume fraction.

    !For 1d case with `kappa > 0.9*` ( cannot remember the exact number), it is equivalent to use
    ! `phi_b` and two cell centers to construct a polynomial fit and obtain dphidn on the boundary.
    ! For kappa = 1, this is exactly same as what we do at domain wall in non-EB solver (with max_order = 3).
    ! The problem of the stencil in PeleC and EB/CNS seems to be that in the case of kappa close to 1,
    !    the cut cell itself is excluded.
    ! That seems to cause problems for flows with steep gradient.
    ! This is probably the reason in one of the ANAG paper the approach based on least square fit was
    ! developed and it was mentioned the original stencil was not stable for low Reynolds flows.
    ! End notes from WZ

    use amrex_mlebabeclap_3d_module, only : amrex_get_dx_eb

    ! Array bounds
    integer,            intent(in   ) :: lo(0:2), hi(0:2)
    integer,            intent(in   ) :: Nebg, Nsten

    ! EB Data
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)

    ! Stencil data to fill
    type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)

    ! Grid spacing
    real(amrex_real),   intent(in   ) :: dx

    ! Local variables
    real(amrex_real) :: fac, AREA
    real(amrex_real) :: dx_eb, vf, dg
    real(amrex_real) :: bctx, bcty, bctz ! Boundary centroids
    real(amrex_real) :: anrmx, anrmy, anrmz ! Normals
    real(amrex_real) :: sten(-1:1, -1:1, -1:1) ! Raw stencil
    real(amrex_real) :: gx, gy, gz, sx, sy, sz
    real(amrex_real) :: gxy, gxz, gyz, gxyz
    real(amrex_real) :: anrm
    real(amrex_real) :: sten_sum
    integer :: ii, jj, kk ! Offsets for stencil

    ! Cell indices
    integer :: i, j, k, L
    AREA = dx**(dim-1)
    fac = AREA / dx


    do L = 0, Nsten-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       k = ebg(L) % iv(2)

       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          dx_eb = amrex_get_dx_eb(ebg(L)%eb_vfrac)

          bctx = ebg(L) % eb_centroid(1)
          bcty = ebg(L) % eb_centroid(2)
          bctz = ebg(L) % eb_centroid(3)

          anrmx = ebg(L) % eb_normal(1)
          anrmy = ebg(L) % eb_normal(2)
          anrmz = ebg(L) % eb_normal(3)

          anrm = sqrt( anrmx*anrmx + anrmy*anrmy + anrmz*anrmz)
          anrmx = anrmx/anrm
          anrmy = anrmy/anrm
          anrmz = anrmz/anrm

          dg = dx_eb / max(abs(anrmx),abs(anrmy),abs(anrmz))

          ! Renormalize normal
          gx = bctx - dg*anrmx
          gy = bcty - dg*anrmy
          gz = bctz - dg*anrmz
          sx =  sign(one,anrmx)
          sy =  sign(one,anrmy)
          sz =  sign(one,anrmz)
          ii =  - int(sx)
          jj =  - int(sy)
          kk =  - int(sz)

          gx = sx*gx
          gy = sy*gy
          gz = sz*gz
          gxy = gx*gy
          gxz = gx*gz
          gyz = gy*gz
          gxyz = gx*gy*gz

          sten = 0.0d0
          sten(0,0,0) = (one+gx+gy+gz+gxy+gxz+gyz+gxyz)
          sten(0,0,kk) = (-gz - gxz - gyz - gxyz)
          sten(0,jj,0) = (-gy - gxy - gyz - gxyz)
          sten(0,jj,kk) = (gyz + gxyz)
          sten(ii,0,0) = (-gx - gxy - gxz - gxyz)
          sten(ii,0,kk) = (gxz + gxyz)
          sten(ii,jj,0) = (gxy + gxyz)
          sten(ii,jj,kk) = (-gxyz)

          grad_stencil(L) % iv = ebg(L) % iv
          grad_stencil(L) % iv_base = grad_stencil(L) % iv - 1
          grad_stencil(L) % bcval = one/dg *fac *ebg(L)%eb_area

          grad_stencil(L) % val(-1:1,-1:1,-1:1) = -one/dg*sten(-1:1,-1:1,-1:1) * fac *ebg(L)%eb_area

          sten_sum = sum(sten)
          if (abs(sten_sum - one) .gt. 1.0e-10) then
             write(*,*) 'Trouble: stencil does not add up to unity!'
             write(*,*) 'L =', L, ' normals = ', anrmx, anrmy, anrmz
             write(*,*) 'L =', L, ' centroid = ', bctx, bcty, bctz
             write(*,*) 'L =', L, ' dg = ', dg, ' sten_sum:', sten_sum
             write(*,'(27(E10.4,2x))') sten
             call amrex_abort("pc_fill_bndry_grad_stencil_amrex trouble with stencil")
          endif

       endif ! Restrict to box
    end do ! Loop over cut cells
  end subroutine pc_fill_bndry_grad_stencil_amrex

    subroutine pc_fill_bndry_grad_stencil_ls(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil_ls")

      ! Work in process - least squares boundary stencil capability. Currently doesn't work.

    integer,            intent(in   ) :: lo(0:2),hi(0:2)
    integer,            intent(in   ) :: Nebg, Nsten
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
    type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)
    real(amrex_real),   intent(in   ) :: dx
    integer :: i, j, k, m, c(dim), s(dim), iv(dim), ivs(dim), sh(dim),  baseiv(dim)
    real(amrex_real) :: cy(-1:1),cz(-1:1), bcs, tsten(-1:1,-1:1,-1:1)
    real(amrex_real), dimension(0:2,0:2,0:2,0:2) :: psten, rsten, sten

    ! Local variables
    real(amrex_real) :: r11, r11sq, r12, r22, r22sq, r13, r23, r33, r33sq, beta, alph(0:2)

    real(amrex_real) :: anrmx, anrmy, anrmz ! Normals
    integer :: ii,jj,kk
    integer :: L
    integer :: nls
    integer :: inls
    integer :: si, sj, sk ! Stencil offsets

    do L = 0, Nsten-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       k = ebg(L) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          anrmx = ebg(L)%eb_normal(1)
          anrmy = ebg(L)%eb_normal(2)
          anrmz = ebg(L)%eb_normal(3)

          do sk = 0, 2
             do sj = 0, 2
                do si = 0, 2
                   rsten(si,sj,sk,0) = -si*sign(1.d0,anrmx) - ebg(L)%eb_centroid(1)
                   rsten(si,sj,sk,1) = -sj*sign(1.d0,anrmy) - ebg(L)%eb_centroid(2)
                   rsten(si,sj,sk,2) = -sk*sign(1.d0,anrmz) - ebg(L)%eb_centroid(3)
                enddo
             enddo
          enddo

          r11sq = (sum(rsten(0:2,0:2,0:2,0)*rsten(0:2,0:2,0:2,0)))
          r11 = sqrt(r11)
          r12 = (sum(rsten(0:2,0:2,0:2,0)*rsten(0:2,0:2,0:2,1)))/r11
          r22 = (sum(rsten(0:2,0:2,0:2,1)*rsten(0:2,0:2,0:2,1))) - r12*r12
          r22sq = r22*r22
          r13 = (sum(rsten(0:2,0:2,0:2,0)*rsten(0:2,0:2,0:2,2)))/r11
          r23 = (sum(rsten(0:2,0:2,0:2,1)*rsten(0:2,0:2,0:2,2))  - &
                 sum(rsten(0:2,0:2,0:2,0)*rsten(0:2,0:2,0:2,2))*r12/r11 )/r22
          r33 = sqrt(sum(rsten(0:2,0:2,0:2,2)*rsten(0:2,0:2,0:2,2)) - (r13*r12 + r23*r23))
          r33sq = r33*r33
          beta = (r12*r23 - r13*r22)/(r11*r22)

          do sk = 0, 2
             do sj = 0, 2
                do si = 0, 2
                   alph(0) = rsten(si,sj,sk,0) / r11sq
                   alph(1) = (rsten(si,sj,sk,1) - r12/r11*rsten(si,sj,sk,0))/r22sq
                   alph(2) = (rsten(si,sj,sk,2) - r23/r22*rsten(si,sj,sk,1) + beta*rsten(si,sj,sk,0))/r33sq
                   sten(si,sj,sk,0) = alph(0) - r12/r11*alph(1) + beta*alph(2)
                   sten(si,sj,sk,1) = alph(1) - r23/r22*alph(2)
                   sten(si,sj,sk,2) = alph(2)
                enddo
             enddo
          enddo

          ! Now, grad phi = sum(sten*phi-phi_bc). We want the component normal to the cut face
          sten(:,:,:,0) = sten(:,:,:,0)*anrmx + sten(:,:,:,1)*anrmy + sten(:,:,:,2)*anrmz

          grad_stencil(L)%iv = ebg(L)%iv
          grad_stencil(L)%iv_base(0) =ebg(L)% iv(0) - sign(1.d0,anrmx) - 1 ! move to center of stencil, then move to lower left of that
          grad_stencil(L)%iv_base(1) =ebg(L)% iv(1) - sign(1.d0,anrmy) - 1 ! move to center of stencil, then move to lower left of that
          grad_stencil(L)%iv_base(2) =ebg(L)% iv(2) - sign(1.d0,anrmz) - 1 ! move to center of stencil, then move to lower left of that
          grad_stencil(L)%bcval = -1.0*sum(sten(:,:,:,0))
          grad_stencil(L)%val = sten(:,:,:,0) ! TODO: check this doesn't need to be normalized by sum(sten(:,:,:,0))

       endif ! inside box to work on
    enddo ! Loop over stencils
  end subroutine pc_fill_bndry_grad_stencil_ls


  subroutine pc_fill_bndry_grad_stencil(lo, hi, ebg, Nebg, grad_stencil, Nsten, dx) &
       bind(C,name="pc_fill_bndry_grad_stencil")

    integer,            intent(in   ) :: lo(0:2),hi(0:2)
    integer,            intent(in   ) :: Nebg, Nsten
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
    type(eb_bndry_sten),intent(inout) :: grad_stencil(0:Nsten-1)
    real(amrex_real),   intent(in   ) :: dx
    integer :: i, j, k, m, c(dim), s(dim), iv(dim), ivs(dim), sh(dim), L, baseiv(dim)
    real(amrex_real) :: n(dim), b(dim), x(2), y(2), z(2), d(2), fac, sten(-1:1,-1:1,-1:1), AREA
    real(amrex_real) :: cy(-1:1),cz(-1:1), bcs, tsten(-1:1,-1:1,-1:1)

    integer :: ii,jj,kk

    AREA = dx**(dim-1)
    fac = AREA / dx

    do L = 0, Nsten-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       k = ebg(L) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          n = ebg(L)%eb_normal
          call mysort(c,n)
          ivs = ebg(L)%iv
          s = SIGN(1.d0,n(c))
          b = ebg(L)%eb_centroid(c) * s

          baseiv = ivs + SIGN(1.d0,n) - 1 ! From ivs, move to center of stencil, then move to lower-left of that

          x(1) = 1
          x(2) = 2
          y = b(2) + (x - b(1))*ABS(n(c(2))/n(c(1)))
          z = b(3) + (x - b(1))*ABS(n(c(3))/n(c(1)))

          sh(:) = 0
          if (y(1)<0 .or. y(2)<0) then
             sh(c(2)) = -s(2) ! Slide stencil down to avoid extrapolating, push up eb, shift down base later
             b(2) = b(2) + 1
             y(:) = b(2) + (x(:) - b(1))*ABS(n(c(2))/n(c(1)))
          endif
          if (z(1)<0 .or. z(2)<0) then
             sh(c(3)) = -s(3) ! Slide stencil down to avoid extrapolating, push up eb, shift down base later
             b(3) = b(3) + 1
             z(:) = b(3) + (x(:) - b(1))*ABS(n(c(3))/n(c(1)))
          endif

          d = SQRT( (x-b(1))**2 + (y-b(2))**2 + (z-b(3))**2)

          sten = 0

          ! The two intersections, that are d(1) and d(2) away from the eb centroid, are both in 
          !  y-z planes, bounded in (0:2)x(0:2) in normalized coordinates
          ! For point m, we interpolate z=0,1,2 lines, to (y(m),0), (y(m),1) and (y(m),2), and
          !  then interpolate along y=y(m) to (y(m),z(m))
          do m=1,2
             cy(-1) = HALF*(y(m)-ONE)*(y(m)-TWO)
             cy( 0) =      -y(m)     *(y(m)-TWO)
             cy( 1) = HALF* y(m)     *(y(m)-ONE)

             cz(-1) = HALF*(z(m)-ONE)*(z(m)-TWO)
             cz( 0) =      -z(m)     *(z(m)-TWO)
             cz( 1) = HALF* z(m)     *(z(m)-ONE)

             do kk=-1,1
                do jj=-1,1
                   sten(m-1,jj,kk) = cy(jj)*cz(kk)
                enddo
             enddo
          enddo

          sten(0,-1:1,-1:1) = sten(0,-1:1,-1:1) * d(2)/(d(1)*(d(2)-d(1)))
          sten(1,-1:1,-1:1) = sten(1,-1:1,-1:1) * d(1)/(d(2)*(d(1)-d(2)))
          bcs = - (d(1)+d(2))/(d(1)*d(2))

          ! Transform stencil into regular stencil structure
          do kk=-1,1
             do jj=-1,1
                do ii=-1,1
                   iv(c(1)) = (ii+1) * s(1) + ivs(c(1)) - baseiv(c(1)) - 1
                   iv(c(2)) = (jj+1) * s(2) + ivs(c(2)) - baseiv(c(2)) - 1
                   iv(c(3)) = (kk+1) * s(3) + ivs(c(3)) - baseiv(c(3)) - 1
                   tsten(iv(1),iv(2),iv(3)) = sten(ii,jj,kk)
                enddo
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

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: Nsten, Nvals, Nflux, nc
    type(eb_bndry_sten), intent(in) :: sten(0:Nsten-1)
    real(amrex_real), intent(in   ) :: bcval(0:Nvals-1,1:nc)
    real(amrex_real), intent(inout) :: bcflux(0:Nflux-1,1:nc)
    integer,          intent(in)  :: slo(0:2), shi(0:2)
    integer,          intent(in)  :: Dlo(0:2), Dhi(0:2)
    real(amrex_real), intent(in)  :: s(slo(0):shi(0),slo(1):shi(1),slo(2):shi(2),1:nc)
    real(amrex_real), intent(in)  :: D(Dlo(0):Dhi(0),Dlo(1):Dhi(1),Dlo(2):Dhi(2),1:nc)
    integer :: i,j,k,L,n,ii,jj,kk


    do L = 0, Nsten-1
       i = sten(L) % iv(0)
       j = sten(L) % iv(1)
       k = sten(L) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          ii = sten(L)%iv_base(0)
          jj = sten(L)%iv_base(1)
          kk = sten(L)%iv_base(2)


          do n=1,nc
             bcflux(L,n) = D(i,j,k,n) * (bcval(L,n) * sten(L)%bcval + &
                  sum(sten(L)%val(-1:1,-1:1,-1:1) * s(ii:ii+2,jj:jj+2,kk:kk+2,n)) )
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

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: Nsten, Nebg, Nvals, Nflux, nc
    type(eb_bndry_sten), intent(in) :: sten(0:Nsten-1)
    type(eb_bndry_geom),intent(in   ) :: ebg(0:Nebg-1)
    real(amrex_real), intent(in   ) :: bcval(0:Nvals-1,1:nc)
    real(amrex_real), intent(inout) :: bcflux(0:Nflux-1,1:nc)
    integer,          intent(in)  ::  slo(0:2),  shi(0:2)
    integer,          intent(in)  :: mulo(0:2), muhi(0:2)
    integer,          intent(in)  :: xilo(0:2), xihi(0:2)
    real(amrex_real), intent(in)  ::  s( slo(0):shi(0),  slo(1):shi(1),  slo(2):shi(2) ,1:nc)
    real(amrex_real), intent(in)  :: mu(mulo(0):muhi(0),mulo(1):muhi(1),mulo(2):muhi(2))
    real(amrex_real), intent(in)  :: xi(xilo(0):xihi(0),xilo(1):xihi(1),xilo(2):xihi(2))
    integer :: i,j,k,L,M,ii,jj,kk,iii,jjj,kkk
    real(amrex_real) :: Nmag, n(dim), t1(dim), t2(dim), denom, ndota
    real(amrex_real) :: Uo(-1:1,-1:1,-1:1,dim), Ut(-1:1,-1:1,-1:1,dim), alpha(dim)
    real(amrex_real) :: Qt(dim,dim), dUtdn(dim), tauDotN(dim), bco(dim), bct(dim)

    do L = 0, Nsten-1
       i = sten(L) % iv(0)
       j = sten(L) % iv(1)
       k = sten(L) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          ii = sten(L)%iv_base(0)
          jj = sten(L)%iv_base(1)
          kk = sten(L)%iv_base(2)

          Nmag = SQRT(ebg(L)%eb_normal(1)**2 + ebg(L)%eb_normal(2)**2 + ebg(L)%eb_normal(3)**2)
          n(1) = ebg(L)%eb_normal(1) / Nmag
          n(2) = ebg(L)%eb_normal(2) / Nmag
          n(3) = ebg(L)%eb_normal(3) / Nmag

          alpha = 0.d0
          alpha(MINLOC(ABS(n))) = 1.d0

          ndota = n(1)*alpha(1) + n(2)*alpha(2) + n(3)*alpha(3)
          t1(1) = alpha(1) - ndota*n(1)
          t1(2) = alpha(2) - ndota*n(2)
          t1(3) = alpha(3) - ndota*n(3)
          denom = 1.d0 / SQRT(t1(1)**2 + t1(2)**2 + t1(3)**2)
          t1(1) = t1(1) * denom
          t1(2) = t1(2) * denom
          t1(3) = t1(3) * denom

          t2(1) = n(2)*t1(3) - n(3)*t1(2)
          t2(2) = n(3)*t1(1) - n(1)*t1(3)
          t2(3) = n(1)*t1(2) - n(2)*t1(1)

          Qt(1,1) = n(1)
          Qt(1,2) = n(2)
          Qt(1,3) = n(3)
          Qt(2,1) = t1(1)
          Qt(2,2) = t1(2)
          Qt(2,3) = t1(3)
          Qt(3,1) = t2(1)
          Qt(3,2) = t2(2)
          Qt(3,3) = t2(3)

          ! Transform velocities at all stencil points to coordinates aligned with EB
          Uo(-1:1,-1:1,-1:1,1:dim) = s(ii:ii+2,jj:jj+2,kk:kk+2,1:dim)
          do kkk=-1,1
             do jjj=-1,1
                do iii=-1,1
                   do M=1,dim
                      Ut(iii,jjj,kkk,M) = Qt(M,1) * Uo(iii,jjj,kkk,1) &
                           +              Qt(M,2) * Uo(iii,jjj,kkk,2) &
                           +              Qt(M,3) * Uo(iii,jjj,kkk,3)
                   enddo
                enddo
             enddo
          enddo

          ! Transform eb boundary velocities to coordinates aligned with EB
          bco(1:dim) = bcval(L,1:dim)
          bct(1) = Qt(1,1) * bco(1) + Qt(1,2)*bco(2) + Qt(1,3)*bco(3)
          bct(2) = Qt(2,1) * bco(1) + Qt(2,2)*bco(2) + Qt(2,3)*bco(3)
          bct(3) = Qt(3,1) * bco(1) + Qt(3,2)*bco(2) + Qt(3,3)*bco(3)

          ! Compute normal derivative (times eb area) using precomputed stencil
          dUtdn(1) = sum(sten(L)%val(-1:1,-1:1,-1:1) * Ut(-1:1,-1:1,-1:1,1)) + bct(1) * sten(L)%bcval
          dUtdn(2) = sum(sten(L)%val(-1:1,-1:1,-1:1) * Ut(-1:1,-1:1,-1:1,2)) + bct(2) * sten(L)%bcval
          dUtdn(3) = sum(sten(L)%val(-1:1,-1:1,-1:1) * Ut(-1:1,-1:1,-1:1,3)) + bct(3) * sten(L)%bcval

          tauDotN(1) = (FOUR3RD*mu(i,j,k) + xi(i,j,k)) * dUtdn(1)
          tauDotN(2) =          mu(i,j,k)              * dUtdn(2)
          tauDotN(3) =          mu(i,j,k)              * dUtdn(3)

          bcflux(L,1) = Qt(1,1) * tauDotN(1) + Qt(2,1) * tauDotN(2) + Qt(3,1) * tauDotN(3)
          bcflux(L,2) = Qt(1,2) * tauDotN(1) + Qt(2,2) * tauDotN(2) + Qt(3,2) * tauDotN(3)
          bcflux(L,3) = Qt(1,3) * tauDotN(1) + Qt(2,3) * tauDotN(2) + Qt(3,3) * tauDotN(3)
       endif
    enddo

  end subroutine pc_apply_eb_boundry_visc_flux_stencil

  subroutine pc_fill_flux_interp_stencil(lo, hi, slo, shi, sten, Nsten, idir, &
       fc, fclo, fchi, fa, falo, fahi) bind(C,name="pc_fill_flux_interp_stencil")

    integer,          intent(in)  ::  lo(0:2),  hi(0:2)
    integer,          intent(in)  :: slo(0:2), shi(0:2)
    integer,          intent(in)  :: Nsten, idir
    type(face_sten),intent(inout) :: sten(0:Nsten-1)
    integer,          intent(in)  :: fclo(0:2), fchi(0:2)
    integer,          intent(in)  :: falo(0:2), fahi(0:2)
    real(amrex_real), intent(in)  :: fc( fclo(0):fchi(0),fclo(1):fchi(1),fclo(2):fchi(2),2)
    real(amrex_real), intent(in)  :: fa( falo(0):fahi(0),falo(1):fahi(1),falo(2):fahi(2))
    integer :: i,j,k,n,in,jn,kn
    real(amrex_real) :: cx,cy,cz,acx,acy,acz

    if (idir .eq. 0) then
       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             sten(n)%val = 0.d0
             cy = fc(i,j,k,1)
             cz = fc(i,j,k,2)
             jn = SIGN(1.d0, cy)
             kn = SIGN(1.d0, cz)
             acy = ABS(cy)
             acz = ABS(cz)
             sten(n)%val(0, 0 ) = fa(i,j,k) * (1.d0 - acy) * (1.d0 - acz)
             sten(n)%val(jn,0 ) = fa(i,j,k) *      acy     * (1.d0 - acz)
             sten(n)%val(0, kn) = fa(i,j,k) * (1.d0 - acy) *      acz
             sten(n)%val(jn,kn) = fa(i,j,k) *      acy     *      acz
          endif
       enddo
    else if (idir.eq.1) then
       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             sten(n)%val = 0.d0
             cx = fc(i,j,k,1)
             cz = fc(i,j,k,2)
             in = SIGN(1.d0, cx)
             kn = SIGN(1.d0, cz)
             acx = ABS(cx)
             acz = ABS(cz)
             sten(n)%val(0, 0 ) = fa(i,j,k) * (1.d0 - acx) * (1.d0 - acz)
             sten(n)%val(in,0 ) = fa(i,j,k) *      acx     * (1.d0 - acz)
             sten(n)%val(0, kn) = fa(i,j,k) * (1.d0 - acx) *      acz
             sten(n)%val(in,kn) = fa(i,j,k) *      acx     *      acz
          endif
       enddo
    else
       do n = 0, Nsten-1
          i = sten(n) % iv(0)
          j = sten(n) % iv(1)
          k = sten(n) % iv(2)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) &
               .and. k.ge.lo(2) .and. k.le.hi(2) ) then
             sten(n)%val = 0.d0
             cx = fc(i,j,k,1)
             cy = fc(i,j,k,2)
             in = SIGN(1.d0, cx)
             jn = SIGN(1.d0, cy)
             acx = ABS(cx)
             acy = ABS(cy)
             sten(n)%val(0, 0 ) = fa(i,j,k) * (1.d0 - acx) * (1.d0 - acy)
             sten(n)%val(in,0 ) = fa(i,j,k) *      acx     * (1.d0 - acy)
             sten(n)%val(0, jn) = fa(i,j,k) * (1.d0 - acx) *      acy
             sten(n)%val(in,jn) = fa(i,j,k) *      acx     *      acy
          endif
       enddo
    endif
  end subroutine pc_fill_flux_interp_stencil

  subroutine pc_apply_face_stencil(lo, hi, slo, shi, sten, Nsten, idir, vin, vin_lo, vin_hi, &
    vout, vout_lo, vout_hi, nc, in_place) bind(C,name="pc_apply_face_stencil")

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: slo(0:2), shi(0:2)
    integer,          intent(in   ) :: Nsten, idir, nc, in_place
    type(face_sten),  intent(in   ) :: sten(0:Nsten-1)
    integer,          intent(in   ) ::  vin_lo(0:2),  vin_hi(0:2)
    integer,          intent(in   ) :: vout_lo(0:2), vout_hi(0:2)
    real(amrex_real), intent(in   ) ::  vin( vin_lo(0):vin_hi(0),  vin_lo(1):vin_hi(1),  vin_lo(2):vin_hi(2),  1:nc)
    real(amrex_real), intent(inout) :: vout(vout_lo(0):vout_hi(0),vout_lo(1):vout_hi(1),vout_lo(2):vout_hi(2), 1:nc)
    integer :: i,j,k,L,n
    real(amrex_real), allocatable :: newval(:)

    if (in_place .eq. 1) then

       allocate(newval(0:Nsten-1))

       if (idir.eq.0) then
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   newval(L) = sum(sten(L)%val(-1:1,-1:1) * vin(i,j-1:j+1,k-1:k+1,n))
                endif
             enddo

             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = newval(L)
                endif
             enddo

          enddo
       else if (idir.eq.1) then
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   newval(L) = sum(sten(L)%val(-1:1,-1:1) * vin(i-1:i+1,j,k-1:k+1,n))
                endif
             enddo

             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = newval(L)
                endif
             enddo

          enddo
       else
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   newval(L) = sum(sten(L)%val(-1:1,-1:1) * vin(i-1:i+1,j-1:j+1,k,n))
                endif
             enddo

             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = newval(L)
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
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = sum(sten(L)%val(-1:1,-1:1) * vin(i,j-1:j+1,k-1:k+1,n))
                endif
             enddo
          enddo
       else if (idir.eq.1) then
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = sum(sten(L)%val(-1:1,-1:1) * vin(i-1:i+1,j,k-1:k+1,n))
                endif
             enddo
          enddo
       else
          do n = 1,nc
             do L = 0, Nsten-1
                i = sten(L) % iv(0)
                j = sten(L) % iv(1)
                k = sten(L) % iv(2)
                if (i.ge.lo(0) .and. i.le.hi(0) &
                     .and. j.ge.lo(1) .and. j.le.hi(1) &
                     .and. k.ge.lo(2) .and. k.le.hi(2) ) then
                   vout(i,j,k,n) = sum(sten(L)%val(-1:1,-1:1) * vin(i-1:i+1,j-1:j+1,k,n))
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
       f2,   f2lo,   f2hi, &
       ebflux, nebflux,    &
       DC,   DClo,   DChi, &
       W,    Wlo,    Whi,  &
       vf,   vflo,   vfhi, VOL, nc, &
       as_crse, rr_drho_crse, rdclo, rdchi, rr_flag_crse, rfclo, rfchi, &
       as_fine, dm_as_fine, dflo, dfhi, &
       levmsk, lmlo, lmhi,dt) bind(C,name="pc_fix_div_and_redistribute")

    use amrex_eb_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, &
         covered_by_fine=>fine_cell, reredistribution_threshold

    use meth_params_module, only: levmsk_notcovered, eb_small_vfrac

    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer,          intent(in   ) :: nc, Ncut, nebflux
    type(eb_bndry_geom), intent(in   ) :: sv_ebg(0:Ncut-1)
    integer,          intent(in   ) :: fglo(0:2),fghi(0:2)
    integer,          intent(in   ) :: f0lo(0:2),f0hi(0:2)
    integer,          intent(in   ) :: f1lo(0:2),f1hi(0:2)
    integer,          intent(in   ) :: f2lo(0:2),f2hi(0:2)
    integer,          intent(in   ) :: DClo(0:2),DChi(0:2)
    integer,          intent(in   ) :: vflo(0:2),vfhi(0:2)
    integer,          intent(in   ) ::  Wlo(0:2), Whi(0:2)
    integer,          intent(in   ) :: flag(fglo(0):fghi(0),fglo(1):fghi(1),fglo(2):fghi(2))
    real(amrex_real), intent(in   ) ::   f0(f0lo(0):f0hi(0),f0lo(1):f0hi(1),f0lo(2):f0hi(2),1:nc)
    real(amrex_real), intent(in   ) ::   f1(f1lo(0):f1hi(0),f1lo(1):f1hi(1),f1lo(2):f1hi(2),1:nc)
    real(amrex_real), intent(in   ) ::   f2(f2lo(0):f2hi(0),f2lo(1):f2hi(1),f2lo(2):f2hi(2),1:nc)
    real(amrex_real), intent(inout) ::   ebflux(0:nebflux-1,1:nc)
    real(amrex_real), intent(inout) ::   DC(DClo(0):DChi(0),DClo(1):DChi(1),DClo(2):DChi(2),1:nc)
    real(amrex_real), intent(in   ) ::   vf(vflo(0):vfhi(0),vflo(1):vfhi(1),vflo(2):vfhi(2))
    real(amrex_real), intent(in   ) ::    W( Wlo(0):Whi(0) , Wlo(1):Whi(1) , Wlo(2):Whi(2))
    real(amrex_real), intent(in   ) :: VOL,dt
    real(amrex_real) :: VOLINV, sum_kappa, sum_div, kappa_inv, DNC
    integer :: i,j,k,L,n
    integer :: nbr(-1:1,-1:1,-1:1)
    real(amrex_real) :: dM(0:Ncut-1), HD(0:Ncut-1), sum_kappa_inv
    
    integer, intent(in) :: as_crse, as_fine
    integer, intent(in), dimension(3) :: rdclo,rdchi,rfclo,rfchi,dflo,dfhi,lmlo,lmhi
    real(amrex_real), intent(inout) :: rr_drho_crse(rdclo(1):rdchi(1),rdclo(2):rdchi(2),rdclo(3):rdchi(3),nc)
    real(amrex_real), intent(out) :: dm_as_fine(dflo(1):dfhi(1),dflo(2):dfhi(2),dflo(3):dfhi(3),nc)
    integer,  intent(in) ::  levmsk (lmlo(1):lmhi(1),lmlo(2):lmhi(2),lmlo(3):lmhi(3))
    integer,  intent(in) ::  rr_flag_crse(rfclo(1):rfchi(1),rfclo(2):rfchi(2),rfclo(3):rfchi(3))
    real(amrex_real) :: drho

    integer :: ii,jj,kk,iii,jjj,kkk
    logical :: valid_dst_cell
    logical :: as_crse_crse_cell, as_crse_covered_cell, as_fine_valid_cell, as_fine_ghost_cell

    dm_as_fine=0.d0
    VOLINV = 1.d0 / VOL

    do n=1,nc

       ! Recompute conservative divergence, DC, on cut cells...need DC in 2 grow cells for final result
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (       i.ge.lo(0)-2 .and. i.le.hi(0)+2 &
               .and. j.ge.lo(1)-2 .and. j.le.hi(1)+2 &
               .and. k.ge.lo(2)-2 .and. k.le.hi(2)+2 ) then
             kappa_inv = 1.d0 / MAX(vf(i,j,k),1.d-12)
             DC(i,j,k,n) = - ( f0(i+1,j,k,n) - f0(i,j,k,n) &
                  +            f1(i,j+1,k,n) - f1(i,j,k,n) &
                  +            f2(i,j,k+1,n) - f2(i,j,k,n) + ebflux(L,n)) * VOLINV * kappa_inv
          endif
       enddo

       ! Compute non-conservative and hybrid divergence, DNC and HD, and redistribution mass dM in cut cells
       ! Will need in 1 grow cells (see below), so it depends on having a conservative div in 2 grow cells
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (       i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 &
               .and. k.ge.lo(2)-1 .and. k.le.hi(2)+1 ) then
             call get_neighbor_cells(flag(i,j,k),nbr)
             sum_kappa = sum(nbr(-1:1,-1:1,-1:1) * vf(i-1:i+1,j-1:j+1,k-1:k+1))
             sum_div =   sum(nbr(-1:1,-1:1,-1:1) * vf(i-1:i+1,j-1:j+1,k-1:k+1) * DC(i-1:i+1,j-1:j+1,k-1:k+1,n))
             DNC = sum_div / sum_kappa
             if (sv_ebg(L) % eb_vfrac < eb_small_vfrac) then
                 dM(L) = vf(i,j,k)*(DC(i,j,k,n))
                 HD(L) = 0.0d0
             else
                 dM(L) = vf(i,j,k)*(1.d0 - vf(i,j,k))*(DC(i,j,k,n) - DNC)
                 HD(L) = vf(i,j,k)*DC(i,j,k,n) + (1.d0 - vf(i,j,k))*DNC
             endif

          endif
       enddo

       ! Now that we finished computing HD and dM everywhere, it is safe to increment DC to hold HD
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (       i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 &
               .and. k.ge.lo(2)-1 .and. k.le.hi(2)+1 ) then
             DC(i,j,k,n) = HD(L)
          endif
       enddo

       ! Redistribute dM - THIS REQUIRES THAT DC BE GOOD IN 1 GROW CELL
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (       i.ge.lo(0)-1 .and. i.le.hi(0)+1 &
               .and. j.ge.lo(1)-1 .and. j.le.hi(1)+1 &
               .and. k.ge.lo(2)-1 .and. k.le.hi(2)+1 ) then

             call get_neighbor_cells(flag(i,j,k),nbr)
             nbr(0,0,0) = 0.d0 ! redistribute to all neighbors but me and those that are really small (for which we've elsewhere adjusted HD)
             do kk=-1,1
                do jj=-1,1
                   do ii=-1,1
                      if(vf(ii+i,jj+j,kk+k) .lt. eb_small_vfrac) nbr(ii,jj,kk) = 0.0d0
                   enddo
                enddo
             enddo

             sum_kappa = sum(nbr(-1:1,-1:1,-1:1) * vf(i-1:i+1,j-1:j+1,k-1:k+1) * W(i-1:i+1,j-1:j+1,k-1:k+1))
             sum_kappa_inv = 1.d0 / sum_kappa
             DC(i-1:i+1,j-1:j+1,k-1:k+1,n) = DC(i-1:i+1,j-1:j+1,k-1:k+1,n) &
                  + dM(L) * nbr(-1:1,-1:1,-1:1) * W(i-1:i+1,j-1:j+1,k-1:k+1) * sum_kappa_inv
        
          !re redistribution book keeping
          as_crse_crse_cell = .false.
          as_crse_covered_cell = .false.
          if (as_crse .eq. 1) then
                as_crse_crse_cell = is_inside(i,j,k,lo,hi) .and. &
                           rr_flag_crse(i,j,k) .eq. crse_fine_boundary_cell
                as_crse_covered_cell = rr_flag_crse(i,j,k) .eq. covered_by_fine
          end if

          as_fine_valid_cell = .false.  ! valid cells near box boundary
          as_fine_ghost_cell = .false.  ! ghost cells just outside valid region
          if (as_fine .eq. 1) then
            as_fine_valid_cell = is_inside(i,j,k,lo,hi)
            as_fine_ghost_cell = levmsk(i,j,k) .eq. levmsk_notcovered ! not covered by other grids
          end if

         do kk = -1,1
            do jj = -1,1
               do ii = -1,1
                  if((ii.ne. 0 .or. jj.ne.0 .or. kk.ne. 0) .and. nbr(ii,jj,kk).eq.1) then

                      iii = i + ii
                      jjj = j + jj
                      kkk = k + kk

                      drho = dM(L)*sum_kappa_inv*W(iii,jjj,kkk)

                      valid_dst_cell = is_inside(iii,jjj,kkk,lo,hi)

                      if (as_crse_crse_cell) then
                            if (rr_flag_crse(iii,jjj,kkk).eq.covered_by_fine &
                                       .and. vf(i,j,k).gt.reredistribution_threshold) then
                                rr_drho_crse(i,j,k,n) = rr_drho_crse(i,j,k,n) &
                                          + dt*drho*(vf(iii,jjj,kkk)/vf(i,j,k))
                            end if
                      end if

                      if (as_crse_covered_cell) then
                            if (valid_dst_cell) then
                               if (rr_flag_crse(iii,jjj,kkk).eq.crse_fine_boundary_cell &
                                          .and. vf(iii,jjj,kkk).gt.reredistribution_threshold) then
                                        ! the recipient is a crse/fine boundary cell
                                    rr_drho_crse(iii,jjj,kkk,n) = rr_drho_crse(iii,jjj,kkk,n) &
                                             - dt*drho
                                end if
                            end if
                     end if

                    if (as_fine_valid_cell) then
                            if (.not.valid_dst_cell) then
                                dm_as_fine(iii,jjj,kkk,n) = dm_as_fine(iii,jjj,kkk,n) &
                                          + dt*drho*vf(iii,jjj,kkk)
                            end if
                    end if

                    if (as_fine_ghost_cell) then
                            if (valid_dst_cell) then
                                     dm_as_fine(i,j,k,n) = dm_as_fine(i,j,k,n) &
                                         - dt*drho*vf(iii,jjj,kkk)
                            end if
                    end if

                    endif
                enddo
             enddo
           enddo !end neighbor loop (re redistribution)

        endif
       enddo !redistribution loop

    enddo !component loop

  end subroutine pc_fix_div_and_redistribute

  subroutine pc_set_body_state(lo, hi, S, Slo, Shi, mask, mlo, mhi, b, nc, bval) &
       bind(C,name="pc_set_body_state")

    integer,          intent(in   ) :: nc, bval
    integer,          intent(in   ) :: lo(1:3),hi(1:3)
    integer,          intent(in   ) :: Slo(1:3),Shi(1:3)
    integer,          intent(in   ) :: mlo(1:3),mhi(1:3)
    integer,          intent(in   ) :: mask(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(inout) :: S(Slo(1):Shi(1),Slo(2):Shi(2),Slo(3):Shi(3),1:nc)
    real(amrex_real), intent(in   ) :: b(1:nc)
    integer :: i,j,k,n

    do n=1,nc
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                if (mask(i,j,k).eq.bval) S(i,j,k,n)=b(n)
             enddo
          enddo
       enddo
    enddo

  end subroutine pc_set_body_state

  subroutine pc_fill_sv_ebg(lo, hi, ebg, Nebg, vfrac, vflo, vfhi, bcent, blo, bhi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi) &
       bind(C,name="pc_fill_sv_ebg")


    integer,          intent(in   ) ::  lo(0:2),  hi(0:2)
    integer, dimension(3), intent(in) :: vflo, vfhi, blo, bhi, axlo, axhi, aylo, ayhi, azlo, azhi
    integer,          intent(in   ) :: Nebg
    type(eb_bndry_geom), intent(inout) :: ebg(0:Nebg-1)
    integer :: i,j,k,L,n,ii,jj,kk


    real(amrex_real), intent(in) :: vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    real(amrex_real), intent(in) :: bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(amrex_real), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    real(amrex_real) :: axm, axp, aym, ayp, azm, azp, apnorm, apnorminv
    do L = 0, Nebg-1
       i = ebg(L) % iv(0)
       j = ebg(L) % iv(1)
       k = ebg(L) % iv(2)
       if (i.ge.lo(0) .and. i.le.hi(0) &
            .and. j.ge.lo(1) .and. j.le.hi(1) &
            .and. k.ge.lo(2) .and. k.le.hi(2) ) then

          axm = apx(i,j,k)
          axp = apx(i+1,j,k)
          aym = apy(i,j,k)
          ayp = apy(i,j+1,k)
          azm = apz(i,j,k)
          azp = apz(i,j,k+1)

          apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)
          if (apnorm .eq. 0.d0 ) then
             write(0,*) 'Cell id: ', i, ',',j,',',k
             write(0,*) ' box: ', lo(0), ',', lo(1), ',', lo(2), '; ', hi(0), ',', hi(1), ',', hi(2)
             write(0,*) 'Volume fraction: ', vfrac(i,j,k)
             write(0,*) axm, axp, aym, ayp, azm, azp
             write(0,*) 'L = ', L, ' out of ', Nebg-1
             call amrex_abort("pc_fill_sv_ebg: zero apnorm")
          end if
          apnorminv = -1.d0 / apnorm
          ebg(L) % eb_normal(1) = (axm-axp) * apnorminv  ! pointing to the wall
          ebg(L) % eb_normal(2) = (aym-ayp) * apnorminv
          ebg(L) % eb_normal(3) = (azm-azp) * apnorminv

          ebg(L) % eb_area = apnorm

          ebg(L) % eb_centroid(1) = bcent(i,j,k,1)
          ebg(L) % eb_centroid(2) = bcent(i,j,k,2)
          ebg(L) % eb_centroid(3) = bcent(i,j,k,3)

          ebg(L) % eb_vfrac = vfrac(i,j,k)
       endif
    enddo

  end subroutine pc_fill_sv_ebg

  subroutine pc_set_synthetic_data(lo,hi,M,m_lo,m_hi,prob_lo,dx) &
       bind(C, name="pc_set_synthetic_data")

    use amrex_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3),hi(3)
    integer         , intent(in   ) :: m_lo(3),m_hi(3)
    double precision, intent(inout) :: M(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    double precision, intent(in) :: prob_lo(3), dx
    integer          :: i,j,k
    double precision:: x, y

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             x = (i-prob_lo(1))*dx
             y = (j-prob_lo(2))*dx
             M(i,j,k) = 4.0*(x*x + y*y)
          enddo
       enddo
    enddo

  end subroutine pc_set_synthetic_data

end module nbrsTest_nd_module
