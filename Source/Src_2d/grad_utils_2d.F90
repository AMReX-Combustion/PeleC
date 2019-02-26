module grad_utils_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim

#ifdef PELEC_USE_EB
  use amrex_ebcellflag_module, only : get_neighbor_cells
  use pelec_eb_stencil_types_module, only : eb_bndry_geom
#endif

  implicit none

  private

  public :: pc_compute_tangential_vel_derivs

#ifdef PELEC_USE_EB
  public :: pc_compute_tangential_vel_derivs_eb
#endif

contains

  subroutine pc_compute_tangential_vel_derivs(lo,  hi, dlo, dhi,&
       Q,   Qlo,   Qhi,&
       td,  tdlo,  tdhi,&
       deltax, idir) bind(C, name = "pc_compute_tangential_vel_derivs")

    use prob_params_module, only : physbc_lo, physbc_hi
    use meth_params_module, only : QVAR, QU, QV
    use amrex_constants_module

    implicit none

    integer, intent(in   ) ::   lo(2),  hi(2)
    integer, intent(in   ) ::  dlo(2), dhi(2)
    integer, intent(in   ) ::  Qlo(2), Qhi(2)
    integer, intent(inout) :: tdlo(2),tdhi(2)
    integer, intent(in   ) :: idir
    double precision, intent(in   ) ::  Q( Qlo(1): Qhi(1), Qlo(2): Qhi(2),QVAR)
    double precision, intent(inout) :: td(tdlo(1):tdhi(1),tdlo(2):tdhi(2),dim)
    double precision, intent(in   ) :: deltax(2)
    
    integer :: i, j
    double precision :: dxinv(2)

    dxinv = 1.d0/deltax
    if (idir .eq. 0) then
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             td(i,j,1:dim) = FOURTH*dxinv(2)*(Q(i,j+1,QU:QV)+Q(i-1,j+1,QU:QV)-Q(i,j-1,QU:QV)-Q(i-1,j-1,QU:QV))
          enddo
       enddo
    else
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             td(i,j,1:dim) = FOURTH*dxinv(1)*(Q(i+1,j,QU:QV)+Q(i+1,j-1,QU:QV)-Q(i-1,j,QU:QV)-Q(i-1,j-1,QU:QV))
          enddo
       enddo
    endif
    
  end subroutine pc_compute_tangential_vel_derivs

#ifdef PELEC_USE_EB
  subroutine pc_compute_tangential_vel_derivs_eb(lo, hi, dlo, dhi,&
       sv_ebg, Ncut, &
       Q,   Qlo,   Qhi,&
       td,  tdlo,  tdhi,&
       flag, fglo, fghi, &
       deltax, idir) bind(C, name = "pc_compute_tangential_vel_derivs_eb")

    use prob_params_module, only : physbc_lo, physbc_hi
    use meth_params_module, only : QVAR, QU, QV
    use amrex_constants_module
    use amrex_ebcellflag_module, only : get_neighbor_cells

    implicit none

    integer, intent(in   ) ::   lo(0:1),  hi(0:1)
    integer, intent(in   ) ::  dlo(0:1), dhi(0:1)
    integer, intent(in   ) ::  Qlo(0:1), Qhi(0:1)
    integer, intent(inout) :: tdlo(0:1),tdhi(0:1)
    integer, intent(in   ) :: fglo(0:1),fghi(0:1)
    integer, intent(in   ) :: idir, Ncut
    type(eb_bndry_geom), intent(in   ) :: sv_ebg(0:Ncut-1)
    double precision, intent(in   ) ::  Q( Qlo(0): Qhi(0), Qlo(1): Qhi(1),QVAR)
    double precision, intent(inout) :: td(tdlo(0):tdhi(0),tdlo(1):tdhi(1),dim)
    integer,          intent(in   ) :: flag(fglo(0):fghi(0),fglo(1):fghi(1))
    double precision, intent(in   ) :: deltax(0:1)
    
    integer :: i, j, L, ii, jj, cnt
    double precision :: dxinv(0:1)
    integer :: nbrLO(-1:1,-1:1), nbrHI(-1:1,-1:1), nxy(-1:0,-1:1), nyx(-1:1, -1:0), g(-1:0,-1:0,dim)

    ! Find 6 cells surround each face that are connected
    !   to the cells on both sides of interface
    dxinv = 1.d0/deltax
    if (idir .eq. 0) then
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) ) then

             call get_neighbor_cells(flag(i-1,j),nbrLO)
             call get_neighbor_cells(flag(i  ,j),nbrHI)
             nxy = nbrLO(0:1,-1:1) * nbrHI(-1:0,-1:1)

             cnt = 0
             g = 0
             do jj=-1,0
                do ii=-1,0
                   g(ii,jj,:) = g(ii,jj,:) + &
                        (        Q(i+ii,j+jj+1,QU:QV) -   Q(i+ii,j+jj,QU:QV))*dxinv(1) &
                        *      nxy(  ii,  jj+1)       * nxy(  ii,  jj)
                   cnt = cnt + nxy(  ii,  jj+1)       * nxy(  ii,  jj)
                enddo
             enddo
             td(i,j,1:dim) = sum(g) / MAX(1,cnt)
          endif
       enddo
    else
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          if (i.ge.lo(0) .and. i.le.hi(0) &
               .and. j.ge.lo(1) .and. j.le.hi(1) ) then

             call get_neighbor_cells(flag(i,j-1),nbrLO)
             call get_neighbor_cells(flag(i,j  ),nbrHI)
             nyx = nbrLO(-1:1,0:1) * nbrHI(-1:1,-1:0)

             cnt = 0
             g = 0
             do jj=-1,0
                do ii=-1,0
                   g(ii,jj,:) = g(ii,jj,:) + &
                        (        Q(i+ii+1,j+jj,QU:QV) -   Q(i+ii,j+jj,QU:QV))*dxinv(0) &
                        *      nyx(  ii+1,  jj)       * nyx(  ii,  jj)
                   cnt = cnt + nyx(  ii+1,  jj)       * nyx(  ii,  jj)
                enddo
             enddo
             td(i,j,1:dim) = sum(g) / MAX(1,cnt)
          endif
       enddo
    endif
  end subroutine pc_compute_tangential_vel_derivs_eb
#endif

end module grad_utils_module
