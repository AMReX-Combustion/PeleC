module grad_utils_module

  use amrex_fort_module, only : rt=>amrex_real, dim=>bl_spacedim

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
    use meth_params_module, only : QVAR, QU, QV, QW
    use amrex_constants_module

    implicit none

    integer, intent(in   ) ::   lo(3),  hi(3)
    integer, intent(in   ) ::  dlo(3), dhi(3)
    integer, intent(in   ) ::  Qlo(3), Qhi(3)
    integer, intent(inout) :: tdlo(3),tdhi(3)
    integer, intent(in   ) :: idir
    real(rt), intent(in   ) ::  Q( Qlo(1): Qhi(1), Qlo(2): Qhi(2), Qlo(3): Qhi(3),QVAR)
    real(rt), intent(inout) :: td(tdlo(1):tdhi(1),tdlo(2):tdhi(2),tdlo(3):tdhi(3),6)
    real(rt), intent(in   ) :: deltax(3)

    integer :: i, j, k
    real(rt) :: dxinv(3)

    dxinv = 1.d0/deltax
    if (idir .eq. 0) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                td(i,j,k,1) = FOURTH*dxinv(2)*(Q(i,j+1,k,QU)+Q(i-1,j+1,k,QU)-Q(i,j-1,k,QU)-Q(i-1,j-1,k,QU))
                td(i,j,k,2) = FOURTH*dxinv(2)*(Q(i,j+1,k,QV)+Q(i-1,j+1,k,QV)-Q(i,j-1,k,QV)-Q(i-1,j-1,k,QV))
                td(i,j,k,3) = FOURTH*dxinv(2)*(Q(i,j+1,k,QW)+Q(i-1,j+1,k,QW)-Q(i,j-1,k,QW)-Q(i-1,j-1,k,QW))
                td(i,j,k,4) = FOURTH*dxinv(3)*(Q(i,j,k+1,QU)+Q(i-1,j,k+1,QU)-Q(i,j,k-1,QU)-Q(i-1,j,k-1,QU))
                td(i,j,k,5) = FOURTH*dxinv(3)*(Q(i,j,k+1,QV)+Q(i-1,j,k+1,QV)-Q(i,j,k-1,QV)-Q(i-1,j,k-1,QV))
                td(i,j,k,6) = FOURTH*dxinv(3)*(Q(i,j,k+1,QW)+Q(i-1,j,k+1,QW)-Q(i,j,k-1,QW)-Q(i-1,j,k-1,QW))
             enddo
          enddo
       enddo

    else if (idir .eq. 1) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                td(i,j,k,1) = FOURTH*dxinv(1)*(Q(i+1,j,k,QU)+Q(i+1,j-1,k,QU)-Q(i-1,j,k,QU)-Q(i-1,j-1,k,QU))
                td(i,j,k,2) = FOURTH*dxinv(1)*(Q(i+1,j,k,QV)+Q(i+1,j-1,k,QV)-Q(i-1,j,k,QV)-Q(i-1,j-1,k,QV))
                td(i,j,k,3) = FOURTH*dxinv(1)*(Q(i+1,j,k,QW)+Q(i+1,j-1,k,QW)-Q(i-1,j,k,QW)-Q(i-1,j-1,k,QW))
                td(i,j,k,4) = FOURTH*dxinv(3)*(Q(i,j,k+1,QU)+Q(i,j-1,k+1,QU)-Q(i,j,k-1,QU)-Q(i,j-1,k-1,QU))
                td(i,j,k,5) = FOURTH*dxinv(3)*(Q(i,j,k+1,QV)+Q(i,j-1,k+1,QV)-Q(i,j,k-1,QV)-Q(i,j-1,k-1,QV))
                td(i,j,k,6) = FOURTH*dxinv(3)*(Q(i,j,k+1,QW)+Q(i,j-1,k+1,QW)-Q(i,j,k-1,QW)-Q(i,j-1,k-1,QW))
             enddo
          enddo
       enddo

    else
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                td(i,j,k,1) = FOURTH*dxinv(1)*(Q(i+1,j,k,QU)+Q(i+1,j,k-1,QU)-Q(i-1,j,k,QU)-Q(i-1,j,k-1,QU))
                td(i,j,k,2) = FOURTH*dxinv(1)*(Q(i+1,j,k,QV)+Q(i+1,j,k-1,QV)-Q(i-1,j,k,QV)-Q(i-1,j,k-1,QV))
                td(i,j,k,3) = FOURTH*dxinv(1)*(Q(i+1,j,k,QW)+Q(i+1,j,k-1,QW)-Q(i-1,j,k,QW)-Q(i-1,j,k-1,QW))
                td(i,j,k,4) = FOURTH*dxinv(2)*(Q(i,j+1,k,QU)+Q(i,j+1,k-1,QU)-Q(i,j-1,k,QU)-Q(i,j-1,k-1,QU))
                td(i,j,k,5) = FOURTH*dxinv(2)*(Q(i,j+1,k,QV)+Q(i,j+1,k-1,QV)-Q(i,j-1,k,QV)-Q(i,j-1,k-1,QV))
                td(i,j,k,6) = FOURTH*dxinv(2)*(Q(i,j+1,k,QW)+Q(i,j+1,k-1,QW)-Q(i,j-1,k,QW)-Q(i,j-1,k-1,QW)) 
             enddo
          enddo
       enddo
    endif
    
  end subroutine pc_compute_tangential_vel_derivs

#ifdef PELEC_USE_EB
  subroutine pc_compute_tangential_vel_derivs_eb(lo, hi, dlo, dhi, &
       sv_ebg, Ncut, &
       Q,   Qlo,   Qhi,&
       td,  tdlo,  tdhi,&
       flag, fglo, fghi, &
       deltax, idir) bind(C, name = "pc_compute_tangential_vel_derivs_eb")

    use prob_params_module, only : physbc_lo, physbc_hi
    use meth_params_module, only : QVAR, QU, QV, QW
    use amrex_constants_module
    use amrex_ebcellflag_module, only : get_neighbor_cells

    implicit none

    integer, intent(in   ) ::   lo(3),  hi(3)
    integer, intent(in   ) ::  dlo(3), dhi(3)
    integer, intent(in   ) ::  Qlo(3), Qhi(3)
    integer, intent(inout) :: tdlo(3),tdhi(3)
    integer, intent(in   ) :: fglo(3),fghi(3)
    integer, intent(in   ) :: idir, Ncut
    type(eb_bndry_geom), intent(in   ) :: sv_ebg(0:Ncut-1)
    real(rt), intent(in   ) ::  Q( Qlo(1): Qhi(1), Qlo(2): Qhi(2), Qlo(3): Qhi(3),QVAR)
    real(rt), intent(inout) :: td(tdlo(1):tdhi(1),tdlo(2):tdhi(2),tdlo(3):tdhi(3),6)
    integer,  intent(in   ) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in   ) :: deltax(3)

    integer :: i, j, k, L, ii, jj, kk, cnt
    real(rt) :: dxinv(3)
    integer :: nbrLO(-1:1,-1:1,-1:1), nbrHI(-1:1,-1:1,-1:1)
    integer :: nxy(-1:0,-1:1), nxz(-1:0,-1:1)
    integer :: nyx(-1:1,-1:0), nyz(-1:0,-1:1)
    integer :: nzx(-1:1,-1:0), nzy(-1:1,-1:0)

    integer  :: ihip, ihim, ilop, ilom, jhip, jhim, jlop, jlom, khip, khim, klop, klom
    real(rt) :: wlo, whi
    real(rt), parameter :: weights(0:2) = [0.d0, 1.d0, 0.5d0]

    dxinv = 1.d0/deltax
    if (idir .eq. 0) then
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (i.ge.lo(1) .and. i.le.hi(1) &
               .and. j.ge.lo(2) .and. j.le.hi(2) &
               .and. k.ge.lo(3) .and. k.le.hi(3) ) then

             call get_neighbor_cells(flag(i-1,j,k),nbrLO)
             call get_neighbor_cells(flag(i  ,j,k),nbrHI)

             jhip = j + nbrHI(0, 1,0)
             jhim = j - nbrHI(0,-1,0)
             jlop = j + nbrLO(0, 1,0)
             jlom = j - nbrLO(0,-1,0)
             whi = weights(jhip-jhim)
             wlo = weights(jlop-jlom)
             td(i,j,k,1) = (0.5d0*dxinv(2)) * &
                  ((Q(i  ,jhip,k,QU)-Q(i  ,jhim,k,QU))*whi &
                  +(Q(i-1,jlop,k,QU)-Q(i-1,jlom,k,QU))*wlo)
             td(i,j,k,2) = (0.5d0*dxinv(2)) * &
                  ((Q(i  ,jhip,k,QV)-Q(i  ,jhim,k,QV))*whi &
                  +(Q(i-1,jlop,k,QV)-Q(i-1,jlom,k,QV))*wlo)
             td(i,j,k,3) = (0.5d0*dxinv(2)) * &
                  ((Q(i  ,jhip,k,QW)-Q(i  ,jhim,k,QW))*whi &
                  +(Q(i-1,jlop,k,QW)-Q(i-1,jlom,k,QW))*wlo)

             khip = k + nbrHI(0,0, 1)
             khim = k - nbrHI(0,0,-1)
             klop = k + nbrLO(0,0, 1)
             klom = k - nbrLO(0,0,-1)
             whi = weights(khip-khim)
             wlo = weights(klop-klom)
             td(i,j,k,4) = (0.5d0*dxinv(3)) * &
                  ((Q(i  ,j,khip,QU)-Q(i  ,j,khim,QU))*whi &
                  +(Q(i-1,j,klop,QU)-Q(i-1,j,klom,QU))*wlo)
             td(i,j,k,5) = (0.5d0*dxinv(3)) * &
                  ((Q(i  ,j,khip,QV)-Q(i  ,j,khim,QV))*whi &
                  +(Q(i-1,j,klop,QV)-Q(i-1,j,klom,QV))*wlo)
             td(i,j,k,6) = (0.5d0*dxinv(3)) * &
                  ((Q(i  ,j,khip,QW)-Q(i  ,j,khim,QW))*whi &
                  +(Q(i-1,j,klop,QW)-Q(i-1,j,klom,QW))*wlo)

          endif
       enddo
    else if (idir .eq. 1) then
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (i.ge.lo(1) .and. i.le.hi(1) &
               .and. j.ge.lo(2) .and. j.le.hi(2) &
               .and. k.ge.lo(3) .and. k.le.hi(3) ) then

             call get_neighbor_cells(flag(i,j-1,k),nbrLO)
             call get_neighbor_cells(flag(i,j  ,k),nbrHI)

             ihip = i + nbrHI( 1,0,0)
             ihim = i - nbrHI(-1,0,0)
             ilop = i + nbrLO( 1,0,0)
             ilom = i - nbrLO(-1,0,0)
             whi = weights(ihip-ihim)
             wlo = weights(ilop-ilom)
             td(i,j,k,1) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j  ,k,QU)-Q(ihim,j  ,k,QU))*whi &
                  +(Q(ilop,j-1,k,QU)-Q(ilom,j-1,k,QU))*wlo)
             td(i,j,k,2) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j  ,k,QV)-Q(ihim,j  ,k,QV))*whi &
                  +(Q(ilop,j-1,k,QV)-Q(ilom,j-1,k,QV))*wlo)
             td(i,j,k,3) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j  ,k,QW)-Q(ihim,j  ,k,QW))*whi &
                  +(Q(ilop,j-1,k,QW)-Q(ilom,j-1,k,QW))*wlo)

             khip = k + nbrHI(0,0, 1)
             khim = k - nbrHI(0,0,-1)
             klop = k + nbrLO(0,0, 1)
             klom = k - nbrLO(0,0,-1)
             whi = weights(khip-khim)
             wlo = weights(klop-klom)
             td(i,j,k,4) = (0.5d0*dxinv(3)) * &
                  ((Q(i,j  ,khip,QU)-Q(i,j  ,khim,QU))*whi &
                  +(Q(i,j-1,klop,QU)-Q(i,j-1,klom,QU))*wlo)
             td(i,j,k,5) = (0.5d0*dxinv(3)) * &
                  ((Q(i,j  ,khip,QV)-Q(i,j  ,khim,QV))*whi &
                  +(Q(i,j-1,klop,QV)-Q(i,j-1,klom,QV))*wlo)
             td(i,j,k,6) = (0.5d0*dxinv(3)) * &
                  ((Q(i,j  ,khip,QW)-Q(i,j  ,khim,QW))*whi &
                  +(Q(i,j-1,klop,QW)-Q(i,j-1,klom,QW))*wlo)
          endif
       enddo
    else if (idir .eq. 2) then
       do L = 0, Ncut-1
          i = sv_ebg(L) % iv(0)
          j = sv_ebg(L) % iv(1)
          k = sv_ebg(L) % iv(2)
          if (i.ge.lo(1) .and. i.le.hi(1) &
               .and. j.ge.lo(2) .and. j.le.hi(2) &
               .and. k.ge.lo(3) .and. k.le.hi(3) ) then

             call get_neighbor_cells(flag(i,j,k-1),nbrLO)
             call get_neighbor_cells(flag(i,j,k  ),nbrHI)

             ihip = i + nbrHI( 1,0,0)
             ihim = i - nbrHI(-1,0,0)
             ilop = i + nbrLO( 1,0,0)
             ilom = i - nbrLO(-1,0,0)
             whi = weights(ihip-ihim)
             wlo = weights(ilop-ilom)
             td(i,j,k,1) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j,k  ,QU)-Q(ihim,j,k  ,QU))*whi &
                  +(Q(ilop,j,k-1,QU)-Q(ilom,j,k-1,QU))*wlo)
             td(i,j,k,2) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j,k  ,QV)-Q(ihim,j,k  ,QV))*whi &
                  +(Q(ilop,j,k-1,QV)-Q(ilom,j,k-1,QV))*wlo)
             td(i,j,k,3) = (0.5d0*dxinv(1)) * &
                  ((Q(ihip,j,k  ,QW)-Q(ihim,j,k  ,QW))*whi &
                  +(Q(ilop,j,k-1,QW)-Q(ilom,j,k-1,QW))*wlo)

             jhip = j + nbrHI(0, 1,0)
             jhim = j - nbrHI(0,-1,0)
             jlop = j + nbrLO(0, 1,0)
             jlom = j - nbrLO(0,-1,0)
             whi = weights(jhip-jhim)
             wlo = weights(jlop-jlom)
             td(i,j,k,4) = (0.5d0*dxinv(2)) * &
                  ((Q(i,jhip,k  ,QU)-Q(i,jhim,k  ,QU))*whi &
                  +(Q(i,jlop,k-1,QU)-Q(i,jlom,k-1,QU))*wlo)
             td(i,j,k,5) = (0.5d0*dxinv(2)) * &
                  ((Q(i,jhip,k  ,QV)-Q(i,jhim,k  ,QV))*whi &
                  +(Q(i,jlop,k-1,QV)-Q(i,jlom,k-1,QV))*wlo)
             td(i,j,k,6) = (0.5d0*dxinv(2)) * &
                  ((Q(i,jhip,k  ,QW)-Q(i,jhim,k  ,QW))*whi &
                  +(Q(i,jlop,k-1,QW)-Q(i,jlom,k-1,QW))*wlo)
          endif
       enddo
    endif
  end subroutine pc_compute_tangential_vel_derivs_eb
#endif

end module grad_utils_module
