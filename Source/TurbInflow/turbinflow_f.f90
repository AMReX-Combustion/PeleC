module turbinflow_module

  implicit none

  logical, save :: turbinflow_initialized = .false.

  integer, save :: lenfname
  integer, save :: iturbfile(128)

  integer, save :: npboxcells(3)
  double precision, save :: pboxlo(3), dx(3), dxinv(3)

  integer, parameter :: nplane = 32
  double precision, allocatable, save :: sdata(:,:,:,:)
  double precision, save :: szlo=0.d0, szhi=0.d0

  integer, parameter :: isswirltype = 0  ! periodic
  double precision, save :: units_conversion = 100.d0  ! m --> cm & m/s --> cm/s

  private

  public :: turbinflow_initialized, init_turbinflow, get_turbvelocity

!$omp threadprivate(turbinflow_initialized,lenfname,iturbfile,npboxcells,pboxlo,dx,dxinv)
!$omp threadprivate(sdata,szlo,szhi,units_conversion)

contains

  subroutine init_turbinflow(turbfile, is_cgs)
    character (len=*), intent(in) :: turbfile
    logical, intent(in), optional :: is_cgs

    integer, parameter :: iunit = 20
    integer :: ierr, n, npts(3)
    double precision :: probsize(3), pboxsize(3)

    if (present(is_cgs)) then
       if (is_cgs) then
          units_conversion = 100.d0
       else
          units_conversion = 1.d0
       end if
    end if

    lenfname = len_trim(turbfile)
    do n=1,lenfname
       iturbfile(n) = ichar(turbfile(n:n))
    end do

    ! read header
    open(iunit, file=trim(turbfile)//'/HDR', form='formatted', action='read', &
         status='old', iostat=ierr)
    
    if (ierr .ne. 0) then
       call bl_abort('Problem opening file: ' // trim(turbfile) // '/HDR')
    end if

    read(iunit,*) npts
    read(iunit,*) probsize

    close(iunit)

    probsize = probsize * units_conversion

    dx = probsize / dble(npts-1)
    dxinv = 1.d0/dx

    pboxsize(1:2) = probsize(1:2) - 2.d0*dx(1:2)  ! because there is one ghost point on each side
    pboxsize(3) = probsize(3) ! no ghost point in z-direction

    npboxcells(1:2) = npts(1:2) - 3
    npboxcells(3)   = npts(3) - 1
    
    ! period box covers -0.5*pboxsize(1) <= x <= 0.5*pboxsize(1)
    !                   -0.5*pboxsize(2) <= y <= 0.5*pboxsize(2)
    !                                  0 <= z <= pboxsize(3)
    pboxlo(1:2) = -0.5d0*pboxsize(1:2)
    pboxlo(3) = 0.d0

    allocate(sdata(npts(1),npts(2),nplane,3))

    turbinflow_initialized = .true.

  end subroutine init_turbinflow


  subroutine get_turbvelocity(lo1,lo2,hi1,hi2,x,y,z,v)
    integer, intent(in) :: lo1,lo2,hi1,hi2
    double precision, intent(in) :: x(lo1:hi1), y(lo2:hi2)
    double precision, intent(in) :: z
    double precision, intent(out) :: v(lo1:hi1,lo2:hi2,3)

    integer :: i, j, k, n, i0, j0, k0, ii, jj
    double precision :: xx, yy, zz, zdata(0:2,0:2), ydata(0:2)
    double precision :: cx(0:2), cy(0:2), cz(0:2) 

    if (.not. turbinflow_initialized) call bl_error("turbinflow module uninitialized")

    if (z.lt.szlo+0.5d0*dx(3) .or. z.gt.szhi-0.5d0*dx(3)) then
       call store_planes(z)
    end if

    zz = (z-szlo)*dxinv(3)
    k0 = nint(zz) - 1
    zz = zz - dble(k0)
    cz(0) = 0.5d0*(zz-1.d0)*(zz-2.d0)
    cz(1) = zz*(2.d0-zz)
    cz(2) = 0.5d0*zz*(zz-1.d0)
    k0 = k0 + 1 ! because it's Fortran
    k0 = min(max(k0,1),nplane-2)

    do n=1,3
       do j=lo2,hi2
          yy = (y(j)-pboxlo(2))*dxinv(2)
          j0 = nint(yy) - 1
          yy = yy - dble(j0)
          cy(0) = 0.5d0*(yy-1.d0)*(yy-2.d0)
          cy(1) = yy*(2.d0-yy)
          cy(2) = 0.5d0*yy*(yy-1.d0)
          j0 = modulo(j0, npboxcells(2)) + 2 ! +2 because Fortran index starts with 1 and there is a ghost point
          do i=lo1,hi1
             xx = (x(i)-pboxlo(1))*dxinv(1)
             i0 = nint(xx) - 1
             xx = xx - dble(i0)
             cx(0) = 0.5d0*(xx-1.d0)*(xx-2.d0)
             cx(1) = xx*(2.d0-xx)
             cx(2) = 0.5d0*xx*(xx-1.d0)
             i0 = modulo(i0, npboxcells(1)) + 2 ! +2 as j0
             
             do jj=0,2
                do ii=0,2
                   zdata(ii,jj) = cz(0)*sdata(i0+ii,j0+jj,k0  ,n) &
                        +         cz(1)*sdata(i0+ii,j0+jj,k0+1,n) &
                        +         cz(2)*sdata(i0+ii,j0+jj,k0+2,n)
                end do
             end do

             do ii=0,2
                ydata(ii) = cy(0)*zdata(ii,0) + cy(1)*zdata(ii,1) + cy(2)*zdata(ii,2)
             end do

             v(i,j,n) = cx(0)*ydata(0) + cx(1)*ydata(1) + cx(2)*ydata(2)

          end do
       end do
    end do

  end subroutine get_turbvelocity

  subroutine store_planes(z)
    double precision, intent(in) :: z
    integer :: izlo, iplane, k, n
    izlo = nint(z*dxinv(3)) - 1
    szlo = izlo*dx(3)
    szhi = szlo + dble(nplane-1)*dx(3)
    do n=1,3
       do iplane=1,nplane
          k = modulo(izlo+iplane-1, npboxcells(3)) + 1
          call getplane(iturbfile, lenfname, sdata(:,:,iplane,n), k, n, isswirltype)
          sdata(:,:,iplane,n) = sdata(:,:,iplane,n)*units_conversion
       end do
    end do
  end subroutine store_planes

end module turbinflow_module
