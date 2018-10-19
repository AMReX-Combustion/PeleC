module pmf_module
  implicit none
  character(len=72) :: pmf_filename = ''
  integer :: pmf_init = 0, pmf_M, pmf_N, pmf_do_average = 0
  double precision, allocatable :: pmf_X(:)
  double precision, allocatable ::  pmf_Y(:,:)
  character (len=20), allocatable :: pmf_names(:)

contains

  subroutine read_pmf()
    character (len=20) :: ctmp1, ctmp2
    character (len=1) :: ctmp
    integer index, found, ist, reason, i, j, lsize, pos1, pos2, n

    character(len=:), allocatable :: line

    !     Read 2 header lines, first looks like VARIABLES = NAME1 NAME2 NAME3..., we dont care about second
    open(unit=32,file=pmf_filename,status='old')

    reason = 0

    !  Get length of first line
    lsize = 1
    do while (.true.)
       read(32,'(A)',EOR=5,END=5,ADVANCE='NO') ctmp
       lsize = lsize+1
    enddo
5   allocate(character(len=lsize) :: line)
    rewind(32)
    read(32,'(A)') line

    ! Get number of variables
    ist = INDEX(line, "=")
    line = trim(line(ist+1:))

    pmf_M = 0
    pos1 = 1
    do while (line(pos1:pos1) .eq. ' ')
       pos1 = pos1+1
    enddo
    DO
       pos2 = INDEX(line(pos1:), " ")
       IF (pos2 == 0) THEN
          pmf_M = pmf_M + 1
          EXIT
       END IF
       pmf_M = pmf_M + 1
       pos1 = pos2+pos1
       do while (line(pos1:pos1) .eq. ' ')
          pos1 = pos1+1
       enddo
    END DO
    allocate(pmf_names(pmf_M))
    pmf_M = pmf_M - 1 ! remove the X 

    ! Get names
    n = 0
    pos1 = 1
    do while (line(pos1:pos1) .eq. ' ')
       pos1 = pos1+1
    enddo

    DO
       pos2 = INDEX(line(pos1:), " ")
       IF (pos2 == 0) THEN
          n = n+1
          pmf_names(n) = line(pos1:)
          EXIT
       END IF
       n = n+1
       pmf_names(n) = line(pos1:pos1+pos2-2)
       pos1 = pos2+pos1
       do while (line(pos1:pos1) .eq. ' ')
          pos1 = pos1+1
       enddo
    END DO

    read(32,*) line! Throw away this line

    !  Count the lines of data
    pmf_N = 0
    do while (.true.)
       read(32,*,END=2) line
       pmf_N = pmf_N+1
    enddo
2   continue

    ! Strip quotes off names
    do n=1,pmf_M+1
       if (pmf_names(n)(1:1) .eq. '"') pmf_names(n) = pmf_names(n)(2:)
       ist = len(trim(pmf_names(n)))
       if (pmf_names(n)(ist:ist) .eq. '"') pmf_names(n) = pmf_names(n)(:ist-1)
    enddo

    allocate(pmf_X(pmf_N))
    allocate(pmf_Y(pmf_N,pmf_M))
    
    !  Now read data
    rewind(32)
    read(32,'(A)',advance='yes') line
    read(32,'(A)',advance='yes') line
    do i = 1,pmf_N
       read(32,*) pmf_X(i),(pmf_Y(i,j),j=1,pmf_M)
    enddo

    !  Now mark that we have read the data
    pmf_init = 1
  end subroutine read_pmf

  function pmf_ncomp() result(ncomp)
    integer :: ncomp
    ncomp = pmf_M
  end function pmf_ncomp

  function pmf_npts() result(npts)
    integer :: npts
    npts = pmf_N
  end function pmf_npts

  subroutine interp_pmf(xlo,xhi,y_vector,M_ret)
    double precision xlo,xhi,y_vector(*)
    double precision sum,xmid
    integer i,j,k,lo_loside,lo_hiside       
    integer hi_loside,hi_hiside,loside,hiside
    double precision ylo,yhi,x1,y1,x2,y2,dydx
    integer M_ret

    M_ret = pmf_M
    if (pmf_init .eq. 0) then
       if (trim(pmf_filename) .eq. '') then
          call bl_abort('must set pmf_filename prior to calling interp_pmf')
       endif
       call read_pmf()
    endif

    if (pmf_do_average .eq.1) then
       lo_loside = 0
       lo_hiside = 0
       hi_loside = 0
       hi_hiside = 0
       if (xlo .le. pmf_X(1)) then
          lo_loside = 1
          lo_hiside = 1
       end if
       if (xhi .le. pmf_X(1)) then
          hi_loside = 1
          hi_hiside = 1
       end if
       if (xlo .ge. pmf_X(pmf_N)) then
          lo_loside = pmf_N
          lo_hiside = pmf_N
       end if
       if (xhi .ge. pmf_X(pmf_N)) then
          hi_loside = pmf_N
          hi_hiside = pmf_N
       end if
       if (lo_loside.eq.0) then
          do i = 1, pmf_N-1                           
             if ( (xlo .ge. pmf_X(i)) .and.  &
                  (xlo .le. pmf_X(i+1)) ) then
                lo_loside  = i
                lo_hiside  = i+1
             end if
          end do
       end if
       if (hi_loside.eq.0) then            
          do i = 1, pmf_N-1                           
             if ( (xhi .ge. pmf_X(i)) .and. &
                  (xhi .le. pmf_X(i+1)) ) then
                hi_loside = i
                hi_hiside = i + 1
             end if
          end do
       end if

       do j = 1, pmf_M

          x1 = pmf_X(lo_loside)
          y1 = pmf_Y(lo_loside,j)

          x2 = pmf_X(lo_hiside)
          y2 = pmf_Y(lo_hiside,j)

          if (lo_loside.eq.lo_hiside) then
             dydx = 0.d0
          else
             dydx = (y2-y1)/(x2-x1)
          end if

          ylo = y1 + dydx*(xlo - x1)

          if (lo_loside .eq. hi_loside) then

             yhi = y1 + dydx*(xhi - x1)

             y_vector(j) = 0.5d0*(ylo + yhi)

          else

             sum = (x2 - xlo) * 0.5d0 * (ylo + y2)

             x1 = pmf_X(hi_loside)
             y1 = pmf_Y(hi_loside,j)

             x2 = pmf_X(hi_hiside)
             y2 = pmf_Y(hi_hiside,j)

             if (hi_loside.eq.hi_hiside) then
                dydx = 0.d0
             else
                dydx = (y2-y1)/(x2-x1)
             end if

             yhi = y1 + dydx*(xhi - x1)

             sum = sum + (xhi - x1)*0.5d0*(yhi+y1)

             do k = lo_hiside,hi_loside-1

                sum = sum + (pmf_X(k+1)-pmf_X(k)) * 0.5d0 &
                     * (pmf_Y(k,j) + pmf_Y(k+1,j))

             end do

             y_vector(j) = sum / (xhi - xlo)

          end if
       end do
    else
       xmid = 0.5d0*(xlo + xhi)
       loside = 0
       hiside = 0
       if (xmid .le. pmf_X(1)) then
          loside = 1
          hiside = 1
       end if
       if (xmid .ge. pmf_X(pmf_N)) then
          loside = pmf_N
          hiside = pmf_N
       end if
       if (loside.eq.0) then
          do i = 1, pmf_N-1                           
             if ( (xmid .ge. pmf_X(i)) .and.  &
                  (xmid .le. pmf_X(i+1)) ) then
                loside  = i
                hiside  = i+1
             end if
          end do
       end if

       do j = 1, pmf_M

          x1 = pmf_X(loside)
          y1 = pmf_Y(loside,j)

          x2 = pmf_X(hiside)
          y2 = pmf_Y(hiside,j)

          if (loside.eq.hiside) then
             dydx = 0.d0
          else
             dydx = (y2-y1)/(x2-x1)
          end if

          y_vector(j) = y1 + dydx*(xlo - x1)
       end do
    endif
  end subroutine interp_pmf

end module pmf_module

subroutine initialize_pmf(filename)
  use pmf_module
  use chemistry_module, only : nspecies, get_species_index
  character (len=*) :: filename
  integer :: n
  pmf_filename = filename
  call read_pmf()
  if (pmf_ncomp() .ne. nspecies + 3) then
     print *,'Number of dependent variables in file:',pmf_ncomp()
     print *,'Number expected:',nspecies+3
     stop 'pmf data file not compatible with current chemistry model, wrong number of species'
  endif

  if (pmf_names(1) .ne. "X")  stop 'pmf data file not compatible with pmf data reader, X must be first variable'
  if (pmf_names(2) .ne. "temp")  stop 'pmf data file not compatible with pmf data reader, temp must be second variable'
  if (pmf_names(3) .ne. "u")  stop 'pmf data file not compatible with pmf data reader, u must be third variable'
  do n=1,nspecies
     if (n .ne. get_species_index(pmf_names(4+n))) then
        stop 'pmf data file not compatible with current chemistry model, wrong species'
     endif
  enddo
end subroutine initialize_pmf

subroutine pmf(xlo,xhi,y_vector,M)
  use pmf_module
  double precision :: xlo,xhi,y_vector(*)
  integer :: M
  call interp_pmf(xlo,xhi,y_vector,M)
end subroutine pmf

