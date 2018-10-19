!!$-------------------------------------------------
!!$                                                 
!!$                                                 
!!$     OVERVIEW: call pmf(xlo,xhi,y,Niny)
!!$     computes y(xlo,xhi) by averaging data over  
!!$     the range [xlo,xhi]. Here x is a scalar, y  
!!$     a vector, with Niny components used. Piece- 
!!$     constant Extrapolation is applied to points 
!!$     laying outside the range of the orig data   
!!$                                                 
!!$     Reads in the Cantera tab separated solution
!!$     (X,Y,T,rho)
!!$                                                 
!!$-------------------------------------------------
module InterpCanteraSolution
  implicit none
  integer, save :: N, M

  PRIVATE :: N, M
  
  double precision , save, allocatable :: x_data(:)
  double precision , save, allocatable :: y_data(:,:)

  character(len=128) :: FileName

contains 

  subroutine InitReadCantera
    implicit none
    integer :: iunit,i,io,j
    double precision :: temp(13)
    
    M = 12

    ! Read Cantera premixed flame solution
     FileName = 'PMF_100bar_MultiDiffusion_refine.csv'
    ! FileName = 'PMF_100bar_100K_Scaled.csv'
    iunit = 10
    OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",ACTION="READ")

    N = 0
    do
       read(iunit,*,iostat=io) temp
       if (io/=0) exit
       N = N + 1
    end do

    rewind(iunit)

    allocate(y_data(N,M))
    allocate(x_data(N))
    
    do i=1,N
       read(iunit,*) x_data(i), (y_data(i,j),j=1,M)
    end do
    
    CLOSE(UNIT=iunit)
    
  end subroutine InitReadCantera
  !==============================================!
  ! Given a x point interpolates and returns the !
  !              state space                     !
  !==============================================!
  subroutine pmf(xlo,xhi,y_vector,Niny)
    implicit none
    double precision, intent(in) :: xlo,xhi
    double precision, intent(out) :: y_vector(:)
    integer, intent(out) :: Niny
    integer :: i,j,k,lo_loside,lo_hiside       
    integer :: hi_loside,hi_hiside       
    double precision :: sum    
    double precision :: ylo,yhi,x1,y1,x2,y2,dydx   

    Niny = 12 

!!$    interpolation routine

      lo_loside = 0
      lo_hiside = 0
      hi_loside = 0
      hi_hiside = 0

      if (xlo .le. x_data(1)) then
         lo_loside = 1
         lo_hiside = 1
      end if
      
      if (xhi .le. x_data(1)) then
         hi_loside = 1
         hi_hiside = 1
      end if

      if (xlo .ge. x_data(N)) then
         lo_loside = N
         lo_hiside = N
      end if

      if (xhi .ge. x_data(N)) then
         hi_loside = N
         hi_hiside = N
      end if

      if (lo_loside.eq.0) then
         do i = 1, N-1                           
            if ( (xlo .ge. x_data(i)) .and. (xlo .le. x_data(i+1)) ) then
               lo_loside  = i
               lo_hiside  = i+1
            end if
         end do
      end if
      if (hi_loside.eq.0) then            
         do i = 1, N-1                           
            if ( (xhi .ge. x_data(i)) .and. (xhi .le. x_data(i+1)) ) then
               hi_loside = i
               hi_hiside = i + 1
            end if
         end do
      end if
      
      do j = 1, M
         
         x1 = x_data(lo_loside)
         y1 = y_data(lo_loside,j)
         
         x2 = x_data(lo_hiside)
         y2 = y_data(lo_hiside,j)
         
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
            
            x1 = x_data(hi_loside)
            y1 = y_data(hi_loside,j)
            
            x2 = x_data(hi_hiside)
            y2 = y_data(hi_hiside,j)
            
            if (hi_loside.eq.hi_hiside) then
               dydx = 0.d0
            else
               dydx = (y2-y1)/(x2-x1)
            end if
            
            yhi = y1 + dydx*(xhi - x1)
            
            sum = sum + (xhi - x1)*0.5d0*(yhi+y1)
            
            do k = lo_hiside,hi_loside-1
               
               sum = sum + (x_data(k+1)-x_data(k))* 0.5d0 * (y_data(k,j) + y_data(k+1,j))
               
            end do
            
            y_vector(j) = sum / (xhi - xlo)
            
         end if
      end do
      
    end subroutine pmf

end module InterpCanteraSolution
