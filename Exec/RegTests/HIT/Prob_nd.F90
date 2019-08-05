module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")

    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
    use probdata_module
    use amrex_fort_module
    use amrex_constants_module, only: HALF
    use network, only: nspecies, naux, molec_wt
    use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity
    use prob_params_module, only: dim
    use eos_module
    use forcing_src_module, only: u0, v0, w0, forcing

    implicit none

    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer untin,i

    namelist /fortin/ iname, binfmt, restart, lambda0, reynolds_lambda0, mach_t0, prandtl, inres, uin_norm, u0, v0, w0, forcing

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character probin*(maxlen)

    ! Local
    double precision, dimension(:), allocatable :: data
    integer, parameter :: out_unit=20
    type(eos_t) :: eos_state
    integer(kind=8) :: nx, ny, nz

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults here
    iname = ""
    binfmt = .false.
    restart = .false.
    lambda0 = 0.5_amrex_real
    reynolds_lambda0 = 100.0_amrex_real
    mach_t0 = 0.1_amrex_real
    prandtl = 0.71_amrex_real
    inres = 0
    uin_norm = 1.0_amrex_real
    u0 = 0.d0
    v0 = 0.d0
    w0 = 0.d0
    forcing = 0.d0 ! = dissipation rate / (3 * urms**2) (linear forcing, see Rosales PoF 2015)

    ! Initial pressure and temperature
    p0 = 1.013d6 ! [erg cm^-3]
    T0 = 300.d0

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    L_x = probhi(1) - problo(1)
    L_y = probhi(2) - problo(2)
    L_z = probhi(3) - problo(3)

    ! Define the molecular weight for air
    molec_wt = 28.97

    ! Wavelength associated to Taylor length scale
    k0 = 2.d0/lambda0

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Initial density, velocity, and material properties
    rho0  = eos_state % rho
    eint0 = eos_state % e
    urms0 = mach_t0 * eos_state % cs / sqrt(3.d0)
    tau  = lambda0 / urms0
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = rho0 * urms0 * lambda0 / reynolds_lambda0
    const_conductivity = const_viscosity * eos_state % cp / prandtl

    ! Write this out to file (might be useful for postprocessing)
    if ( amrex_pd_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"lambda0, k0, rho0, urms0, tau, p0, T0, gamma, mu, k, c_s0, Reynolds, Mach, Prandtl, u0, v0, w0, forcing"
       write(out_unit,*) lambda0, "," , k0, "," , eos_state % rho, "," , urms0, "," , tau, "," , &
            eos_state % p, "," , eos_state % T, "," , gamma_const, "," , const_viscosity, "," , &
            const_conductivity, "," , eos_state % cs, "," , &
            reynolds_lambda0, "," , mach_t0, "," , prandtl, "," , &
            u0, "," , v0, "," , w0, "," , forcing
       close(out_unit)
    endif

    ! Load velocity fields from file. Assume data set ordered in Fortran
    ! format and reshape the data accordingly. One thing to keep in mind
    ! is that this contains the entire input data. We will interpolate
    ! this data later to just match our box. Another assumption is that
    ! the input data is a periodic cube. If the input cube is smaller
    ! than our domain size, the cube will be repeated throughout the
    ! domain (hence the mod operations in the interpolation).

    if (restart) then
       if ( amrex_pd_ioprocessor() ) then
          write(*,*)"Skipping input file reading and assuming restart."
       endif
    else
       nx = int8(inres)
       ny = int8(1)
       nz = int8(1)
       if (dim .ge. 2) then
          ny = int8(inres)
          if (dim .ge. 3) then
             nz = int8(inres)
          endif
       endif
       allocate(data(0:nx*ny*nz*6-1))
       allocate(xinput(0:nx-1,0:ny-1,0:nz-1))
       ! allocate(yinput(0:nx-1,0:ny-1,0:nz-1))
       ! allocate(zinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(uinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(vinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(winput(0:nx-1,0:ny-1,0:nz-1))
       if (binfmt) then
          call read_binary(iname, nx, ny, nz, data)
       else
          call read_csv(iname, nx, ny, nz, data)
       endif

       uinput = urms0 / uin_norm * reshape(data(3::6), (/nx, ny, nz/))
       vinput = urms0 / uin_norm * reshape(data(4::6), (/nx, ny, nz/))
       winput = urms0 / uin_norm * reshape(data(5::6), (/nx, ny, nz/))
       xinput = reshape(data(0::6), (/nx, ny, nz/))
       ! yinput = reshape(data(1::6), (/nx, ny, nz/))
       ! zinput = reshape(data(2::6), (/nx, ny, nz/))

       ! Get the xarray table and the differences.
       allocate(xarray(0:nx-1))
       allocate(xdiff(0:nx-1))
       xarray(0:nx-1) = xinput(:,0,0)
       xdiff(:nx-2) = xarray(1:) - xarray(:nx-2)
       xdiff(nx-1) = xarray(nx-1) - xarray(nx-2)

       ! Dimensions of the input box.
       Linput = maxval(xinput(:,0,0)) + HALF*xdiff(nx-1)

       ! Deallocate some stuff
       deallocate(data)
    endif

  end subroutine amrex_probinit


  ! ::: -----------------------------------------------------------
  ! ::: This routine is called at problem setup time and is used
  ! ::: to initialize data on each grid.
  ! :::
  ! ::: NOTE:  all arrays have one cell of ghost zones surrounding
  ! :::        the grid interior.  Values in these cells need not
  ! :::        be set here.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: level     => amr level of grid
  ! ::: time      => time at which to init data
  ! ::: lo,hi     => index limits of grid interior (cell centered)
  ! ::: nstate    => number of state components.  You should know
  ! :::		   this already!
  ! ::: state     <=  Scalar array
  ! ::: delta     => cell size
  ! ::: xlo,xhi   => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name = "pc_initdata")

    use probdata_module
    use network, only: nspecies, naux, molec_wt
    use eos_type_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UFS, UTEMP
    use amrex_constants_module, only: ZERO, HALF, M_PI
    use prob_params_module, only: dim
    use eos_module
    use forcing_src_module, only: u0, v0, w0, forcing

    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3),state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)

    type(eos_t) :: eos_state
    integer :: i, j, k
    double precision :: x, y, z, u, v, w
    double precision :: xmod, ymod, zmod
    integer :: m, mp1, n, np1, p, pp1
    double precision :: rr, s, t, uinterp, vinterp, winterp
    double precision :: f0, f1, f2, f3, f4, f5, f6, f7

    ! Set the equation of state variables
    call build(eos_state)
    eos_state % p = p0
    eos_state % T = T0
    eos_state % massfrac    = 0.d0
    eos_state % massfrac(1) = 1.d0
    call eos_tp(eos_state)

    ! Uniform density, temperature, pressure (internal energy)
    state(:,:,:,URHO)            = rho0
    do i = 1,nspecies
       state(:,:,:,UFS+i-1)      = rho0 * eos_state % massfrac(i)
    enddo
    state(:,:,:,UTEMP)           = T0
    state(:,:,:,UEINT)           = rho0 * eint0

    ! Fill in the velocities and energy.
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
       zmod = mod(z,Linput)
       call locate(xarray, inres, zmod, p)
       pp1 = mod(p+1,inres)
       t = (zmod - xarray(p)) / xdiff(p)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)
          ymod = mod(y,Linput)
          call locate(xarray, inres, ymod, n)
          np1 = mod(n+1,inres)
          s = (ymod - xarray(n)) / xdiff(n)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)
             xmod = mod(x,Linput)
             call locate(xarray, inres, xmod, m)
             mp1 = mod(m+1,inres)
             rr = (xmod - xarray(m)) / xdiff(m)

             if (dim .eq. 1) then
                f0 = (1-rr)
                f1 = rr
                uinterp = uinput(m,0,0) * f0 + &
                     uinput(mp1,0,0) * f1
                vinterp = ZERO
                winterp = ZERO
             elseif (dim .eq. 2) then
                ! Factors for bilinear interpolation
                f0 = (1-rr) * (1-s)
                f1 = rr * (1-s)
                f2 = (1-rr) * s
                f3 = rr * s
                uinterp = uinput(m,n,0) * f0 + &
                     uinput(mp1,n,0) * f1 + &
                     uinput(m,np1,0) * f2 + &
                     uinput(mp1,np1,0) * f3
                vinterp = vinput(m,n,0) * f0 + &
                     vinput(mp1,n,0) * f1 + &
                     vinput(m,np1,0) * f2 + &
                     vinput(mp1,np1,0) * f3
                winterp = ZERO
             elseif (dim .eq. 3) then
                ! Factors for trilinear interpolation
                f0 = (1-rr) * (1-s) * (1-t)
                f1 = rr * (1-s) * (1-t)
                f2 = (1-rr) * s * (1-t)
                f3 = (1-rr) * (1-s) * t
                f4 = rr * (1-s) * t
                f5 = (1-rr) * s * t
                f6 = rr * s * (1-t)
                f7 = rr * s * t
                uinterp = uinput(m,n,p) * f0 + &
                     uinput(mp1,n,p) * f1 + &
                     uinput(m,np1,p) * f2 + &
                     uinput(m,n,pp1) * f3 + &
                     uinput(mp1,n,pp1) * f4 + &
                     uinput(m,np1,pp1) * f5+ &
                     uinput(mp1,np1,p) * f6 + &
                     uinput(mp1,np1,pp1) * f7
                vinterp = vinput(m,n,p) * f0 + &
                     vinput(mp1,n,p) * f1 + &
                     vinput(m,np1,p) * f2 + &
                     vinput(m,n,pp1) * f3 + &
                     vinput(mp1,n,pp1) * f4 + &
                     vinput(m,np1,pp1) * f5+ &
                     vinput(mp1,np1,p) * f6 + &
                     vinput(mp1,np1,pp1) * f7
                winterp = winput(m,n,p) * f0 + &
                     winput(mp1,n,p) * f1 + &
                     winput(m,np1,p) * f2 + &
                     winput(m,n,pp1) * f3 + &
                     winput(mp1,n,pp1) * f4 + &
                     winput(m,np1,pp1) * f5+ &
                     winput(mp1,np1,p) * f6 + &
                     winput(mp1,np1,pp1) * f7
             endif

             u = uinterp + u0
             v = vinterp + v0
             w = winterp + w0
             state(i,j,k,UMX)   = state(i,j,k,URHO) * u
             if (dim .ge. 2) then
                state(i,j,k,UMY)   = state(i,j,k,URHO) * v
                if (dim .ge. 3) then
                   state(i,j,k,UMZ)   = state(i,j,k,URHO) * w
                endif
             endif
             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * state(i,j,k,URHO) * (u**2 + v**2 + w**2)
          enddo
       enddo
    enddo

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")

  end subroutine pc_prob_close


  ! ::: -----------------------------------------------------------
  ! ::: Read a binary file
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: iname => filename
  ! ::: nx    => input resolution
  ! ::: ny    => input resolution
  ! ::: nz    => input resolution
  ! ::: data  <= output data
  ! ::: -----------------------------------------------------------
  subroutine read_binary(iname,nx,ny,nz,data)

    implicit none

    character(len=255), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny, nz
    double precision, intent(out) :: data(0:nx*ny*nz*6-1)

    integer, parameter :: in_unit=1
    integer :: ios = 0
    integer(kind=8) :: chunk
    integer(kind=8) :: i
    integer(kind=8), parameter :: one = 1
    integer(kind=8), parameter :: six = 6

    open(in_unit,file=trim(iname), access='stream', form='unformatted', status='old', action='read')

    ! We have to read the file in chunks in case it is bigger than (2^31-1) bytes
    chunk = nx*ny*six
    do i = 0, nz-one
       read(in_unit,iostat=ios)data(i*chunk:(i+one)*chunk-one)
    enddo
    close(in_unit)

    if (ios .ne. 0) then
       write(*,*)'Error in binary input file read. Exiting with read error', ios
       stop 99
    endif

  end subroutine read_binary


  ! ::: -----------------------------------------------------------
  ! ::: Read a csv file
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: iname => filename
  ! ::: nx    => input resolution
  ! ::: ny    => input resolution
  ! ::: nz    => input resolution
  ! ::: data  <= output data
  ! ::: -----------------------------------------------------------
  subroutine read_csv(iname,nx,ny,nz,data)

    implicit none

    character(len=255), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny, nz
    double precision, intent(out) :: data(0:nx*ny*nz*6-1)

    integer :: i
    integer, parameter :: in_unit=1
    integer :: nlines = 0, ios = 0

    ! Get number of lines in file
    open(in_unit,file=trim(iname), access='sequential', form='formatted', status='old', action='read')
    read(in_unit,*) ! skip header
    nlines = 0
    do
       read(in_unit,*,iostat=ios)
       if (ios .ne. 0) exit
       nlines = nlines + 1
    enddo
    ios = 0

    ! Quick sanity check
    if (nlines .ne. nx*ny*nz) then
       write(*,'("Number of lines in the input file (=",I0,") does not ")')nlines
       write(*,'("  match the input resolution (n=",I0,") in the probin file")')nx
       stop 99
    endif

    ! Read the data from the file
    rewind(in_unit)
    read(in_unit,*) ! skip header
    do i = 0, nlines-1
       read(in_unit, *, iostat=ios)data(i*6:(i+1)*6-1)
       if (ios .ne. 0) then
          write(*,*)'Error in CSV input file read. Exiting with read error', ios
          stop 99
       endif
    enddo
    close(in_unit)

  end subroutine read_csv


  ! ::: -----------------------------------------------------------
  ! ::: Search for the closest index in an array to a given value
  ! ::: using the bisection technique.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: xtable(0:n-1) => array to search in (ascending order)
  ! ::: n             => number of elements in array
  ! ::: x             => x location
  ! ::: idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
  ! ::: -----------------------------------------------------------
  subroutine locate(xtable, n, x, idxlo)

    implicit none

    double precision, intent(in) :: xtable(0:n-1)
    integer, intent(in) :: n
    double precision, intent(in) :: x
    integer, intent(out) :: idxlo

    ! Local variables
    integer :: idxhi, idxmid
    logical :: notdone

    ! If x is out of bounds, return boundary index
    if (x >= xtable(n-1)) then
       idxlo=n-1
       return
    elseif (x <= xtable(0)) then
       idxlo=0
       return
    endif

    ! Make sure the search array is increasing
    if (xtable(0) > xtable(n-1)) then
       write(*,'("Error in locate: non ascending input search array.")')
       stop 99
    endif

    ! Do the bisection
    idxlo = 0
    idxhi = n-1
    notdone = .true.
    do while (notdone)
       if (idxhi-idxlo <= 1) then
          notdone = .false.
       else
          idxmid = (idxhi+idxlo)/2
          if (x >= xtable(idxmid)) then
             idxlo = idxmid
          else
             idxhi = idxmid
          endif
       endif
    enddo
    return
  end subroutine locate


  ! ::: -----------------------------------------------------------
  ! ::: Calculate the Taylor length scale in x-direction.
  ! :::    lambda = <u^2>/<(du/dx)^2>
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: u(n,n,n) => input velocity array
  ! ::: x(n,n,n) => input coordinate array
  ! ::: n        => size of array
  ! ::: lambda  <=> lambda
  ! ::: -----------------------------------------------------------
  subroutine calculate_taylor_microscale(u, x, n, lambda)

    implicit none

    double precision, intent(in) :: u(n,n,n)
    double precision, intent(in) :: x(n,n,n)
    integer, intent(in) :: n
    double precision, intent(out) :: lambda

    ! Local variables
    integer :: i
    double precision :: dudx2 = 0, u2 = 0

    ! Square of velocity
    u2 = sum(u(:,:,:)*u(:,:,:))

    ! Calculate the gradients. Assume periodicity for the edges.
    do i = 2, n-1
       dudx2 = dudx2 + sum((u(i+1,:,:) - u(i-1,:,:))**2 / (x(i+1,:,:) - x(i-1,:,:))**2)
    enddo
    dudx2 = dudx2 + sum((u(2,:,:) - u(n,:,:))**2 / (x(2,:,:) - x(n,:,:))**2)
    dudx2 = dudx2 + sum((u(1,:,:) - u(n-1,:,:))**2 / (x(1,:,:) - x(n-1,:,:))**2)

    lambda = sqrt(u2/dudx2)

    return
  end subroutine calculate_taylor_microscale

end module pc_prob_module
