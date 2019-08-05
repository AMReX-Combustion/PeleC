module pc_prob_module

  implicit none

  private

  public :: amrex_probinit, pc_initdata, pc_prob_close

contains
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name = "amrex_probinit")
    
    use probdata_module
    use prob_params_module, only : center, dim
    use amrex_fort_module

    implicit none
    integer :: init, namlen
    integer :: name(namlen)
    double precision :: problo(3), probhi(3)

    integer :: untin,i

    namelist /fortin/ probtype, p_ambient, dens_ambient, exp_energy, &
         r_init, nsub

    ! Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character :: probin*(maxlen)

    if (namlen .gt. maxlen) then
       call bl_error('probin file name too long')
    end if

    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do

    ! set namelist defaults

    p_ambient = 1.d-5        ! ambient pressure (in erg/cc)
    dens_ambient = 1.d0      ! ambient density (in g/cc)
    exp_energy = 1.d0        ! absolute energy of the explosion (in erg)
    r_init = 0.05d0          ! initial radius of the explosion (in cm)
    nsub = 4

    ! Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! set local variable defaults
    do i=1,dim
       center(i) = (problo(i)+probhi(i))/2.d0
    enddo

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
  ! ::: nvar      => number of state components.
  ! ::: state     <= scalar array
  ! ::: delta     => cell size
  ! ::: xlo, xhi  => physical locations of lower left and upper
  ! :::              right hand corner of grid.  (does not include
  ! :::		   ghost region).
  ! ::: -----------------------------------------------------------

  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name="pc_initdata")
    use network, only: nspecies
    use probdata_module
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use prob_params_module, only : center, dim
    use amrex_constants_module, only: ZERO, M_PI, FOUR3RD, ONE, HALF
    use eos_type_module
    use eos_module
    implicit none

    integer :: level, nvar
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    double precision :: xlo(3), xhi(3), time, delta(3)
    double precision :: state(state_lo(1):state_hi(1), &
         state_lo(2):state_hi(2), &
         state_lo(3):state_hi(3),nvar)


    double precision :: xmin,ymin,zmin
    double precision :: xx, yy, zz, xl, xr
    double precision :: dist
    double precision :: vctr, p_exp, vol_ambient, vol_pert, dx_sub

    integer :: i,j,k, ii, jj, kk
    integer :: npert, nambient

    type (eos_t) :: eos_state

    call build(eos_state)

    eos_state % rho = dens_ambient
    eos_state % e = exp_energy / eos_state % rho
    eos_state % massfrac = ZERO
    eos_state % massfrac(nspecies) = ONE
    call eos_re(eos_state)

    if ( (dim.eq.3 .and. probtype.eq.32) &
         .or. (dim.eq.2 .and. probtype.eq.21) &
         .or. (dim.eq.1 .and. probtype.eq.11) ) then

       ! set explosion pressure -- we will convert the point-explosion
       ! energy into a corresponding pressure distributed throughout the
       ! perturbed volume
       vctr  = M_PI*r_init**2
       p_exp = eos_state % p / vctr

       do k = lo(3), hi(3)
          zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

          do j = lo(2), hi(2)
             ymin = xlo(2) + delta(2)*dble(j-lo(2))

             do i = lo(1), hi(1)
                xmin = xlo(1) + delta(1)*dble(i-lo(1))

                if (dim.eq.1) then

                   vol_pert    = 0.d0
                   vol_ambient = 0.d0
                   dx_sub = delta(1)/dble(nsub)

                   do ii = 0, nsub-1
                      xx = xmin + (dble(ii) + 0.5d0) * dx_sub
                      dist = xx
                      if(dist <= r_init) then
                         vol_pert = vol_pert + dist
                      else
                         vol_ambient = vol_ambient + dist
                      endif
                   enddo

                   eos_state % p = (vol_pert*p_exp + vol_ambient*p_ambient)/(vol_pert+vol_ambient)

                else

                   npert = 0
                   nambient = 0

                   do jj = 0, nsub-1
                      yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)

                      do ii = 0, nsub-1
                         xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

                         dist = (center(1)-xx)**2 + (center(2)-yy)**2

                         if(dist <= r_init**2) then
                            npert = npert + 1
                         else
                            nambient = nambient + 1
                         endif

                      enddo
                   enddo

                   eos_state % p = (dble(npert)*p_exp + dble(nambient)*p_ambient) / &
                        (dble(npert) + dble(nambient))

                endif

                eos_state % rho = dens_ambient

                call eos_rp(eos_state)

                state(i,j,k,URHO) = eos_state % rho
                state(i,j,k,UMX) = 0.d0
                state(i,j,k,UMY) = 0.d0
                state(i,j,k,UMZ) = 0.d0

                state(i,j,k,UEINT) = eos_state % rho * eos_state % e

                state(i,j,k,UEDEN) = state(i,j,k,UEINT) +  &
                     0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                     state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                     state(i,j,k,UMZ)**2/state(i,j,k,URHO))

                state(i,j,k,UFS:UFS+nspecies-1) = eos_state % massfrac(:) * state(i,j,k,URHO)

             enddo
          enddo
       enddo

    else if ( (dim.eq.1 .and. probtype.eq.12) &
         .or. (dim.eq.2 .and. probtype.eq.23) &
         .or. (dim.eq.3 .and. probtype.eq.33) ) then

       ! set explosion pressure -- we will convert the point-explosion energy into
       ! a corresponding pressure distributed throughout the perturbed volume
       vctr  = FOUR3RD*M_PI*r_init**3
       p_exp = eos_state % p / vctr

       do k = lo(3), hi(3)
          zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

          do j = lo(2), hi(2)
             ymin = xlo(2) + delta(2)*dble(j-lo(2))

             do i = lo(1), hi(1)
                xmin = xlo(1) + delta(1)*dble(i-lo(1))

                ! Should make this cleaner....contains direct copy of previous code...
                if (dim.eq.1) then

                   vol_pert    = 0.d0
                   vol_ambient = 0.d0
                   dx_sub = delta(1)/dble(nsub)

                   do ii = 0, nsub-1
                      xl = xmin + (dble(ii)        ) * dx_sub
                      xr = xmin + (dble(ii) + 1.0d0) * dx_sub
                      xx = 0.5d0*(xl + xr)

                      ! the volume of a subzone is (4/3) pi (xr^3 - xl^3).
                      ! we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
                      ! The (4/3) pi dr factor is common, so we can neglect it. 
                      if(xx <= r_init) then
                         vol_pert = vol_pert + (xr*xr + xl*xr + xl*xl)
                      else
                         vol_ambient = vol_ambient + (xr*xr + xl*xr + xl*xl)
                      endif
                   enddo

                   eos_state % p = (vol_pert*p_exp + vol_ambient*p_ambient)/(vol_pert+vol_ambient)

                   print *,i,eos_state%p



                elseif (dim.eq.2) then

                   vol_pert    = 0.d0
                   vol_ambient = 0.d0

                   do jj = 0, nsub-1
                      yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)

                      do ii = 0, nsub-1
                         xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

                         dist = sqrt(xx**2 + yy**2)

                         ! The volume of a cell is a annular cylindrical region.  
                         ! The main thing that matters is the distance from the
                         ! symmetry axis.
                         !   V = pi*dy*(x_r**2 - x_l**2) = pi*dy*dx*HALF*xx
                         ! (where x_r is the coordinate of the x right edge,
                         !        x_l is the coordinate of the x left edge,
                         !    and xx  is the coordinate of the x center of the cell)
                         !
                         ! since dx and dy are constant, they cancel out
                         if (dist <= r_init) then
                            vol_pert    = vol_pert    + xx
                         else
                            vol_ambient = vol_ambient + xx
                         endif

                      end do
                   end do

                   eos_state % p = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

                else

                   npert = 0
                   nambient = 0

                   do kk = 0, nsub-1
                      zz = zmin + (delta(3)/dble(nsub))*(kk + 0.5d0)

                      do jj = 0, nsub-1
                         yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)

                         do ii = 0, nsub-1
                            xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

                            dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2

                            if(dist <= r_init**2) then
                               npert = npert + 1
                            else
                               nambient = nambient + 1
                            endif

                         end do
                      end do
                   end do

                   eos_state % p = (dble(npert)*p_exp + dble(nambient)*p_ambient)/  &
                        dble(npert + nambient)

                end if

                eos_state % rho = dens_ambient

                call eos_rp(eos_state)

                state(i,j,k,URHO) = eos_state % rho
                state(i,j,k,UMX) = 0.d0
                state(i,j,k,UMY) = 0.d0
                state(i,j,k,UMZ) = 0.d0

                state(i,j,k,UEINT) = eos_state % rho * eos_state % e

                state(i,j,k,UEDEN) = state(i,j,k,UEINT) +  &
                     0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                     state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                     state(i,j,k,UMZ)**2/state(i,j,k,URHO))

                state(i,j,k,UFS:UFS+nspecies-1) = eos_state % massfrac(:) * state(i,j,k,URHO)

             enddo
          enddo
       enddo

    else if (dim.eq.2 .and. probtype.eq.22) then

       !  set explosion pressure -- we will convert the point-explosion
       !  energy into a corresponding pressure distributed throughout
       !  the perturbed volume
       vctr  = M_PI*r_init**2
       p_exp = eos_state % p / vctr
       k = lo(3)
       j = lo(2)

       do i = lo(1), hi(1)
          xmin = xlo(1) + delta(1)*dble(i-lo(1))

          vol_pert    = 0.d0
          vol_ambient = 0.d0

          do ii = 0, nsub-1
             xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)

             dist = xx

             if (dist <= r_init) then
                vol_pert    = vol_pert    + dist
             else
                vol_ambient = vol_ambient + dist
             endif

          enddo

          eos_state % rho = dens_ambient
          eos_state % p = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

          call eos_rp(eos_state)

          state(i,j,k,URHO) = eos_state % rho
          state(i,j,k,UMX) = 0.d0
          state(i,j,k,UMY) = 0.d0
          state(i,j,k,UMZ) = 0.d0

          state(i,j,k,UEINT) = eos_state % rho * eos_state % e

          state(i,j,k,UEDEN) = state(i,j,k,UEINT) +  &
               0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
               state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
               state(i,j,k,UMZ)**2/state(i,j,k,URHO))

          state(i,j,k,UFS:UFS+nspecies-1) = eos_state % massfrac(:) * state(i,j,k,URHO)

       enddo

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                state(i,j,k,URHO ) = state(i,lo(2),lo(3),URHO)
                state(i,j,k,UMX  ) = state(i,lo(2),lo(3),UMX)
                state(i,j,k,UMY  ) = state(i,lo(2),lo(3),UMY)
                state(i,j,k,UMZ  ) = state(i,lo(2),lo(3),UMZ)
                state(i,j,k,UEDEN) = state(i,lo(2),lo(3),UEDEN)
                state(i,j,k,UEINT) = state(i,lo(2),lo(3),UEINT)
                state(i,j,k,UFS:UFS+nspecies-1) = state(i,lo(2),lo(3),UFS:UFS+nspecies-1)
             end do
          end do
       enddo

    else 

       call bl_error('Dont know this probtype in initdata')

    end if

    call destroy(eos_state)

  end subroutine pc_initdata

  subroutine pc_prob_close() &
       bind(C, name="pc_prob_close")
  end subroutine pc_prob_close

end module pc_prob_module
