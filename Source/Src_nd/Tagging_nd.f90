module tagging_module

  implicit none

  double precision, save ::    denerr,   dengrad
  double precision, save ::    enterr,   entgrad
  double precision, save ::    velerr,   velgrad
  double precision, save ::    vorterr,  vfracerr
  double precision, save ::   temperr,  tempgrad
  double precision, save ::  presserr, pressgrad
  double precision, save ::    raderr,   radgrad
  double precision, save ::  ftracerr, ftracgrad
  integer         , save ::  max_denerr_lev,   max_dengrad_lev
  integer         , save ::  max_enterr_lev,   max_entgrad_lev
  integer         , save ::  max_velerr_lev,   max_velgrad_lev
  integer         , save ::  max_vorterr_lev
  integer         , save ::  max_vfracerr_lev
  integer         , save ::  max_temperr_lev,  max_tempgrad_lev
  integer         , save ::  max_presserr_lev, max_pressgrad_lev
  integer         , save ::  max_raderr_lev,   max_radgrad_lev
  integer         , save ::  max_ftracerr_lev, max_ftracgrad_lev

  public

contains

  ! All tagging subroutines in this file must be threadsafe because
  ! they are called inside OpenMP parallel regions.

  ! ::: -----------------------------------------------------------
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: tag      <=  integer tag array
  ! ::: lo,hi     => index extent of work region
  ! ::: set       => integer value to tag cell for refinement
  ! ::: clear     => integer value to untag cell
  ! ::: temp      => temperature array
  ! ::: np        => number of components in temp array (should be 1)
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::              corner of work region
  ! ::: problo    => phys loc of lower left corner of prob domain
  ! ::: time      => problem evolution time
  ! ::: level     => refinement level of this array
  ! ::: -----------------------------------------------------------

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the density
  ! ::: -----------------------------------------------------------

  subroutine pc_denerror(tag,taglo,taghi, &
                         set,clear, &
                         den,denlo,denhi, &
                         lo,hi,nd,domlo,domhi, &
                         delta,xlo,problo,time,level) &
                         bind(C, name="pc_denerror")

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nd, level
    integer          :: taglo(3), taghi(3)
    integer          :: denlo(3), denhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3),nd)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high density
    if (level .lt. max_denerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (den(i,j,k,1) .ge. denerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high density gradient
    if (level .lt. max_dengrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(den(i+1*dg(1),j,k,1) - den(i,j,k,1))
                ay = ABS(den(i,j+1*dg(2),k,1) - den(i,j,k,1))
                az = ABS(den(i,j,k+1*dg(3),1) - den(i,j,k,1))
                ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. dengrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine pc_denerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the temperature
  ! ::: -----------------------------------------------------------

  subroutine pc_temperror(tag,taglo,taghi, &
                          set,clear, &
                          temp,templo,temphi, &
                          lo,hi,np,domlo,domhi, &
                          delta,xlo,problo,time,level) &
                          bind(C, name="pc_temperror")

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: templo(3), temphi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (temp(i,j,k,1) .ge. temperr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(temp(i+1*dg(1),j,k,1) - temp(i,j,k,1))
                ay = ABS(temp(i,j+1*dg(2),k,1) - temp(i,j,k,1))
                az = ABS(temp(i,j,k+1*dg(3),1) - temp(i,j,k,1))
                ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. tempgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine pc_temperror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on a flame trac
  ! ::: -----------------------------------------------------------

  subroutine pc_ftracerror(tag,taglo,taghi, &
                           set,clear, &
                           ftrac,ftraclo,ftrachi, &
                           lo,hi,np,domlo,domhi, &
                           delta,xlo,problo,time,level) &
                           bind(C, name="pc_ftracerror")

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: ftraclo(3), ftrachi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: ftrac(ftraclo(1):ftrachi(1),ftraclo(2):ftrachi(2),ftraclo(3):ftrachi(3),np)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high ftracerature
    if (level .lt. max_ftracerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ftrac(i,j,k,1) .ge. ftracerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high ftracerature gradient
    if (level .lt. max_ftracgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(ftrac(i+1*dg(1),j,k,1) - ftrac(i,j,k,1))
                ay = ABS(ftrac(i,j+1*dg(2),k,1) - ftrac(i,j,k,1))
                az = ABS(ftrac(i,j,k+1*dg(3),1) - ftrac(i,j,k,1))
                ax = MAX(ax,ABS(ftrac(i,j,k,1) - ftrac(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(ftrac(i,j,k,1) - ftrac(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(ftrac(i,j,k,1) - ftrac(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. ftracgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine pc_ftracerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the pressure
  ! ::: -----------------------------------------------------------

  subroutine pc_presserror(tag,taglo,taghi, &
                           set,clear, &
                           press,presslo,presshi, &
                           lo,hi,np,domlo,domhi, &
                           delta,xlo,problo,time,level) &
                           bind(C, name="pc_presserror")

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: presslo(3), presshi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high pressure
    if (level .lt. max_presserr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (press(i,j,k,1) .ge. presserr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high pressure gradient
    if (level .lt. max_pressgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
                ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
                az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
                ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. pressgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine pc_presserror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the velocity
  ! ::: -----------------------------------------------------------

  subroutine pc_velerror(tag,taglo,taghi, &
                         set,clear, &
                         vel,vello,velhi, &
                         lo,hi,nv,domlo,domhi, &
                         delta,xlo,problo,time,level) &
                         bind(C, name="pc_velerror")

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nv, level
    integer          :: taglo(3), taghi(3)
    integer          :: vello(3), velhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),nv)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high velocity
    if (level .lt. max_velerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ABS(vel(i,j,k,1)) .ge. velerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high velocity gradient
    if (level .lt. max_velgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(vel(i+1*dg(1),j,k,1) - vel(i,j,k,1))
                ay = ABS(vel(i,j+1*dg(2),k,1) - vel(i,j,k,1))
                az = ABS(vel(i,j,k+1*dg(3),1) - vel(i,j,k,1))
                ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. velgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine pc_velerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the vorticity
  ! ::: -----------------------------------------------------------

  subroutine pc_vorterror(tag,taglo,taghi, &
                          set,clear, &
                          vort,vortlo,vorthi, &
                          lo,hi,nv,domlo,domhi, &
                          delta,xlo,problo,time,level) &
                          bind(C, name="pc_vorterror")

    implicit none

    integer          :: set, clear, nv, level
    integer          :: taglo(3), taghi(3)
    integer          :: vortlo(3), vorthi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: vort(vortlo(1):vorthi(1),vortlo(2):vorthi(2),vortlo(3):vorthi(3),nv)
    double precision :: delta(3), xlo(3), problo(3), time

    integer          :: i, j, k

    !     Tag on regions of high vorticity
    if (level .lt. max_vorterr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ABS(vort(i,j,k,1)).ge.vorterr*2.d0**level) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    end if

  end subroutine pc_vorterror


  ! ::: -----------------------------------------------------------------
  ! ::: This routine will tag high error cells based on the volume fraction
  ! ::: -------------------------------------------------------------------

  subroutine pc_vfracerror(tag,taglo,taghi, &
                          set,clear, &
                          vfrac,vfraclo,vfrachi, &
                          lo,hi,nv,domlo,domhi, &
                          delta,xlo,problo,time,level) &
                          bind(C, name="pc_vfracerror")

    implicit none

    integer          :: set, clear, nv, level
    integer          :: taglo(3), taghi(3)
    integer          :: vfraclo(3), vfrachi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: vfrac(vfraclo(1):vfrachi(1),vfraclo(2):vfrachi(2),vfraclo(3):vfrachi(3),nv)
    double precision :: delta(3), xlo(3), problo(3), time

    integer          :: i, j, k

    !     Tag on regions of high vorticity
    if (level .lt. max_vfracerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (vfrac(i,j,k,1).gt.0.0 .and.vfrac(i,j,k,1).lt.1.0 ) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    end if

  end subroutine pc_vfracerror


end module tagging_module
