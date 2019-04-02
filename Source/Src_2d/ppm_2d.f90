module ppm_module

  implicit none

  private 

  public ppm

contains
  
      ! characteristics based on u
      subroutine ppm(s,s_l1,s_l2,s_h1,s_h2, &
                     u,cspd,qd_l1,qd_l2,qd_h1,qd_h2, &
                     flatn, &
                     Ip,Im, &
                     ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

        use meth_params_module, only : ppm_type, weno_variant
        use amrex_constants_module
        use weno_module

        implicit none

        integer          s_l1,s_l2,s_h1,s_h2
        integer          qd_l1,qd_l2,qd_h1,qd_h2
        integer          ilo1,ilo2,ihi1,ihi2
        double precision s(s_l1:s_h1,s_l2:s_h2)
        double precision     u(qd_l1:qd_h1,qd_l2:qd_h2,1:2)
        double precision  cspd(qd_l1:qd_h1,qd_l2:qd_h2)
        double precision flatn(s_l1:s_h1,s_l2:s_h2)
        double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3)
        double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3)
        double precision dx,dy,dt

        ! local
        integer i,j

        logical extremum, bigp, bigm

        double precision dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
        double precision sgn, sigma, s6, amax, delam, delap
        double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
        double precision dachkm, dachkp
        double precision vl, vr

        ! s_{\ib,+}, s_{\ib,-}
        double precision, allocatable :: sp(:,:)
        double precision, allocatable :: sm(:,:)

        ! \delta s_{\ib}^{vL}
        double precision, allocatable :: dsvl(:,:)

        ! s_{i+\half}^{H.O.}
        double precision, allocatable :: sedge(:,:)
        
        double precision, allocatable :: sedge_weno_p(:,:)
        double precision, allocatable :: sedge_weno_m(:,:)

        ! cell-centered indexing
        allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
        allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

        ! constant used in Colella 2008
        C = 1.25d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! x-direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! cell-centered indexing w/extra x-ghost cell
        allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

        ! edge-centered indexing for x-faces
        if (ppm_type .eq. 1) then
           allocate(sedge(ilo1-1:ihi1+2,ilo2-1:ihi2+1))
        else
           allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))
        end if
        
        if (ppm_type .eq. 3) then
          allocate(sedge_weno_p(ilo1-2:ihi1+2,ilo2-2:ihi2+2))
          allocate(sedge_weno_m(ilo1-2:ihi1+2,ilo2-2:ihi2+2))
        end if
        
        
        ! compute s at x-edges
        if (ppm_type .eq. 1) then

           ! compute van Leer slopes in x-direction (CW Eq. 1.7, 1.8
           ! w/ zone widths (dxi) all equal)
           dsvl = ZERO
           do j=ilo2-1,ihi2+1
              do i=ilo1-2,ihi1+2
                 dsc = HALF * (s(i+1,j) - s(i-1,j))
                 dsl = TWO  * (s(i  ,j) - s(i-1,j))
                 dsr = TWO  * (s(i+1,j) - s(i  ,j))
                 if (dsl*dsr .gt. ZERO) &
                      dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
              end do
           end do

           ! interpolate s to x-edges (CW 1.6)
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+2
                 sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
                 ! make sure sedge lies in between adjacent
                 ! cell-centered values -- this is not part of the
                 ! original CW algorithm, but Colella & Sekora say it
                 ! is automatically imposed by the van Leer limiting
                 ! above
                 sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i-1,j)))
                 sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i-1,j)))
              end do
           end do

           ! copy sedge into sp and sm
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sp(i,j) = sedge(i+1,j)
                 sm(i,j) = sedge(i  ,j)
              end do
           end do

           ! flatten the parabola BEFORE doing the other
           ! monotonization -- this is the method that Flash does
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sm(i,j) = flatn(i,j)*sm(i,j) + (ONE-flatn(i,j))*s(i,j)
                 sp(i,j) = flatn(i,j)*sp(i,j) + (ONE-flatn(i,j))*s(i,j)
              enddo
           enddo


           ! modify using quadratic limiters (CW 1.10) -- with a
           ! slightly different form from Colella & Sekora (Eqs. 14,
           ! 15)
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                    sp(i,j) = s(i,j)
                    sm(i,j) = s(i,j)
                 else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                    sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                 else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                    sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                 end if
              end do
           end do


        else if (ppm_type .eq. 2) then

           ! interpolate s to x-edges
           do j=ilo2-1,ihi2+1
              do i=ilo1-2,ihi1+3
                 sedge(i,j) = SEVEN12TH*(s(i-1,j)+s(i,j)) - TWELFTH*(s(i-2,j)+s(i+1,j))
                 ! limit sedge
                 if ((sedge(i,j)-s(i-1,j))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                    D2  = THREE*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                    D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                    D2R = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                    sgn = sign(ONE,D2)
                    D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                    sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - SIXTH*D2LIM
                 end if
              end do
           end do

           ! use Colella 2008 limiters
           ! This is a new version of the algorithm 
           ! to eliminate sensitivity to roundoff.
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1

                 alphap = sedge(i+1,j)-s(i,j)
                 alpham = sedge(i  ,j)-s(i,j)
                 bigp = abs(alphap).gt.TWO*abs(alpham)
                 bigm = abs(alpham).gt.TWO*abs(alphap)
                 extremum = .false.

                 if (alpham*alphap .ge. ZERO) then
                    extremum = .true.
                 else if (bigp .or. bigm) then
                    ! Possible extremum. We look at cell centered values and face
                    ! centered values for a change in sign in the differences adjacent to
                    ! the cell. We use the pair of differences whose minimum magnitude is the
                    ! largest, and thus least susceptible to sensitivity to roundoff.
                    dafacem = sedge(i,j) - sedge(i-1,j)
                    dafacep = sedge(i+2,j) - sedge(i+1,j)
                    dabarm = s(i,j) - s(i-1,j)
                    dabarp = s(i+1,j) - s(i,j)
                    dafacemin = min(abs(dafacem),abs(dafacep))
                    dabarmin= min(abs(dabarm),abs(dabarp))
                    if (dafacemin.ge.dabarmin) then
                       dachkm = dafacem
                       dachkp = dafacep
                    else
                       dachkm = dabarm
                       dachkp = dabarp
                    endif
                    extremum = (dachkm*dachkp .le. ZERO)
                 end if

                 if (extremum) then
                    D2  = SIX*(alpham + alphap)
                    D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                    D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                    D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                    sgn = sign(ONE,D2)
                    D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                    alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                    alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                 else
                    if (bigp) then
                       sgn = sign(ONE,alpham)
                       amax = -alphap**2 / (4*(alpham + alphap))
                       delam = s(i-1,j) - s(i,j)
                       if (sgn*amax .ge. sgn*delam) then
                          if (sgn*(delam - alpham).ge.1.d-10) then
                             alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                          else 
                             alphap = -TWO*alpham
                          endif
                       endif
                    end if
                    if (bigm) then
                       sgn = sign(ONE,alphap)
                       amax = -alpham**2 / (4*(alpham + alphap))
                       delap = s(i+1,j) - s(i,j)
                       if (sgn*amax .ge. sgn*delap) then
                          if (sgn*(delap - alphap).ge.1.d-10) then
                             alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                          else
                             alpham = -TWO*alphap
                          endif
                       endif
                    end if
                 end if

                 sm(i,j) = s(i,j) + alpham
                 sp(i,j) = s(i,j) + alphap

              end do
           end do

           ! flatten the parabola AFTER doing the monotonization --
           ! (ppm_type = 2 is here)
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sm(i,j) = flatn(i,j)*sm(i,j) + (ONE-flatn(i,j))*s(i,j)
                 sp(i,j) = flatn(i,j)*sp(i,j) + (ONE-flatn(i,j))*s(i,j)
              enddo
           enddo

        else if (ppm_type .eq. 3) then
    
          do j=ilo2-2,ihi2+2
            do i=ilo1-2,ihi1+2

              select case (weno_variant)
              case (0)
                call weno5js_face(s(i-3:i+2,j), vl, vr)
              case (1)
                call weno5z_face(s(i-3:i+2,j), vl, vr)
              case (2)
                call weno3z_face(s(i-3:i+2,j), vl, vr)
              case (3)
                call weno7z_face(s(i-4:i+3,j), vl, vr)
              end select     
              sedge_weno_p(i,j) = vr
              sedge_weno_m(i,j) = vl
            enddo
          enddo
       
          do j=ilo2-1,ihi2+1
            do i=ilo1-1,ihi1+1
              sm(i,j) = sedge_weno_p(i,j)
              sp(i,j) = sedge_weno_m(i+1,j) 
            enddo
          enddo
                  
        end if

       ! compute x-component of Ip and Im
        do j=ilo2-1,ihi2+1
           do i=ilo1-1,ihi1+1
              ! Ip/m is the integral under the parabola for the extent
              ! that a wave can travel over a timestep
              !
              ! Ip integrates to the right edge of a cell
              ! Im integrates to the left edge of a cell

              s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))

              ! u-c wave
              sigma = abs(u(i,j,1)-cspd(i,j))*dt/dx

              if (u(i,j,1)-cspd(i,j) <= ZERO) then
                 Ip(i,j,1,1) = sp(i,j)
              else
                 Ip(i,j,1,1) = sp(i,j) - &
                      HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,1)-cspd(i,j) >= ZERO) then
                 Im(i,j,1,1) = sm(i,j)
              else
                 Im(i,j,1,1) = sm(i,j) + &
                      HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif

              ! u wave
              sigma = abs(u(i,j,1))*dt/dx

              if (u(i,j,1) <= ZERO) then
                 Ip(i,j,1,2) = sp(i,j)
              else
                 Ip(i,j,1,2) = sp(i,j) - &
                      HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,1) >= ZERO) then
                 Im(i,j,1,2) = sm(i,j)
              else
                 Im(i,j,1,2) = sm(i,j) + &
                      HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif

              ! u + c wave
              sigma = abs(u(i,j,1)+cspd(i,j))*dt/dx

              if (u(i,j,1) + cspd(i,j) <= ZERO) then
                 Ip(i,j,1,3) = sp(i,j)
              else
                 Ip(i,j,1,3) = sp(i,j) - &
                      HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,1) + cspd(i,j) >= ZERO) then
                 Im(i,j,1,3) = sm(i,j)
              else
                 Im(i,j,1,3) = sm(i,j) + &
                      HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif
           end do
        end do

        deallocate(sedge,dsvl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! y-direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! cell-centered indexing w/extra y-ghost cell
        allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

        ! edge-centered indexing for y-faces
        if (ppm_type .eq. 1) then
           allocate(sedge(ilo1-1:ihi1+1,ilo2-1:ihi2+2))
        else
           allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))
        end if

        ! compute s at y-edges
        if (ppm_type .eq. 1) then

           ! compute van Leer slopes in y-direction
           dsvl = ZERO
           do j=ilo2-2,ihi2+2
              do i=ilo1-1,ihi1+1
                 dsc = HALF * (s(i,j+1) - s(i,j-1))
                 dsl = TWO  * (s(i,j  ) - s(i,j-1))
                 dsr = TWO  * (s(i,j+1) - s(i,j  ))
                 if (dsl*dsr .gt. ZERO) &
                      dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
              end do
           end do

           ! interpolate s to y-edges
           do j=ilo2-1,ihi2+2
              do i=ilo1-1,ihi1+1
                 sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
                 ! make sure sedge lies in between adjacent cell-centered values
                 sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i,j-1)))
                 sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i,j-1)))
              end do
           end do

           ! copy sedge into sp and sm
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sp(i,j) = sedge(i,j+1)
                 sm(i,j) = sedge(i,j  )
              end do
           end do

           ! flatten the parabola BEFORE doing the other
           ! monotonization -- this is the method that Flash does
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sm(i,j) = flatn(i,j)*sm(i,j) + (ONE-flatn(i,j))*s(i,j)
                 sp(i,j) = flatn(i,j)*sp(i,j) + (ONE-flatn(i,j))*s(i,j)
              enddo
           enddo

           ! modify using quadratic limiters
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                    sp(i,j) = s(i,j)
                    sm(i,j) = s(i,j)
                 else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                    sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                 else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                    sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                 end if
              end do
           end do

        else if (ppm_type .eq. 2) then

           ! interpolate s to y-edges
           do j=ilo2-2,ihi2+3
              do i=ilo1-1,ihi1+1
                 sedge(i,j) = SEVEN12TH*(s(i,j-1)+s(i,j)) - TWELFTH*(s(i,j-2)+s(i,j+1))
                 ! limit sedge
                 if ((sedge(i,j)-s(i,j-1))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                    D2  = THREE*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                    D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                    D2R = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                    sgn = sign(ONE,D2)
                    D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                    sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - SIXTH*D2LIM
                 end if
              end do
           end do

           ! use Colella 2008 limiters
           ! This is a new version of the algorithm 
           ! to eliminate sensitivity to roundoff.
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1

                 alphap = sedge(i,j+1)-s(i,j)
                 alpham = sedge(i,j  )-s(i,j)
                 bigp = abs(alphap).gt.TWO*abs(alpham)
                 bigm = abs(alpham).gt.TWO*abs(alphap)
                 extremum = .false.

                 if (alpham*alphap .ge. ZERO) then
                    extremum = .true.
                 else if (bigp .or. bigm) then
                    ! Possible extremum. We look at cell centered values and face
                    ! centered values for a change in sign in the differences adjacent to
                    ! the cell. We use the pair of differences whose minimum magnitude is the
                    ! largest, and thus least susceptible to sensitivity to roundoff.
                    dafacem = sedge(i,j) - sedge(i,j-1)
                    dafacep = sedge(i,j+2) - sedge(i,j+1)
                    dabarm = s(i,j) - s(i,j-1)
                    dabarp = s(i,j+1) - s(i,j)
                    dafacemin = min(abs(dafacem),abs(dafacep))
                    dabarmin= min(abs(dabarm),abs(dabarp))
                    if (dafacemin.ge.dabarmin) then
                       dachkm = dafacem
                       dachkp = dafacep
                    else
                       dachkm = dabarm
                       dachkp = dabarp
                    endif
                    extremum = (dachkm*dachkp .le. ZERO)
                 end if

                 if (extremum) then
                    D2  = SIX*(alpham + alphap)
                    D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                    D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                    D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                    sgn = sign(ONE,D2)
                    D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                    alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                    alphap = alphap*D2LIM/max(abs(D2),1.d-10)
                 else
                    if (bigp) then
                       sgn = sign(ONE,alpham)
                       amax = -alphap**2 / (4*(alpham + alphap))
                       delam = s(i,j-1) - s(i,j)
                       if (sgn*amax .ge. sgn*delam) then
                          if (sgn*(delam - alpham).ge.1.d-10) then
                             alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                          else 
                             alphap = -TWO*alpham
                          endif
                       endif
                    end if
                    if (bigm) then
                       sgn = sign(ONE,alphap)
                       amax = -alpham**2 / (4*(alpham + alphap))
                       delap = s(i,j+1) - s(i,j)
                       if (sgn*amax .ge. sgn*delap) then
                          if (sgn*(delap - alphap).ge.1.d-10) then
                             alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                          else
                             alpham = -TWO*alphap
                          endif
                       endif
                    end if
                 end if

                 sm(i,j) = s(i,j) + alpham
                 sp(i,j) = s(i,j) + alphap

              end do
           end do

           ! flatten the parabola AFTER doing the monotonization --
           ! (ppm_type = 2 is here)
           do j=ilo2-1,ihi2+1
              do i=ilo1-1,ihi1+1
                 sm(i,j) = flatn(i,j)*sm(i,j) + (ONE-flatn(i,j))*s(i,j)
                 sp(i,j) = flatn(i,j)*sp(i,j) + (ONE-flatn(i,j))*s(i,j)
              enddo
           enddo

        else if (ppm_type .eq. 3) then
    
          do j=ilo2-2,ihi2+2
            do i=ilo1-2,ihi1+2

              select case (weno_variant)
              case (0)
                call weno5js_face(s(i,j-3:j+2), vl, vr)
              case (1)
                call weno5z_face(s(i,j-3:j+2), vl, vr)
              case (2)
                call weno3z_face(s(i,j-3:j+2), vl, vr)
              case (3)
                call weno7z_face(s(i,j-4:j+3), vl, vr)
              end select     
              sedge_weno_p(i,j) = vr
              sedge_weno_m(i,j) = vl
            enddo
          enddo
       
          do j=ilo2-1,ihi2+1
            do i=ilo1-1,ihi1+1
              sm(i,j) = sedge_weno_p(i,j)
              sp(i,j) = sedge_weno_m(i,j+1) 
            enddo
          enddo
        
        end if

        ! compute y-component of Ip and Im
        do j=ilo2-1,ihi2+1
           do i=ilo1-1,ihi1+1

              s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))

              ! v-c wave
              sigma = abs(u(i,j,2)-cspd(i,j))*dt/dy

              if (u(i,j,2)-cspd(i,j) <= ZERO) then
                 Ip(i,j,2,1) = sp(i,j)
              else
                 Ip(i,j,2,1) = sp(i,j) - &
                   HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,2)-cspd(i,j) >= ZERO) then
                 Im(i,j,2,1) = sm(i,j) 
              else
                 Im(i,j,2,1) = sm(i,j) + &
                   HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif

              ! v wave
              sigma = abs(u(i,j,2))*dt/dy

              if (u(i,j,2) <= ZERO) then
                 Ip(i,j,2,2) = sp(i,j) 
              else
                 Ip(i,j,2,2) = sp(i,j) - &
                   HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,2) >= ZERO) then
                 Im(i,j,2,2) = sm(i,j) 
              else
                 Im(i,j,2,2) = sm(i,j) + &
                   HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif

              ! v+c wave
              sigma = abs(u(i,j,2)+cspd(i,j))*dt/dy

              if (u(i,j,2)+cspd(i,j) <= ZERO) then
                 Ip(i,j,2,3) = sp(i,j) 
              else
                 Ip(i,j,2,3) = sp(i,j) - &
                   HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
              endif

              if (u(i,j,2)+cspd(i,j) >= ZERO) then
                 Im(i,j,2,3) = sm(i,j) 
              else
                 Im(i,j,2,3) = sm(i,j) + &
                   HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
              endif

           end do
        end do

        deallocate(sp,sm,dsvl,sedge)
        if (ppm_type .eq. 3) then
         deallocate(sedge_weno_m)
         deallocate(sedge_weno_p)
       end if
       
        
      end subroutine ppm

end module ppm_module
