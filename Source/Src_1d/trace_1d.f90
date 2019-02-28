module trace_module

  implicit none

  private

  public trace

contains

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::
  subroutine trace(q,dq,c,flatn,qd_l1,qd_h1, &
       dloga,dloga_l1,dloga_h1, &
       srcQ,src_l1,src_h1,&
       qxm,qxp,qpd_l1,qpd_h1, &
       ilo,ihi,domlo,domhi,dx,dt)

    use meth_params_module, only : plm_iorder, QVAR, QRHO, QU, QREINT, QPRES, &
         npassive, qpass_map, small_dens, ppm_type, fix_mass_flux, use_pslope
    use prob_params_module, only : physbc_lo, physbc_hi
    use slope_module, only : uslope, pslope
    use amrex_constants_module

    implicit none

    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer    qd_l1,   qd_h1
    integer dloga_l1,dloga_h1
    integer   qpd_l1,  qpd_h1
    integer   src_l1,  src_h1
    double precision dx, dt
    double precision     q( qd_l1: qd_h1,QVAR)
    double precision  srcQ(src_l1:src_h1,QVAR)
    double precision flatn(qd_l1:qd_h1)
    double precision     c(qd_l1:qd_h1)
    double precision dloga(dloga_l1:dloga_h1)

    double precision   dq( qpd_l1: qpd_h1,QVAR)
    double precision  qxm( qpd_l1: qpd_h1,QVAR)
    double precision  qxp( qpd_l1: qpd_h1,QVAR)

    !     Local variables
    integer          :: i
    integer          :: n, ipassive

    double precision :: hdt,dtdx
    double precision :: cc, csq, rho, u, p, rhoe
    double precision :: drho, du, dp, drhoe

    double precision :: enth, alpham, alphap, alpha0r, alpha0e
    double precision :: spminus, spplus, spzero
    double precision :: apright, amright, azrright, azeright
    double precision :: apleft, amleft, azrleft, azeleft
    double precision :: acmprght, acmpleft
    double precision :: sourcr,sourcp,source,courn,eta,dlogatmp

    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in trace with ppm_type != 0'
       call bl_error("Error:: PeleC_1d.f90 :: trace")
    end if

    hdt = HALF * dt
    dtdx = dt/dx

    ! Compute slopes
    if (plm_iorder .eq. 1) then

       dq(ilo-1:ihi+1,1:QVAR) = ZERO

    else

       call uslope(q,flatn,qd_l1,qd_h1, &
            dq,qpd_l1,qpd_h1, &
            ilo,ihi,QVAR)

       if (use_pslope .eq. 1) &
            call pslope(q(:,QPRES),q(:,QRHO), &
            flatn      , qd_l1, qd_h1, &
            dq(:,QPRES),qpd_l1,qpd_h1, &
            srcQ       ,src_l1,src_h1, &
            ilo,ihi,dx)

    endif

    ! Compute left and right traced states
    do i = ilo-1, ihi+1

       cc = c(i)
       csq = cc**2
       rho = q(i,QRHO)
       u = q(i,QU)
       p = q(i,QPRES)
       rhoe = q(i,QREINT)
       enth = ( (rhoe+p)/rho )/csq

       drho  = dq(i,QRHO)
       du    = dq(i,QU)
       dp    = dq(i,QPRES)
       drhoe = dq(i,QREINT)

       alpham = HALF*(dp/(rho*cc) - du)*rho/cc
       alphap = HALF*(dp/(rho*cc) + du)*rho/cc
       alpha0r = drho - dp/csq
       alpha0e = drhoe - dp*enth

       if (u-cc .gt. ZERO) then
          spminus = -ONE
       else
          spminus = (u-cc)*dtdx
       endif
       if (u+cc .gt. ZERO) then
          spplus = -ONE
       else
          spplus = (u+cc)*dtdx
       endif
       if (u .gt. ZERO) then
          spzero = -ONE
       else
          spzero = u*dtdx
       endif

       apright = HALF*(-ONE - spplus )*alphap
       amright = HALF*(-ONE - spminus)*alpham
       azrright = HALF*(-ONE - spzero )*alpha0r
       azeright = HALF*(-ONE - spzero )*alpha0e

       if (i .ge. ilo) then
          qxp(i,QRHO  ) = rho + apright + amright + azrright
          qxp(i,QRHO  ) = max(small_dens,qxp(i,QRHO))
          qxp(i,QU    ) = u + (apright - amright)*cc/rho
          qxp(i,QPRES ) = p + (apright + amright)*csq
          qxp(i,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

          ! add source term
          qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
          qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
          qxp(i  ,QU    ) = qxp(i,QU    ) + hdt*srcQ(i,QU)
          qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
          qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)
       end if

       if (u-cc .ge. ZERO) then
          spminus = (u-cc)*dtdx
       else
          spminus = ONE
       endif
       if (u+cc .ge. ZERO) then
          spplus = (u+cc)*dtdx
       else
          spplus = ONE
       endif
       if (u .ge. ZERO) then
          spzero = u*dtdx
       else
          spzero = ONE
       endif

       apleft = HALF*(ONE - spplus )*alphap
       amleft = HALF*(ONE - spminus)*alpham
       azrleft = HALF*(ONE - spzero )*alpha0r
       azeleft = HALF*(ONE - spzero )*alpha0e

       if (i .le. ihi) then
          qxm(i+1,QRHO  ) = rho + apleft + amleft + azrleft
          qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QU    ) = u + (apleft - amleft)*cc/rho
          qxm(i+1,QPRES ) = p + (apleft + amleft)*csq
          qxm(i+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

          ! add source terms
          qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
          qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QU    ) = qxm(i+1,QU    ) + hdt*srcQ(i,QU)
          qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
          qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)
       end if

       if(dloga(i).ne.0)then
          courn = dtdx*(cc+abs(u))
          eta = (ONE-courn)/(cc*dt*abs(dloga(i)))
          dlogatmp = min(eta,ONE)*dloga(i)
          sourcr = -HALF*dt*rho*dlogatmp*u
          sourcp = sourcr*csq
          source = sourcp*enth
          if (i .le. ihi) then
             qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
             qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
             qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
             qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
          end if
          if (i .ge. ilo) then
             qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
             qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
             qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
             qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
          end if
       endif
    enddo

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Right state
       do i = ilo,ihi+1
          u = q(i,QU)
          if (u .gt. ZERO) then
             spzero = -ONE
          else
             spzero = u*dtdx
          endif
          acmprght = HALF*(-ONE - spzero )*dq(i,n)
          qxp(i,n) = q(i,n) + acmprght
       enddo

       ! Left state
       do i = ilo-1,ihi
          u = q(i,QU)
          if (u .ge. ZERO) then
             spzero = u*dtdx
          else
             spzero = ONE
          endif
          acmpleft = HALF*(ONE - spzero )*dq(i,n)
          qxm(i+1,n) = q(i,n) + acmpleft
       enddo
    enddo

  end subroutine trace

end module trace_module
