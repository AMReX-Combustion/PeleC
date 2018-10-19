  module probdata_module
  use eos_module
  implicit none

  integer, save :: TestNum, numPoints
  character(len=128) :: OutputFile

contains

  subroutine RunRStests

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, UTEMP, UFS, UFA, NADV, UFX
    use eos_type_module
    use riemann_util_module

    integer :: iN2, iO2, iH2, iH2O
    type(eos_t) :: eos_state
    type(eos_t) :: gdnv_state
    character(len=128) :: FileName
    integer :: iunit,i

    integer, parameter :: R_RHO = 1
    integer, parameter :: R_UN  = 2
    integer, parameter :: R_UT1 = 3
    integer, parameter :: R_UT2 = 4
    integer, parameter :: R_P   = 5
    integer, parameter :: R_Y   = 6

    integer, parameter :: VECLEN = 16



    ! Other left and right state arrays
    double precision :: qtempl(VECLEN,1:5+nspec)
    double precision :: qtempr(VECLEN,1:5+nspec)
    double precision :: rhoe_l(VECLEN)
    double precision :: rhoe_r(VECLEN)
    double precision :: cspeed(VECLEN)
    double precision :: gamc_l(VECLEN)
    double precision :: gamc_r(VECLEN)
    double precision :: cavg(VECLEN)
    double precision :: csmall(VECLEN)

    ! Riemann solve work arrays
    double precision, dimension(VECLEN) :: u_gd, v_gd, w_gd, &
                                           p_gd, game_gd, re_gd, &
                                           r_gd, ustar
    double precision :: flux_tmp_single(numPoints, NVAR)
    double precision :: flux_tmp_vec(numPoints, NVAR)

    integer :: vis, vie, vic, nspc, vii, nsp


    double precision :: uvel, press, vvel, wvel, Temp, pamb
    double precision :: cs, rho, rhoe, gamc
    integer, parameter :: coord_type = 0
    integer, parameter :: bc_test_val = 1
    integer, parameter :: idir = 1

    double precision, parameter :: small = 1.d-8

    call build(eos_state)
    call build(gdnv_state)



    select case(TestNum)

    case (1)

       uvel = 1.0d0
       vvel = 1.0d0
       wvel = 1.0d0
       Temp = 300.0
       pamb = 1013250.0

       iunit = 10

       eos_state % massfrac(1:nspec) = 0.0
       eos_state % massfrac(1) = 1.0d0
       eos_state % p = pamb
       eos_state % T = Temp
       call eos_tp(eos_state)
       cs = eos_state % cs
       rho = eos_state % rho
       write(*,*) "Computed rho as ", rho
       gamc = eos_state % gam1
       rhoe = eos_state % rho * eos_state%e
       flux_tmp_single = -9.73625d0

       do i  = 1, numPoints, VECLEN
             vis = i
             vie = min(vis+VECLEN-1, numPoints)
             vic = vie - vis + 1

             cspeed(1:vic      ) = cs
             qtempl(1:vic,R_UN ) = uvel
             qtempl(1:vic,R_P  ) = press
             qtempl(1:vic,R_UT1) = vvel
             qtempl(1:vic,R_UT2) = wvel
             qtempl(1:vic,R_RHO) = rho
             do nsp = 1,nspec
               qtempl(1:vic,R_Y-1+nsp) = 0.0
             enddo
             qtempl(1:vic,R_Y) = 1.0

             rhoe_l(1:vic) = rhoe
             gamc_l(1:vic) = gamc

             qtempr(1:vic,R_UN ) = -1.0*uvel
             qtempr(1:vic,R_P  ) = press
             qtempr(1:vic,R_UT1) = vvel
             qtempr(1:vic,R_UT2) = wvel
             qtempr(1:vic,R_RHO) = rho
             do nsp = 1,nspec
               qtempr(1:vic,R_Y-1+nsp) = 0.0
             enddo
             qtempr(1:vic,R_Y) = 1.0

             rhoe_r(1:vic) = rhoe
             gamc_r(1:vic) = gamc

             cavg = cspeed
             csmall = max(small, small * cspeed)



             do vii = 1, VECLEN
                call riemann_md_singlepoint( &
                     qtempl(vii,R_RHO), qtempl(vii,R_UN), qtempl(vii,R_UT1), qtempl(vii,R_UT2), qtempl(vii,R_P), rhoe_l(vii), qtempl(vii,R_Y:R_Y-1+nspec), gamc_l(vii),&
                     qtempr(vii,R_RHO), qtempr(vii,R_UN), qtempr(vii,R_UT1), qtempr(vii,R_UT2), qtempr(vii,R_P), rhoe_r(vii), qtempr(vii,R_Y:R_Y-1+nspec), gamc_r(vii),&
                     u_gd(vii), v_gd(vii), w_gd(vii), p_gd(vii), game_gd(vii), re_gd(vii), r_gd(vii), ustar(vii),&
                     eos_state, gdnv_state, nspec,&
                     flux_tmp_single(vis+vii-1,URHO), flux_tmp_single(vii+vis-1,UMX), flux_tmp_single(vii+vis-1,UMY), flux_tmp_single(vii+vis-1,UMZ), flux_tmp_single(vii+vis-1,UEDEN), flux_tmp_single(vii+vis-1,UEINT), &
                     idir, coord_type, bc_test_val, csmall(vii), cavg(vii) )
                     ! Compute species flux like passive scalar from intermediate state
                     do nsp = 0, nspec-1
                        if (ustar(vii) .gt. ZERO) then
                           flux_tmp_single(vii,UFS+nsp) = flux_tmp_single(vii,URHO)*qtempl(vii,R_Y+nsp)
                        else if (ustar(vii) .lt. ZERO) then
                           flux_tmp_single(vii,UFS+nsp) = flux_tmp_single(vii,URHO)*qtempr(vii,R_Y+nsp)
                        else
                           flux_tmp_single(vii,UFS+nsp) = flux_tmp_single(vii,URHO)&
                                *HALF*(qtempl(vii,R_Y+nsp) + qtempr(vii,R_Y+nsp) )
                        endif
                     enddo

                    ! Clear unused flux slots
                    flux_tmp_single(vii+vis-1, UTEMP) = 0.0
                    if (naux .gt. 0) then
                       flux_tmp_single(vii+vis-1, UFX:UFX+naux) = 0.0
                    endif
                    if (nadv .gt. 0) then
                       flux_tmp_single(vii+vis-1, UFA:UFA+nadv) = 0.0
                    endif
             enddo

             call riemann_md_vec( &
                  qtempl(1:vic,R_RHO), qtempl(1:vic,R_UN), qtempl(1:vic,R_UT1), qtempl(1:vic,R_UT2), qtempl(1:vic,R_P), rhoe_l(1:vic), qtempl(1:vic,R_Y:R_Y-1+nspec), gamc_l(1:vic),&
                  qtempr(1:vic,R_RHO), qtempr(1:vic,R_UN), qtempr(1:vic,R_UT1), qtempr(1:vic,R_UT2), qtempr(1:vic,R_P), rhoe_r(1:vic), qtempr(1:vic,R_Y:R_Y-1+nspec), gamc_r(1:vic),&
                  u_gd(1:vic), v_gd(1:vic), w_gd(1:vic), p_gd(1:vic), game_gd(1:vic), re_gd(1:vic), r_gd(1:vic), ustar(1:vic),&
                  eos_state, nspec,&
                  flux_tmp_vec(vis:vie,URHO), flux_tmp_vec(vis:vie,UMX), flux_tmp_vec(vis:vie,UMY), flux_tmp_vec(vis:vie,UMZ), flux_tmp_vec(vis:vie,UEDEN), flux_tmp_vec(vis:vie,UEINT), &
                  bc_test_val, csmall(1:vic), cavg(1:vic), VECLEN)
             
             ! Compute species flux like passive scalar from intermediate state
             !do nsp = 0, nspec-1
             !   if (ustar(1:vic) .gt. ZERO) then
             !      flux_tmp_vec(1:vic,UFS+nsp) = flux_tmp_vec(1:vic,URHO)*qtempl(1:vic,R_Y+nsp)
             !   else if (ustar(1:vic) .lt. ZERO) then
             !      flux_tmp_vec(1:vic,UFS+nsp) = flux_tmp_vec(1:vic,URHO)*qtempr(1:vic,R_Y+nsp)
             !   else
             !      flux_tmp_vec(1:vic,UFS+nsp) = flux_tmp_vec(1:vic,URHO)&
             !           *HALF*(qtempl(1:vic,R_Y+nsp) + qtempr(1:vic,R_Y+nsp) )
             !   endif
             !enddo
!
             !! Clear unused flux slots
             !flux_tmp_vec(1:vic, UTEMP) = 0.0
             !if (naux .gt. 0) then
             !   flux_tmp_vec(1:vic, UFX:UFX+naux) = 0.0
             !endif
             !if (nadv .gt. 0) then
             !   flux_tmp_vec(1:vic, UFA:UFA+nadv) = 0.0
             !endif



       enddo

       FileName = trim(OutputFile)//'_vec.txt'
       OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
       ! Output to verify is flux_tmp_single vs flux_tmp_vec (vs. flux_tmp_analytic
       do i = 1, numPoints
           write(iunit,'(6(E20.8,4x))') flux_tmp_vec(i,URHO), flux_tmp_vec(i,UMX), flux_tmp_vec(i,UMY),&
                                        flux_tmp_vec(i,UMZ), flux_tmp_vec(i,UEDEN), flux_tmp_vec(i,UEINT)
       enddo
       CLOSE(UNIT=iunit)

       FileName = trim(OutputFile)//'_single.txt'
       OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
       do i = 1, numPoints
           write(iunit,'(6(E20.8,4x))') flux_tmp_single(i,URHO), flux_tmp_single(i,UMX), flux_tmp_single(i,UMY), &
                                        flux_tmp_single(i,UMZ), flux_tmp_single(i,UEDEN), flux_tmp_single(i,UEINT)
       enddo
       CLOSE(UNIT=iunit)

       FileName = trim(OutputFile)//'_diff.txt'
       OPEN(UNIT=iunit,FILE=trim(FileName),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
       do i = 1, numPoints
           write(iunit,'(6(E20.8,4x))') flux_tmp_single(i,URHO) -   flux_tmp_vec(i,URHO), &
                                        flux_tmp_single(i,UMX)  -   flux_tmp_vec(i,UMX),  &
                                        flux_tmp_single(i,UMY)  -   flux_tmp_vec(i,UMY),   &
                                        flux_tmp_single(i,UMZ)  -   flux_tmp_vec(i,UMZ), &
                                        flux_tmp_single(i,UEDEN)-   flux_tmp_vec(i,UEDEN), &
                                        flux_tmp_single(i,UEINT)-   flux_tmp_vec(i,UEINT)
       enddo
       CLOSE(UNIT=iunit)








    end select
    
    call destroy(eos_state)
 
  end subroutine RunRStests


end module probdata_module
