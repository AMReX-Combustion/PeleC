! problem-specific Fortran derive routines go here
module problem_derive_module

  implicit none

  public

contains

  subroutine pc_derrhommserror(rhommserror,m_lo,m_hi,nm, &
                               dat,d_lo,d_hi,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no) &
                               bind(C, name="pc_derrhommserror")
    !
    ! This routine will calculate the error between the density
    ! solution and the MMS solution
    !

    use meth_params_module, only : URHO

    use amrex_constants_module
#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: m_lo(3), m_hi(3), nm
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: rhommserror(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nm)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z,rho

#ifdef USE_MASA
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             ! Get error
             rho = masa_eval_3d_exact_rho(x,y,z)
             rhommserror(i,j,k,1) = dat(i,j,k,URHO) - rho

          end do
       end do
    end do

#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

  end subroutine pc_derrhommserror


  subroutine pc_derummserror(ummserror,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_derummserror")
    !
    ! This routine will calculate the error between the velocity
    ! solution and the MMS solution
    !
    use meth_params_module, only : URHO, UMX
    use amrex_constants_module
#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: ummserror(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z, u, rhoinv

#ifdef USE_MASA
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             ! Get error
             u = masa_eval_3d_exact_u(x,y,z)
             rhoinv = 1.d0/dat(i,j,k,URHO)
             ummserror(i,j,k,1) = dat(i,j,k,UMX) * rhoinv - u

          end do
       end do
    end do
#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

  end subroutine pc_derummserror

  subroutine pc_dervmmserror(vmmserror,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_dervmmserror")
    !
    ! This routine will calculate the error between the velocity
    ! solution and the MMS solution
    !
    use meth_params_module, only : URHO, UMY
    use amrex_constants_module
#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: vmmserror(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z, v, rhoinv

#ifdef USE_MASA
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             ! Get error
             v = masa_eval_3d_exact_v(x,y,z)
             rhoinv = 1.d0/dat(i,j,k,URHO)
             vmmserror(i,j,k,1) = dat(i,j,k,UMY) * rhoinv - v

          end do
       end do
    end do
#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

  end subroutine pc_dervmmserror

  subroutine pc_derwmmserror(wmmserror,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_derwmmserror")
    !
    ! This routine will calculate the error between the velocity
    ! solution and the MMS solution
    !
    use meth_params_module, only : URHO, UMZ
    use amrex_constants_module
#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: wmmserror(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z, w, rhoinv

#ifdef USE_MASA
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             ! Get error
             w = masa_eval_3d_exact_w(x,y,z)
             rhoinv = 1.d0/dat(i,j,k,URHO)
             wmmserror(i,j,k,1) = dat(i,j,k,UMZ) * rhoinv - w

          end do
       end do
    end do
#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

  end subroutine pc_derwmmserror

  subroutine pc_derpmmserror(pmmserror,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_derpmmserror")
    !
    ! This routine will calculate the error between the velocity
    ! solution and the MMS solution
    !

    use amrex_constants_module
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use eos_module
#ifdef USE_MASA
    use masa
#endif

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: pmmserror(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z, rhoinv, eint, p
    type(eos_t) :: eos_state

    call build(eos_state)

#ifdef USE_MASA
    do k = lo(3), hi(3)
       z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

       do j = lo(2), hi(2)
          y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

          do i = lo(1), hi(1)
             x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

             ! Get error
             rhoinv = ONE / dat(i,j,k,URHO)
             eos_state % e        = dat(i,j,k,UEINT) * rhoinv
             eos_state % T        = dat(i,j,k,UTEMP)
             eos_state % rho      = dat(i,j,k,URHO)
             eos_state % massfrac = dat(i,j,k,UFS:UFS+nspecies-1) * rhoinv
             eos_state % aux      = dat(i,j,k,UFX:UFX+naux-1) * rhoinv
             call eos_re(eos_state)
             p = masa_eval_3d_exact_p(x,y,z)
             pmmserror(i,j,k,1) = eos_state % p - p

          end do
       end do
    end do
#else
    call bl_error('MASA is not turned on. Turn on with USE_MASA=TRUE.')
#endif

    call destroy(eos_state)

  end subroutine pc_derpmmserror

end module problem_derive_module
