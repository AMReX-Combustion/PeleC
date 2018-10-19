module derive_module

  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine pc_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="pc_dervel")
    !
    ! This routine will derive the velocity from the momentum.
    !
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine pc_dervel



  subroutine pc_deruplusc(vel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="pc_deruplusc")

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: rhoInv

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / dat(i,j,k,URHO)

             eos_state % e        = dat(i,j,k,UEINT) * rhoInv
             eos_state % T        = dat(i,j,k,UTEMP)
             eos_state % rho      = dat(i,j,k,URHO)
             eos_state % massfrac = dat(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = dat(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1) + eos_state % cs

          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_deruplusc



  subroutine pc_deruminusc(vel,v_lo,v_hi,nv, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="pc_deruminusc")

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: rhoInv

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / dat(i,j,k,URHO)

             eos_state % e        = dat(i,j,k,UEINT) * rhoInv
             eos_state % T        = dat(i,j,k,UTEMP)
             eos_state % rho      = dat(i,j,k,URHO)
             eos_state % massfrac = dat(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = dat(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1) - eos_state % cs
          end do
       end do
    end do

    call build(eos_state)

  end subroutine pc_deruminusc



  subroutine pc_dermagvel(magvel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="pc_dermagvel")
    !
    ! This routine will derive magnitude of velocity.
    !
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: magvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: dat1inv

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat1inv = 1.d0/dat(i,j,k,1)
             magvel(i,j,k,1) = sqrt( (dat(i,j,k,2) * dat1inv)**2 + &
                  (dat(i,j,k,3) * dat1inv)**2 + &
                  (dat(i,j,k,4) * dat1inv)**2 )
          end do
       end do
    end do

  end subroutine pc_dermagvel



  subroutine pc_derradialvel(radvel,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_derradialvel")
    !
    ! This routine will derive the radial velocity.
    !
    use bl_constants_module
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: radvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: x, y, z, r

    do k = lo(3), hi(3)
       z = xlo(3) + (dble(k-lo(3))+HALF) * delta(3) - center(3)
       do j = lo(2), hi(2)
          y = xlo(2) + (dble(j-lo(2))+HALF) * delta(2) - center(2)
          do i = lo(1), hi(1)
             x = xlo(1) + (dble(i-lo(1))+HALF) * delta(1) - center(1)
             r = sqrt(x*x+y*y+z*z)
             radvel(i,j,k,1) = ( dat(i,j,k,2)*x + &
                  dat(i,j,k,3)*y + &
                  dat(i,j,k,4)*z ) / ( dat(i,j,k,1)*r )
          end do
       end do
    end do

  end subroutine pc_derradialvel



  subroutine pc_dermagmom(magmom,m_lo,m_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="pc_dermagmom")
    !
    ! This routine will derive magnitude of momentum.
    !
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: m_lo(3), m_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: magmom(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nv)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
          end do
       end do
    end do

  end subroutine pc_dermagmom



  subroutine pc_derpres(p,p_lo,p_hi,ncomp_p, &
                        u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="pc_derpres")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: p_lo(3), p_hi(3), ncomp_p
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    double precision :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho      = u(i,j,k,URHO)
             eos_state % T        = u(i,j,k,UTEMP)
             eos_state % e        = u(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             p(i,j,k,1) = eos_state % p
          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_derpres



  subroutine pc_dereint1(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="pc_dereint1")

    use bl_constants_module
    use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN 

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: e_lo(3), e_hi(3), ncomp_e
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    double precision :: rhoInv, ux, uy, uz
    integer          :: i, j, k
    !
    ! Compute internal energy from (rho E).
    !
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoInv = ONE/u(i,j,k,URHO)
             ux = u(i,j,k,UMX)*rhoInv
             uy = u(i,j,k,UMY)*rhoInv
             uz = u(i,j,k,UMZ)*rhoInv
             e(i,j,k,1) = u(i,j,k,UEDEN)*rhoInv-HALF*(ux**2+uy**2+uz**2)
          enddo
       enddo
    enddo

  end subroutine pc_dereint1



  subroutine pc_dereint2(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="pc_dereint2")

    use meth_params_module, only: URHO, UEINT

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: e_lo(3), e_hi(3), ncomp_e
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    !
    ! Compute internal energy from (rho e).
    !
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             e(i,j,k,1) = u(i,j,k,UEINT) / u(i,j,k,URHO)
          enddo
       enddo
    enddo

  end subroutine pc_dereint2



  subroutine pc_dersoundspeed(c,c_lo,c_hi,ncomp_c, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no) &
                              bind(C, name="pc_dersoundspeed")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: c_lo(3), c_hi(3), ncomp_c
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),ncomp_c)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    double precision :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho      = u(i,j,k,URHO)
             eos_state % T        = u(i,j,k,UTEMP)
             eos_state % e        = u(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             c(i,j,k,1) = eos_state % cs
          enddo
       enddo
    enddo

    call build(eos_state)

  end subroutine pc_dersoundspeed



  subroutine pc_dermachnumber(mach,m_lo,m_hi,ncomp_mach, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no) &
                              bind(C, name="pc_dermachnumber")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: URHO, UMX, UMZ, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: m_lo(3), m_hi(3), ncomp_mach
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),ncomp_mach)
    double precision ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u   )
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    double precision :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho      = u(i,j,k,URHO)
             eos_state % T        = u(i,j,k,UTEMP)
             eos_state % e        = u(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             mach(i,j,k,1) = sum(u(i,j,k,UMX:UMZ)**2)**0.5 / u(i,j,k,URHO) / eos_state % cs
          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_dermachnumber



  subroutine pc_derentropy(s,s_lo,s_hi,ncomp_s, &
                           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                           domhi,dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="pc_derentropy")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), ncomp_s
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    double precision :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho      = u(i,j,k,URHO)
             eos_state % T        = u(i,j,k,UTEMP)
             eos_state % e        = u(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             s(i,j,k,1) = eos_state % s
          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_derentropy



  subroutine pc_derenuctimescale(t,t_lo,t_hi,ncomp_t, &
                                 u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                                 domhi,dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C, name="pc_derenuctimescale")

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use network, only: nspec, naux
    use prob_params_module, only: dim
    use eos_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: t_lo(3), t_hi(3), ncomp_t
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    double precision :: t(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3),ncomp_t)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u) ! NVAR, enuc
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    double precision :: rhoInv, eint, enuc, t_s, t_e

    type (eos_t)     :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             enuc = abs(u(i,j,k,ncomp_u))

             if (enuc > 1.d-100) then

                rhoInv = ONE / u(i,j,k,URHO)

                eint = u(i,j,k,UEINT) * rhoInv

                t_e = eint / enuc

                ! Calculate sound-speed

                eos_state % rho      = u(i,j,k,URHO)
                eos_state % T        = u(i,j,k,UTEMP)
                eos_state % e        = eint
                eos_state % massfrac = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                call eos_re(eos_state)

                t_s = minval(dx(1:dim)) / eos_state % cs

                t(i,j,k,1) = t_s / t_e

             else

                t(i,j,k,1) = ZERO

             endif

          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_derenuctimescale



  subroutine pc_derspec(spec,s_lo,s_hi,nv, &
                        dat,d_lo,d_hi,nc,lo,hi,domlo, &
                        domhi,delta,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="pc_derspec")
    !
    ! This routine derives the mass fractions of the species.
    !
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine pc_derspec

  subroutine pc_dermolefrac(spec,s_lo,s_hi,nv, &
                        dat,d_lo,d_hi,nc,lo,hi,domlo, &
                        domhi,delta,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="pc_dermolefrac")
    !
    ! This routine derives the mole fractions of the species.
    !
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use network, only: nspec, naux
    use eos_module
    use bl_constants_module
    
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    double precision :: rhoInv
    integer          :: i, j, k

    type (eos_t)     :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
          
             rhoInv = ONE / dat(i,j,k,URHO)

             eos_state % rho      = dat(i,j,k,URHO)
             eos_state % T        = dat(i,j,k,UTEMP)
             eos_state % e        = dat(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = dat(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux      = dat(i,j,k,UFX:UFX+naux-1) * rhoInv
          
             call eos_ytx(eos_state)

             spec(i,j,k,:) = eos_state % molefrac

          end do
       end do
    end do
 
    call destroy(eos_state)

  end subroutine pc_dermolefrac

  subroutine pc_derlogden(logden,l_lo,l_hi,nd, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="pc_derlogden")
    
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: l_lo(3), l_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3), level, grid_no
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: logden(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),nd)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             logden(i,j,k,1) = dlog10(dat(i,j,k,1))
          end do
       end do
    end do

  end subroutine pc_derlogden



  subroutine pc_dermagvort(vort,v_lo,v_hi,nv, & 
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="pc_dermagvort")
    
    !
    ! This routine will calculate vorticity
    !     

    use bl_constants_module
    use prob_params_module, only: dg

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3), level, grid_no
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: vort(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    double precision :: uy, uz, vx, vz, wx, wy, v1, v2, v3
    double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)

    ldat = ZERO

    uy = ZERO
    uz = ZERO
    vx = ZERO
    vz = ZERO
    wx = ZERO
    wy = ZERO

    !
    ! Convert momentum to velocity.
    !
    do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
             ldat(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,1)
             ldat(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
             ldat(i,j,k,4) = dat(i,j,k,4) / dat(i,j,k,1)
          end do
       end do
    end do
    !
    ! Calculate vorticity.
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             vx = HALF * (ldat(i+1,j,k,3) - ldat(i-1,j,k,3)) / delta(1)
             wx = HALF * (ldat(i+1,j,k,4) - ldat(i-1,j,k,4)) / delta(1)

             if (delta(2) > ZERO) then
                uy = HALF * (ldat(i,j+1,k,2) - ldat(i,j-1,k,2)) / delta(2)
                wy = HALF * (ldat(i,j+1,k,4) - ldat(i,j-1,k,4)) / delta(2)
             endif

             if (delta(3) > ZERO) then
                uz = HALF * (ldat(i,j,k+1,2) - ldat(i,j,k-1,2)) / delta(3)
                vz = HALF * (ldat(i,j,k+1,3) - ldat(i,j,k-1,3)) / delta(3)
             endif

             v1 = wy - vz
             v2 = uz - wx
             v3 = vx - uy
             vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)

          end do
       end do
    end do

  end subroutine pc_dermagvort



  subroutine pc_derdivu(divu,u_lo,u_hi,nd, &
                        dat,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="pc_derdivu")
    !
    ! This routine will calculate the divergence of velocity.
    !

    use bl_constants_module
    use prob_params_module, only: dg

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: divu(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: ulo, uhi, vlo, vhi, wlo, whi

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uhi = dat(i+1*dg(1),j,k,2) / dat(i+1*dg(1),j,k,1)
             ulo = dat(i-1*dg(1),j,k,2) / dat(i-1*dg(1),j,k,1)
             vhi = dat(i,j+1*dg(2),k,3) / dat(i,j+1*dg(2),k,1)
             vlo = dat(i,j-1*dg(2),k,3) / dat(i,j-1*dg(2),k,1)
             whi = dat(i,j,k+1*dg(3),4) / dat(i,j,k+1*dg(3),1)
             wlo = dat(i,j,k-1*dg(3),4) / dat(i,j,k-1*dg(3),1)
             divu(i,j,k,1) = HALF * (uhi-ulo) / delta(1)
             if (delta(2) > ZERO) then
                divu(i,j,k,1) = divu(i,j,k,1) + HALF * (vhi-vlo) / delta(2)
             endif
             if (delta(3) > ZERO) then
                divu(i,j,k,1) = divu(i,j,k,1) + HALF * (whi-wlo) / delta(3)
             endif
          end do
       end do
    end do

  end subroutine pc_derdivu



  subroutine pc_derkineng(kineng,k_lo,k_hi,nk, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="pc_derkineng")
    !
    ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2 + w^2)
    !

    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: k_lo(3), k_hi(3), nk
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             kineng(i,j,k,1) = HALF / dat(i,j,k,1) * ( dat(i,j,k,2)**2 + &
                  dat(i,j,k,3)**2 + &
                  dat(i,j,k,4)**2 )
          end do
       end do
    end do

  end subroutine pc_derkineng



  subroutine pc_derenstrophy(enstrophy,e_lo,e_hi,ne, &
                             dat,d_lo,d_hi,nc, &
                             lo,hi,domlo,domhi,delta, &
                             xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="pc_derenstrophy")
    !
    ! This routine will derive enstrophy  = 1/2 rho (x_vorticity^2 + y_vorticity^2 + z_vorticity^2)
    !

    use bl_constants_module
    use prob_params_module, only: dg
    use meth_params_module, only: URHO, UMX, UMY, UMZ

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: e_lo(3), e_hi(3), ne
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: enstrophy(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ne)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    double precision :: uy, uz, vx, vz, wx, wy, v1, v2, v3
    double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:5)

    ldat = ZERO

    uy = ZERO
    uz = ZERO
    vx = ZERO
    vz = ZERO
    wx = ZERO
    wy = ZERO

    !
    ! Convert momentum to velocity.
    !
    do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
             ldat(i,j,k,2) = dat(i,j,k,UMX) / dat(i,j,k,URHO)
             ldat(i,j,k,3) = dat(i,j,k,UMY) / dat(i,j,k,URHO)
             ldat(i,j,k,4) = dat(i,j,k,UMZ) / dat(i,j,k,URHO)
          end do
       end do
    end do
    !
    ! Calculate vorticity.
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             vx = HALF * (ldat(i+1,j,k,3) - ldat(i-1,j,k,3)) / delta(1)
             wx = HALF * (ldat(i+1,j,k,4) - ldat(i-1,j,k,4)) / delta(1)

             if (delta(2) > ZERO) then
                uy = HALF * (ldat(i,j+1,k,2) - ldat(i,j-1,k,2)) / delta(2)
                wy = HALF * (ldat(i,j+1,k,4) - ldat(i,j-1,k,4)) / delta(2)
             endif

             if (delta(3) > ZERO) then
                uz = HALF * (ldat(i,j,k+1,2) - ldat(i,j,k-1,2)) / delta(3)
                vz = HALF * (ldat(i,j,k+1,3) - ldat(i,j,k-1,3)) / delta(3)
             endif

             v1 = wy - vz
             v2 = uz - wx
             v3 = vx - uy
             ldat(i,j,k,5) = sqrt(v1*v1 + v2*v2 + v3*v3)

          end do
       end do
    end do
    !
    ! Calculate enstrophy
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             enstrophy(i,j,k,1) = HALF * dat(i,j,k,URHO) * ldat(i,j,k,5)**2
          end do
       end do
    end do

  end subroutine pc_derenstrophy



  subroutine pc_dernull(kineng,k_lo,k_hi,nk, &
                        dat,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="pc_dernull")
    !
    ! This routine is used by particle_count.  Yes it does nothing.
    !
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: k_lo(3), k_hi(3), nk
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    double precision :: delta(3), xlo(3), time, dt
    double precision :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
    double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

  end subroutine pc_dernull

end module derive_module
