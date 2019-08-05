module timestep_module

  implicit none

  public

contains

  ! Courant-condition limited timestep

  subroutine pc_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="pc_estdt")

    use network, only: nspecies, naux
    use eos_module
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX
    use prob_params_module, only: dim
    use amrex_constants_module
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3
    integer          :: i, j, k

    type (eos_t) :: eos_state

    ! Call EOS for the purpose of computing sound speed

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho      = u(i,j,k,URHO )
             eos_state % T        = u(i,j,k,UTEMP)
             eos_state % e        = u(i,j,k,UEINT) * rhoInv
             eos_state % massfrac = u(i,j,k,UFS:UFS+nspecies-1) * rhoInv
             eos_state % aux      = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos_re(eos_state)

             ! Compute velocity and then calculate CFL timestep.

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

             c = eos_state % cs

             dt1 = dx(1)/(c + abs(ux))
             if (dim .ge. 2) then
                dt2 = dx(2)/(c + abs(uy))
             else
                dt2 = dt1
             endif
             if (dim .eq. 3) then
                dt3 = dx(3)/(c + abs(uz))
             else
                dt3 = dt1
             endif

             dt  = min(dt,dt1,dt2,dt3)

          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_estdt


  ! Diffusion-limited timestep

  subroutine pc_estdt_vel_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) &
       bind(C, name="pc_estdt_vel_diffusion")

    use network, only: nspecies, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use prob_params_module, only: dim
    use amrex_constants_module
    use transport_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: dt1, dt2, dt3
    integer          :: i, j, k, n, np
    double precision,allocatable :: rho_inv(:), D(:)

    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    ! dt < 0.5 dx**2 / (dim*D)
    ! where D = mu/rho, and mu is the dynamic viscosity

    np = hi(1)-lo(1)+1
    call build(coeff,np)
    allocate(rho_inv(np))
    allocate(D(np))

    which_trans % wtr_get_mu = .true.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)

          rho_inv(:) = ONE / state(lo(1):hi(1),j,k,URHO)
          do n=1,nspecies
             do i=1,np
                coeff%eos_state(i)%massfrac(n) = state(lo(1)+i-1,j,k,UFS+n-1) * rho_inv(i)
             end do
          end do
          coeff%eos_state(1:np)%T   = state(lo(1):hi(1),j,k,UTEMP)
          coeff%eos_state(1:np)%rho = state(lo(1):hi(1),j,k,URHO) 

          call transport(which_trans, coeff)

          D(:) = coeff%mu(:) * rho_inv(:)

          dt1 = minval(HALF*dx(1)**2/(dim*D(:)))
          dt2 = dt1
          dt3 = dt1

          if (dim .gt. 1) then
             dt2 = minval(HALF*dx(2)**2/(dim*D(:)))
             if (dim .gt. 2) then
                dt3 = minval(HALF*dx(3)**2/(dim*D(:)))
             endif
          endif

          dt  = min(dt,dt1,dt2,dt3)

       enddo
    enddo

    deallocate(D)
    deallocate(rho_inv)
    call destroy(coeff)

  end subroutine pc_estdt_vel_diffusion

  subroutine pc_estdt_temp_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) &
       bind(C, name="pc_estdt_temp_diffusion")

    use network, only: nspecies, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use prob_params_module, only: dim
    use amrex_constants_module
    use transport_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: dt1, dt2, dt3
    integer          :: i, j, k, n, np
    double precision,allocatable :: rho_inv(:), D(:)

    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    ! dt < 0.5 dx**2 / (dim*D)
    ! where D = k/(rho c_v), and k is the conductivity

    np = hi(1)-lo(1)+1
    call build(coeff,np)
    allocate(rho_inv(np))
    allocate(D(np))

    which_trans % wtr_get_lam = .true.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)

          rho_inv(:) = ONE / state(lo(1):hi(1),j,k,URHO)
          do n=1,nspecies
             do i=1,np
                coeff%eos_state(i)%massfrac(n) = state(lo(1)+i-1,j,k,UFS+n-1) * rho_inv(i)
             end do
          end do

          do i=1,np
             coeff%eos_state(i)%T = state(lo(1)+i-1,j,k,UTEMP) ! Note: Assumes this is a good T!

! added by jbb for nonideal . . .should be no op for ideal

             coeff%eos_state(i)%rho = state(lo(1)+i-1,j,k,URHO) ! Note: Assumes this is a good T!
          end do

          call transport(which_trans, coeff)

          do i=1,np
             call eos_cv(coeff % eos_state(i))
          end do

          D(:) = coeff%lam(:) * rho_inv(:) / coeff % eos_state(:) % cv

          dt1 = minval(HALF*dx(1)**2/(dim*D(:)))
          dt2 = dt1
          dt3 = dt1

          if (dim .gt. 1) then
             dt2 = minval(HALF*dx(2)**2/(dim*D(:)))
             if (dim .gt. 2) then
                dt3 = minval(HALF*dx(3)**2/(dim*D(:)))
             endif
          endif

          dt  = min(dt,dt1,dt2,dt3)

       enddo
    enddo

    deallocate(D)
    deallocate(rho_inv)
    call destroy(coeff)

  end subroutine pc_estdt_temp_diffusion

  subroutine pc_estdt_enth_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) &
       bind(C, name="pc_estdt_enth_diffusion")

    use network, only: nspecies, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
         diffuse_cutoff_density
    use prob_params_module, only: dim
    use amrex_constants_module
    use transport_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: dt1, dt2, dt3, rho_inv
    integer          :: i, j, k, np
    double precision :: cond, D

    type (wtr_t) :: which_trans
    type (trv_t) :: coeff

    ! dt < 0.5 dx**2 / (dim*D)
    ! where D = k/(rho c_p), and k is the conductivity

    np = 1
    call build(coeff,np)
    which_trans % wtr_get_lam = .true.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) > diffuse_cutoff_density) then

                rho_inv = ONE/state(i,j,k,URHO)

                ! we need cv, going to assume we trust eint and will compute T
                coeff % eos_state(1) % rho      = state(i,j,k,URHO )
                coeff % eos_state(1) % T        = state(i,j,k,UTEMP)
                coeff % eos_state(1) % e        = state(i,j,k,UEINT) * rho_inv

                coeff % eos_state(1) % massfrac = state(i,j,k,UFS:UFS+nspecies-1) * rho_inv
                if (naux .gt. 0) then
                   coeff % eos_state(1) % aux   = state(i,j,k,UFX:UFX+naux-1)  * rho_inv
                endif

                call eos_re(coeff % eos_state(1))

                ! we also need the conductivity
                if (coeff%eos_state(1)%rho > diffuse_cutoff_density) then
                   call transport(which_trans, coeff)
                   cond = coeff % lam(1)
                else
                   cond = ZERO
                endif

                D = cond*rho_inv/coeff % eos_state(1) % cp

                dt1 = HALF*dx(1)**2/(dim*D)

                if (dim >= 2) then
                   dt2 = HALF*dx(2)**2/(dim*D)
                else
                   dt2 = dt1
                endif

                if (dim == 3) then
                   dt3 = HALF*dx(3)**2/(dim*D)
                else
                   dt3 = dt1
                endif

                dt  = min(dt,dt1,dt2,dt3)

             endif

          enddo
       enddo
    enddo

    call destroy(coeff)

  end subroutine pc_estdt_enth_diffusion


  ! Check whether the last timestep violated any of our stability criteria.
  ! If so, suggest a new timestep which would not.

  subroutine pc_check_timestep(s_old, so_lo, so_hi, &
                               s_new, sn_lo, sn_hi, &
                               lo, hi, &
                               dx, dt_old, dt_new) &
                               bind(C, name="pc_check_timestep")

    use amrex_constants_module, only: HALF, ONE
    use meth_params_module, only: NVAR, URHO, UTEMP, UEINT, UFS, UFX, UMX, UMZ, &
                                  cfl, do_hydro
    use prob_params_module, only: dim
    use network, only: nspecies, naux
    use eos_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: so_lo(3), so_hi(3)
    integer          :: sn_lo(3), sn_hi(3)
    double precision :: s_old(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),NVAR)
    double precision :: s_new(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
    double precision :: dx(3), dt_old, dt_new

    integer          :: i, j, k
    double precision :: v(3), c, rhoinv, tau_CFL
    type (eos_t)     :: eos_state

    call build(eos_state)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoinv = ONE / s_new(i,j,k,URHO)

             ! CFL hydrodynamic stability criterion

             if (do_hydro .eq. 1) then

                eos_state % rho      = s_new(i,j,k,URHO )
                eos_state % T        = s_new(i,j,k,UTEMP)
                eos_state % e        = s_new(i,j,k,UEINT) * rhoinv
                eos_state % massfrac = s_new(i,j,k,UFS:UFS+nspecies-1) * rhoinv
                eos_state % aux      = s_new(i,j,k,UFX:UFX+naux-1) * rhoinv

                call eos_re(eos_state)

                v = HALF * (s_old(i,j,k,UMX:UMZ) * rhoinv + s_new(i,j,k,UMX:UMZ) * rhoinv)

                c = eos_state % cs

                tau_CFL = minval(dx(1:dim) / (c + abs(v(1:dim))))

                dt_new = min(dt_new, cfl * tau_CFL)

             endif

          enddo
       enddo
    enddo

    call destroy(eos_state)

  end subroutine pc_check_timestep

end module timestep_module
