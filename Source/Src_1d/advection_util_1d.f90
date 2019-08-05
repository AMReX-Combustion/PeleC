module advection_util_1d_module

  implicit none

  private

  public normalize_species_fluxes

contains

  subroutine normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

    use network, only : nspecies
    use meth_params_module, only : NVAR, URHO, UFS
    use amrex_constants_module

    implicit none
    
    integer          :: lo(1),hi(1)
    integer          :: flux_l1,flux_h1
    double precision :: flux(flux_l1:flux_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: sum,fac
    
    do i = lo(1),hi(1)+1
       sum = ZERO
       do n = UFS, UFS+nspecies-1
          sum = sum + flux(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = flux(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspecies-1
          flux(i,n) = flux(i,n) * fac
       end do
    end do
    
  end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_h1,lo,hi)

    use network, only : nspecies
    use meth_params_module, only : NVAR, URHO, UFS
    use amrex_constants_module    

    implicit none

    integer          :: lo(1), hi(1)
    integer          :: u_l1,u_h1
    double precision :: u(u_l1:u_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: fac,sum
    
    do i = lo(1),hi(1)
       sum = ZERO
       do n = UFS, UFS+nspecies-1
          sum = sum + u(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = u(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspecies-1
          u(i,n) = u(i,n) * fac
       end do
    end do
    
  end subroutine normalize_new_species
  
end module advection_util_1d_module
