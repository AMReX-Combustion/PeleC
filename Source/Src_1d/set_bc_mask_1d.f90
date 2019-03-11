module bc_mask_module

  private
  public set_bc_mask

contains

  subroutine set_bc_mask(lo, hi, domlo, domhi, &
                         x_bcMask, x_bcMask_l1, x_bcMask_h1) &
                         bind(C, name="set_bc_mask")
      
  
    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi
   
    implicit none
     
  
    integer, intent(in) :: lo(1), hi(1)
    integer, intent(in) :: domlo(1), domhi(1)
  
    integer, intent(in) :: x_bcMask_l1, x_bcMask_h1
  
    integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1)

  
    if (x_bcMask_l1 == domlo(1)) then
    ! Left x face
      x_bcMask(domlo(1)) = physbc_lo(1)
    end if
    
    if (x_bcMask_h1 == domhi(1)+1) then
    ! Right x face
      x_bcMask(domhi(1)+1) = physbc_hi(1)
    end if

   
  end subroutine set_bc_mask

end module bc_mask_module


