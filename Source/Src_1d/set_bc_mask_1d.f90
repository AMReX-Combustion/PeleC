module bc_mask_module

  private
  public set_bc_mask

contains

  subroutine set_bc_mask(lo, hi, domlo, domhi, &
                         x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2) &
                         bind(C, name="set_bc_mask")
      
  
    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi
   
    implicit none
     
  
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: domlo(2), domhi(2)
  
    integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2
  
    integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2)

  
    x_bcMask(domlo(1),  x_bcMask_l2:x_bcMask_h2) = physbc_lo(1)
    x_bcMask(domhi(1)+1,x_bcMask_l2:x_bcMask_h2) = physbc_hi(1)
    

   
  end subroutine set_bc_mask

end module bc_mask_module


