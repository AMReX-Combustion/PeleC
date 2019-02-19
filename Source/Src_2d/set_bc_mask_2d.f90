module bc_mask_module

  private
  public set_bc_mask

contains

  subroutine set_bc_mask(lo, hi, domlo, domhi, &
                         x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2, &
                         y_bcMask, y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2) &
                         bind(C, name="set_bc_mask")
      
  
    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi
   
    implicit none
     
  
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: domlo(2), domhi(2)
  
    integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2
    integer, intent(in) :: y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2
  
    integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2)
    integer, intent(inout) :: y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2)
  
    if (x_bcMask_l1 == domlo(1)) then
      x_bcMask(domlo(1),  x_bcMask_l2:x_bcMask_h2) = physbc_lo(1)
    end if
    
    if (x_bcMask_h1 == domhi(1)+1) then
      x_bcMask(domhi(1)+1,x_bcMask_l2:x_bcMask_h2) = physbc_hi(1)
    end if
    
    if (y_bcMask_l2 == domlo(2)) then
      y_bcMask(y_bcMask_l1:y_bcMask_h1,domlo(2)) = physbc_lo(2)
    end if

    if (y_bcMask_h2 == domhi(2)+1) then
      y_bcMask(y_bcMask_l1:y_bcMask_h1,domhi(2)+1) = physbc_hi(2)
    end if
   
  end subroutine set_bc_mask

end module bc_mask_module


