module bc_mask_module

  private
  public set_bc_mask

contains

  subroutine set_bc_mask(lo, hi, domlo, domhi, &
                         x_bcMask, x_bcMask_l1, x_bcMask_l2, x_bcMask_l3, x_bcMask_h1, x_bcMask_h2, x_bcMask_h3, &
                         y_bcMask, y_bcMask_l1, y_bcMask_l2, y_bcMask_l3, y_bcMask_h1, y_bcMask_h2, y_bcMask_h3, &
                         z_bcMask, z_bcMask_l1, z_bcMask_l2, z_bcMask_l3, z_bcMask_h1, z_bcMask_h2, z_bcMask_h3) &
                         bind(C, name="set_bc_mask")
      
  
    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi
   
    implicit none
     
  
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
  
    integer, intent(in) :: x_bcMask_l1, x_bcMask_l2, x_bcMask_h1, x_bcMask_h2, x_bcMask_l3, x_bcMask_h3
    integer, intent(in) :: y_bcMask_l1, y_bcMask_l2, y_bcMask_h1, y_bcMask_h2, y_bcMask_l3, y_bcMask_h3
    integer, intent(in) :: z_bcMask_l1, z_bcMask_l2, z_bcMask_h1, z_bcMask_h2, z_bcMask_l3, z_bcMask_h3
  
    integer, intent(inout) :: x_bcMask(x_bcMask_l1:x_bcMask_h1,x_bcMask_l2:x_bcMask_h2,x_bcMask_l3:x_bcMask_h3)
    integer, intent(inout) :: y_bcMask(y_bcMask_l1:y_bcMask_h1,y_bcMask_l2:y_bcMask_h2,y_bcMask_l3:y_bcMask_h3)
    integer, intent(inout) :: z_bcMask(z_bcMask_l1:z_bcMask_h1,z_bcMask_l2:z_bcMask_h2,z_bcMask_l3:z_bcMask_h3)
  
    if (x_bcMask_l1 == domlo(1)) then
    ! Left x face
      x_bcMask(domlo(1),  x_bcMask_l2:x_bcMask_h2, x_bcMask_l3:x_bcMask_h3) = physbc_lo(1)
    endif
    
    if (x_bcMask_h1 == domhi(1)+1) then
    ! Right x face
      x_bcMask(domhi(1)+1,x_bcMask_l2:x_bcMask_h2, x_bcMask_l3:x_bcMask_h3) = physbc_hi(1)
    end if
    
    if (y_bcMask_l2 == domlo(2)) then
    ! Left y face
      y_bcMask(y_bcMask_l1:y_bcMask_h1,domlo(2),y_bcMask_l3:y_bcMask_h3) = physbc_lo(2)
    end if
    
    if (y_bcMask_h2 == domhi(2)+1) then
      y_bcMask(y_bcMask_l1:y_bcMask_h1,domhi(2)+1,y_bcMask_l3:y_bcMask_h3) = physbc_hi(2)
    end if
    
    if (z_bcMask_l3 == domlo(3)) then
      z_bcMask(z_bcMask_l1:z_bcMask_h1, z_bcMask_l2:z_bcMask_h2, domlo(3)) = physbc_lo(3)
    end if
    
    if (z_bcMask_h3 == domhi(3)+1) then
      z_bcMask(z_bcMask_l1:z_bcMask_h1, z_bcMask_l2:z_bcMask_h2, domhi(3)+1) = physbc_hi(3)
    end if

  end subroutine set_bc_mask

end module bc_mask_module

