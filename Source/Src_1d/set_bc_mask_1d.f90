   subroutine set_bc_mask(lo, hi, domlo, domhi, &
                          bcMask, bcMask_l1, bcMask_h1) &
      bind(C, name="set_bc_mask")
 
      use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, &
                                   Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall
    
      implicit none
   
     integer, intent(in   ) :: lo(1), hi(1), domlo(1), domhi(1)
     integer, intent(in   ) :: bcMask_l1, bcMask_h1
     integer, intent(inout) :: bcMask(bcMask_l1:bcMask_h1,2)

     ! comp=1 holds physical bc type (Inflow, Outflow, NoSlipWall
     ! comp=2 hold flag for whether is adiabatic (=0) or isothermal (=1)   TODO: check this

     if (bcMask_l1 < domlo(1)) then
     ! Left x face
        bcMask(domlo(1) ,1) = physbc_lo(1)
     end if

     if (bcMask_h1 > domhi(1)) then
     ! Right x face
        bcMask(domhi(1)+1, 1) = physbc_hi(1)
     end if

   end subroutine set_bc_mask
