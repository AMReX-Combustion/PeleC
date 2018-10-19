   subroutine set_bc_mask(lo, hi, domlo, domhi, &
                          bcMask, bcMask_l1, bcMask_l2, bcMask_l3, bcMask_h1, bcMask_h2, bcMask_h3) &
      bind(C, name="set_bc_mask")
 
      use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi, &
                                   Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall
    
      implicit none
   
     integer, intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
     integer, intent(in   ) :: bcMask_l1, bcMask_l2, bcMask_l3, bcMask_h1, bcMask_h2, bcMask_h3
     integer, intent(inout) :: bcMask(bcMask_l1:bcMask_h1,bcMask_l2:bcMask_h2,bcMask_l3:bcMask_h3, 2)

     ! comp=1 holds physical bc type (Inflow, Outflow, NoSlipWall
     ! comp=2 hold flag for whether is adiabatic (=0) or isothermal (=1)   TODO: check this

     if (bcMask_l1 < domlo(1)) then
     ! Left x face
        bcMask(domlo(1), bcMask_l2:bcMask_h2, bcMask_l3:bcMask_h3  ,1) = physbc_lo(1)
     end if

     if (bcMask_h1 > domhi(1)) then
     ! Right x face
        bcMask(domhi(1)+1, bcMask_l2:bcMask_h2, bcMask_l3:bcMask_h3, 1) = physbc_hi(1)
     end if

     if (bcMask_l2 < domlo(2)) then
      ! Left y face
        bcMask(bcMask_l1:bcMask_h1,domlo(2), bcMask_l3:bcMask_h3 ,1) = physbc_lo(2)
     end if

     if (bcMask_h2 > domhi(2)) then
        bcMask(bcMask_l1:bcMask_h1,domhi(2)+1, bcMask_l3:bcMask_h3,1) = physbc_hi(2)
     end if

    if (bcMask_l3 < domlo(3)) then
        bcMask(bcMask_l1:bcMask_h1, bcMask_l2:bcMask_h2, domlo(3)  ,1) = physbc_lo(3)
     end if

     if (bcMask_h3 > domhi(3)) then
        bcMask(bcMask_l1:bcMask_h1, bcMask_l2:bcMask_h2, domhi(3)+1,1) = physbc_hi(3)
     end if
   end subroutine set_bc_mask
