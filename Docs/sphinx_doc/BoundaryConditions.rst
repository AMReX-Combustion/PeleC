
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _BCs:

Boundary Conditions
-------------------

PeleC manages boundary conditions in a form consistent with many AMReX codes. This means that the user is responsible for providing functions that will appropriately fill ghost zones with appropriate state data. This is generally done in the case specific directory by defining the `pc_hypfill` and `pc_denfill` functions. Internally, these functions are registered as the boundary filling functions for all of the components of the `State_Type` type. 

The function prototype is:

.. highlight:: fortran

::

 	subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) bind(C, name="pc_hypfill")


In the inputs files, the user may specify boundary condition flags that are labeled 'Interior, Inflow, Outflow, SlipWall, NoSlipWall', and it is left to the user to define what the correct boundary condition treatment is for the problem at hand. In the specific example in the PMF directory, the routine `bcnormal` is used to fill ghost zones with appropriate state data based on the boundary condition flag for that face. Specification of boundary conditions at EB faces is a special case that will be discussed in the EB section of this guide.  The general philosophy in PeleC is that the user is to provide routines that load *ghost cells outside the domain* with boundary conditions *at the domain edge*. Within PeleC, interpolation/extrapolation or other appropriate treatment is used to apply the boundary condition at the appropriate location. 