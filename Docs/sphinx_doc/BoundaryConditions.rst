
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _BCs:

Boundary Conditions
-------------------

PeleC manages boundary conditions in a form consistent with many AMReX codes. Ghost cell data is updated at each AMR level and 
fluxes are computed the same way as it is done for interior cells.
There are some pre-specified boundary conditions where user intervention is not required. These are:

* *Interior* - same as periodic boundary conditions. One has to make sure the is_periodic flag is also set correctly in the inputs file
* *Outflow*  - First order extrapolation for all conserved quantities is used to populate the ghost cells
* *SlipWall* - All conserved quantities and the tangential momentum component is reflected from interior cells without 
  sign change (REFLECT_EVEN) while the normal component is reflected with a sign change (REFLECT_ODD)
* *NoSlipWall* - REFLECT_EVEN is applied to all conserved quantities except for both tangential and normal momentum components which are updated 
  using REFLECT_ODD  

The boundary condition, *Inflow*, requires the user to provide functions that will 
appropriately fill ghost zones with appropriate state data. 
This is generally done in the *case* specific directory by defining the `pc_hypfill` and `pc_denfill` functions in the bc_fill_nd.f90 file.
The former is used to update all conserved variables while the latter is used for single component (only density) fills during fillpatch operations.
Internally, these functions are registered as the boundary filling functions for all of the components of the `State_Type` type. 
The bc_fill_nd.f90 file in most cases can be reused for a custom problem with modifications to `bcnormal` function.
This function operates on a per cell basis. The user updates the set of conserved variables at the ghost cell, `u_ext` 
given the physical coordinates `x`, values of conserved variables in the interior, `u_int`, coordinate direction `dir` and 
outward normal direction `sgn`.

.. Don't know how useful it is show the function prototype here

.. The function prototype is:

.. .. highlight:: fortran

.. 	subroutine pc_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) bind(C, name="pc_hypfill")

Below is an example to achieve a moving wall using an *Inflow* boundary condition specification.

::

    subroutine bcnormal(x,u_int,u_ext,dir,sgn,rho_only)

        use probdata_module, only : vel_wall
        use network, only: nspec
        use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP

        implicit none

        double precision :: x(3)
        double precision :: u_int(*),u_ext(*)
        logical rho_only
        integer :: dir,sgn

        double precision :: vcell    !internal velocity
        double precision :: vghost   !velocity in ghost cell

        if (rho_only .EQV. .TRUE. ) then
            u_ext(URHO) = u_int(URHO)
        else

        !moving wall
        vcell  = u_int(UMX)/u_int(URHO)
        vghost = 2.d0*vel_wall-vcell  !reflection

        u_ext(URHO)   =  u_int(URHO)
        u_ext(UMX)    =  u_int(URHO)*vghost
        u_ext(UMY)    = -u_int(UMY)
        u_ext(UMZ)    = -u_int(UMZ)
        u_ext(UEINT)  =  u_int(UEINT)
        u_ext(UEDEN)  =  u_int(UEINT)+0.5d0*(u_ext(UMX)**2+u_ext(UMY)**2+u_ext(UMZ)**2)/u_ext(URHO)
        u_ext(UFS:UFS+nspec-1)  = u_int(UFS:UFS+nspec-1)

        endif

    end subroutine bcnormal


There are specific examples in many of the test cases that the user may adapt to their requirements.
In the *VIF* case, `bcnormal` function updates the ghost cells with a turbulent subsonic velocity inlet condition with interpolation of pressure from the interior. 
In the *Oblqshock* case, a supersonic inflow boundary is used with velocity specification at a user defined angle 


Specification of boundary conditions at EB faces is a special case that will be discussed in the EB section of this guide.  
The general philosophy in PeleC is that the user is to provide routines that load *ghost cells outside the domain* with boundary conditions *at the domain edge*. Within PeleC, interpolation/extrapolation or other appropriate treatment is used to apply the boundary condition at the appropriate location. 

