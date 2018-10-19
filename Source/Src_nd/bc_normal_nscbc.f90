  subroutine bcnormal_nscbc(x,u_int,u_ext,dir,sgn,bc_type,bc_params,rho_only)

    use meth_params_module, only: nb_nscbc_params
  
    implicit none

    double precision :: x(3)
    double precision :: u_int(*),u_ext(*)
    logical rho_only
    integer :: dir,sgn
    integer, intent(out) :: bc_type
    double precision, intent(out) :: bc_params(nb_nscbc_params)


  end subroutine bcnormal_nscbc