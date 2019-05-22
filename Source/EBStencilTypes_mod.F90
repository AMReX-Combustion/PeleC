module pelec_eb_stencil_types_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none

#if BL_SPACEDIM == 2

  type, bind(c) :: vol_sten
     real(amrex_real) :: val(-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type vol_sten

  type, bind(c) :: face_sten
     real(amrex_real) :: val(-1:1)
     integer :: iv(0:dim-1)
  end type face_sten

  type, bind(c) :: eb_bndry_sten
     real(amrex_real) :: val(-1:1,-1:1)
     real(amrex_real) :: bcval
     integer          :: iv(0:dim-1)
     integer          :: iv_base(0:dim-1)
  end type eb_bndry_sten

#elif BL_SPACEDIM == 3

  type, bind(c) :: vol_sten
     real(amrex_real) :: val(-1:1,-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type vol_sten

  type, bind(c) :: face_sten
     real(amrex_real) :: val(-1:1,-1:1)
     integer :: iv(0:dim-1)
  end type face_sten

  type, bind(c) :: eb_bndry_sten
     real(amrex_real) :: val(-1:1,-1:1, -1:1)
     real(amrex_real) :: bcval
     integer          :: iv(0:dim-1)
     integer          :: iv_base(0:dim-1)
  end type eb_bndry_sten

#endif

  type, bind(c) :: eb_bndry_geom
     real(amrex_real) :: eb_normal(BL_SPACEDIM)
     real(amrex_real) :: eb_centroid(BL_SPACEDIM)
     real(amrex_real) :: eb_area
     real(amrex_real) :: eb_vfrac
!    real(amrex_real) :: eb_apertureX(2)
!    real(amrex_real) :: eb_apertureY(2)
!    real(amrex_real) :: eb_apertureZ(2)
     integer :: iv(0:dim-1)
  end type eb_bndry_geom

end module pelec_eb_stencil_types_module
