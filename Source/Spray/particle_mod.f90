module particle_mod

  use amrex_fort_module, only: c_real => amrex_real
  use amrex_fort_module, only: bl_spacedim
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t
  
  type, bind(C)  :: particle_t
     real(c_real)    :: pos(bl_spacedim)     !< Position
     real(c_real)    :: vel(bl_spacedim)     !< Particle velocity
     real(c_real)    :: temp                 !< Particle temperature
     real(c_real)    :: diam                 !< Particle diameter
     real(c_real)    :: density              !< Particle density
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type particle_t
  
end module
