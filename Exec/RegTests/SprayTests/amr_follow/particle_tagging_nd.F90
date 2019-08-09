module particle_tagging_module

  implicit none

  public

contains

  ! All tagging subroutines in this file must be threadsafe because
  ! they are called inside OpenMP parallel regions.

  ! ::: -----------------------------------------------------------
  ! ::: INPUTS/OUTPUTS:
  ! ::: 
  ! ::: tag      <=  integer tag array
  ! ::: lo,hi     => index extent of work region
  ! ::: set       => integer value to tag cell for refinement
  ! ::: clear     => integer value to untag cell
  ! ::: temp      => temperature array
  ! ::: np        => number of components in temp array (should be 1)
  ! ::: domlo,hi  => index extent of problem domain
  ! ::: delta     => cell spacing
  ! ::: xlo       => physical location of lower left hand
  ! :::              corner of work region
  ! ::: problo    => phys loc of lower left corner of prob domain
  ! ::: time      => problem evolution time
  ! ::: level     => refinement level of this array
  ! ::: -----------------------------------------------------------

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag cells where particles are
  ! ::: -----------------------------------------------------------

  subroutine pc_particle_tag(tag,taglo,taghi, &
                          set,clear, &
                          var,varlo,varhi, &
                          lo,hi,nd,domlo,domhi, &
                          dx,xlo,problo,time,level) &
                          bind(C, name="pc_particle_tag")

    implicit none
    integer          :: set, clear, nd, level
    integer          :: taglo(3), taghi(3)
    integer          :: varlo(3), varhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: var(varlo(1):varhi(1),varlo(2):varhi(2),varlo(3):varhi(3))
    double precision :: dx(3), xlo(3), problo(3), time
    integer          :: i,j,k
   
    do k = varlo(3),varhi(3)
    do j = varlo(2),varhi(2)
    do i = varlo(1),varhi(1)
       if (var(i,j,k) .ge. 1) then
          tag(i,j,k) = set
       end if
    end do
    end do
    end do

  end subroutine pc_particle_tag

end module particle_tagging_module
