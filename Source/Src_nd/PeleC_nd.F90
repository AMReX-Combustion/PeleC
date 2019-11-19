subroutine pc_network_init() bind(C, name="pc_network_init")

  use network, only: network_init
  implicit none

  call network_init()

end subroutine pc_network_init


! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_network_close() bind(C, name="pc_network_close")

  use network, only: network_close
  implicit none

  call network_close()

end subroutine pc_network_close


! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_transport_init() bind(C, name="pc_transport_init")

  use transport_module, only: transport_init
  implicit none

  call transport_init()

end subroutine pc_transport_init


! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_transport_close() bind(C, name="pc_transport_close")

  use transport_module, only: transport_close
  implicit none

  call transport_close()

end subroutine pc_transport_close


! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_extern_init(name,namlen) bind(C, name="pc_extern_init")

  ! initialize the external runtime parameters in
  ! extern_probin_module

  implicit none
  integer :: namlen
  integer :: name(namlen)

  call runtime_init(name,namlen)

end subroutine pc_extern_init

! :::
! ::: ----------------------------------------------------------------
! :::
#ifdef REACTIONS
subroutine pc_reactor_init() bind(C, name="pc_reactor_init")

#ifdef USE_SUNDIALS_PP
  use cvode_module, only : reactor_init 
#else
  use reactor_module, only: reactor_init
#endif

  implicit none

#ifdef _OPENMP
!$omp parallel
#endif
#ifdef USE_SUNDIALS_PP
  call reactor_init(1,1)
#else
  call reactor_init(1)
#endif
#ifdef _OPENMP
!$omp end parallel
#endif

end subroutine pc_reactor_init

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_reactor_close() bind(C, name="pc_reactor_close")

#ifdef USE_SUNDIALS_PP
  use cvode_module, only: reactor_close
#else
  use reactor_module, only: reactor_close
#endif
  implicit none

#ifdef _OPENMP
!$omp parallel
#endif
  call reactor_close()
#ifdef _OPENMP
!$omp end parallel
#endif

end subroutine pc_reactor_close
#endif

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_num_spec(nspecies_out) bind(C, name="get_num_spec")

  use network, only : nspecies

  implicit none

  integer, intent(out) :: nspecies_out

  nspecies_out = nspecies

end subroutine get_num_spec

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_num_aux(naux_out) bind(C, name="get_num_aux")

  use network, only : naux

  implicit none

  integer, intent(out) :: naux_out

  naux_out = naux

end subroutine get_num_aux

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_spec_names(spec_names_out,ispec,len) &
     bind(C, name="get_spec_names")

  use network, only : spec_names

  implicit none

  integer, intent(in   ) :: ispec
  integer, intent(inout) :: len
  integer, intent(inout) :: spec_names_out(len)

  ! Local variables
  integer   :: i

  len = len_trim(spec_names(ispec+1))

  do i = 1,len
     spec_names_out(i) = ichar(spec_names(ispec+1)(i:i))
  end do

end subroutine get_spec_names

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_aux_names(aux_names_out,iaux,len) &
     bind(C, name="get_aux_names")

  use network, only : aux_names

  implicit none

  integer, intent(in   ) :: iaux
  integer, intent(inout) :: len
  integer, intent(inout) :: aux_names_out(len)

  ! Local variables
  integer   :: i

  len = len_trim(aux_names(iaux+1))

  do i = 1,len
     aux_names_out(i) = ichar(aux_names(iaux+1)(i:i))
  end do

end subroutine get_aux_names

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine set_amr_info(level_in, iteration_in, ncycle_in, time_in, dt_in) &
     bind(C, name="set_amr_info")

  use amrinfo_module, only: amr_level, amr_iteration, amr_ncycle, amr_time, amr_dt

  implicit none

  integer, intent(in) :: level_in, iteration_in, ncycle_in
  double precision, intent(in) :: time_in, dt_in

  if (level_in .ge. 0) then
     amr_level = level_in
  endif

  if (iteration_in .ge. 0) then
     amr_iteration = iteration_in
  endif

  if (ncycle_in .ge. 0) then
     amr_ncycle = ncycle_in
  endif

  if (time_in .ge. 0.0) then
     amr_time = time_in
  endif

  if (dt_in .ge. 0.0) then
     amr_dt = dt_in
  endif

end subroutine set_amr_info

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_method_params(nGrowHyp,cQTHERM,cQVAR,cQRHO,cQU,cQV,cQW,cQGAME,cQPRES,&
     cQREINT,cQTEMP,cQFA,cQFS,cQFX,cNQAUX,cQGAMC,cQC,cQCSML,cQDPDR,cQDPDE,cQRSPEC) &
     bind(C, name="get_method_params")

  ! Passing data from f90 back to C++

  use meth_params_module

  implicit none

  integer, intent(out) :: ngrowHyp, cQTHERM, cQVAR, cQRHO, cQU, cQV, cQW, cQGAME, cQPRES,&
                          cQREINT, cQTEMP, cQFA, cQFS, cQFX, cNQAUX, cQGAMC, cQC, cQCSML,&
                          cQDPDR, cQDPDE, cQRSPEC

  nGrowHyp = NHYP
  cQTHERM = QTHERM
  cQVAR = QVAR
  cQRHO = QRHO - 1
  cQU = QU - 1
  cQV = QV - 1
  cQW = QW - 1
  cQGAME = QGAME - 1
  cQPRES = QPRES - 1
  cQREINT = QREINT - 1
  cQTEMP = QTEMP - 1
  cQFA = QFA - 1
  cQFS = QFS - 1
  cQFX = QFX - 1
  cNQAUX = NQAUX
  cQGAMC = QGAMC - 1
  cQC = QC - 1
  cQCSML = QCSML - 1
  cQDPDR = QDPDR - 1
  cQDPDE = QDPDE - 1
  cQRSPEC = QRSPEC - 1

end subroutine get_method_params

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine allocate_outflow_data(np,nc) &
     bind(C, name="allocate_outflow_data")

  use meth_params_module, only: outflow_data_old, outflow_data_new, outflow_data_allocated

  implicit none

  integer, intent(in) :: np,nc

  if (.not. outflow_data_allocated) then
     allocate(outflow_data_old(nc,np))
     allocate(outflow_data_new(nc,np))
  end if

  outflow_data_allocated = .true.

end subroutine allocate_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
subroutine set_old_outflow_data(radial,time,np,nc) &
     bind(C, name="set_old_outflow_data")

  ! Passing data from C++ to f90

  use meth_params_module, only: outflow_data_old, outflow_data_old_time

  implicit none

  double precision, intent(in) :: radial(nc,np)
  double precision, intent(in) :: time
  integer         , intent(in) :: np,nc

  ! Do this so the routine has the right size
  deallocate(outflow_data_old)
  allocate(outflow_data_old(nc,np))

  outflow_data_old(1:nc,1:np) = radial(1:nc,1:np)

  outflow_data_old_time = time

end subroutine set_old_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
subroutine set_new_outflow_data(radial,time,np,nc) &
     bind(C, name="set_new_outflow_data")

  ! Passing data from C++ to f90

  use meth_params_module, only: outflow_data_new, outflow_data_new_time

  implicit none

  double precision, intent(in) :: radial(nc,np)
  double precision, intent(in) :: time
  integer         , intent(in) :: np,nc

  ! Do this so the routine has the right size
  deallocate(outflow_data_new)
  allocate(outflow_data_new(nc,np))

  outflow_data_new(1:nc,1:np) = radial(1:nc,1:np)

  outflow_data_new_time = time

end subroutine set_new_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
subroutine swap_outflow_data() bind(C, name="swap_outflow_data")

  use meth_params_module, only: outflow_data_new, outflow_data_new_time, &
       outflow_data_old, outflow_data_old_time

  implicit none

  integer                       :: np,nc

  nc = size(outflow_data_new,dim=1)
  np = size(outflow_data_new,dim=2)

  if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
     ! Do this so the routine has the right size
     deallocate(outflow_data_old)
     allocate(outflow_data_old(nc,np))
  end if

  if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
     print *,'size of old and new dont match in swap_outflow_data '
     call bl_error("Error:: PeleC_nd.f90 :: swap_outflow_data")
  end if

  outflow_data_old(1:nc,1:np) = outflow_data_new(1:nc,1:np)

  if (outflow_data_new_time .ge. 0.d0) &
       outflow_data_old_time = outflow_data_new_time
  outflow_data_new_time = -1.d0

end subroutine swap_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine set_method_params(dm,Density,Xmom,Eden,Eint,Temp, &
     FirstAdv,FirstSpec,FirstAux,numadv, &
     diffuse_cutoff_density_in, &
     pstate_loc, pstate_vel, pstate_T, pstate_dia, pstate_rho, pstate_ys, &
     pfld_vel, pfld_rho, pfld_T, pfld_p, pfld_ys) &
     bind(C, name="set_method_params")

  use meth_params_module
  use network, only : nspecies, naux
  use eos_module, only : eos_init, eos_get_small_dens, eos_get_small_temp
  use transport_module, only : transport_init
  use amrex_constants_module, only : ZERO, ONE

  implicit none

  integer, intent(in) :: dm
  integer, intent(in) :: Density, Xmom, Eden, Eint, Temp, &
       FirstAdv, FirstSpec, FirstAux
  integer, intent(in) :: pstate_loc, pstate_vel, pstate_T, pstate_dia, &
                         pstate_rho, pstate_ys, pfld_vel, pfld_rho, pfld_T, &
                         pfld_p, pfld_ys
  integer, intent(in) :: numadv
  double precision, intent(in) :: diffuse_cutoff_density_in
  integer :: iadv, ispec

  integer :: QLAST

  integer :: ioproc

  !---------------------------------------------------------------------
  ! conserved state components
  !---------------------------------------------------------------------

  ! NTHERM: number of thermodynamic variables (rho, 3 momenta, rho*e, rho*E, T)
  ! NVAR  : number of total variables in initial system
  NTHERM = 7

  NVAR = NTHERM + nspecies + naux + numadv

  nadv = numadv

  ! We use these to index into the state "U"
  URHO  = Density   + 1
  UMX   = Xmom      + 1
  UMY   = Xmom      + 2
  UMZ   = Xmom      + 3
  UEDEN = Eden      + 1
  UEINT = Eint      + 1
  UTEMP = Temp      + 1

  if (numadv .ge. 1) then
     UFA   = FirstAdv  + 1
  else
     UFA = 1
  end if

  UFS   = FirstSpec + 1

  if (naux .ge. 1) then
     UFX = FirstAux  + 1
  else
     UFX = 1
  end if

  USHK  = -1
  !---------------------------------------------------------------------
  ! primitive state components
  !---------------------------------------------------------------------

  ! QTHERM: number of primitive variables: rho, p, (rho e), T + 3 velocity components 
  ! QVAR  : number of total variables in primitive form

  QTHERM = NTHERM + 1 ! the + 1 is for QGAME which is always defined in primitive mode

  QVAR = QTHERM + nspecies + naux + numadv
  
  ! NQ will be the number of hydro + radiation variables in the primitive
  ! state.  Initialize it just for hydro here
  NQ = QVAR

  ! We use these to index into the state "Q"
  !QRHO  = 1

  !QU    = 2
  !QV    = 3
  !QW    = 4

  !QGAME = 5

  QLAST   = QGAME

  !QPRES   = QLAST + 1
  !QREINT  = QLAST + 2

  !QTEMP   = QTHERM ! = QLAST + 3

  if (numadv >= 1) then
     QFA = QTHERM + 1
     !QFS = QFA + numadv

  else
     QFA = 1   ! density
     !QFS = QTHERM + 1

  end if

  if (naux >= 1) then
     QFX = QFS + nspecies

  else
     QFX = 1

  end if

  ! The NQAUX here are auxiliary quantities (game, gamc, c, csml, dpdr, dpde, Rspecific)
  ! that we create in the primitive variable call but that do not need to
  ! participate in tracing.

  NQAUX = 6
  QGAMC   = 1
  QC      = 2
  QCSML   = 3
  QDPDR   = 4
  QDPDE   = 5
  QRSPEC   = 6
  ! easy indexing for the passively advected quantities.  This
  ! lets us loop over all groups (advected, species, aux)
  ! in a single loop.
  allocate(qpass_map(QVAR))
  allocate(upass_map(NVAR))

  ! Transverse velocities

  if (dm == 1) then
     upass_map(1) = UMY
     qpass_map(1) = QV

     upass_map(2) = UMZ
     qpass_map(2) = QW

     npassive = 2

  else if (dm == 2) then
     upass_map(1) = UMZ
     qpass_map(1) = QW

     npassive = 1
  else
     npassive = 0
  endif

  do iadv = 1, nadv
     upass_map(npassive + iadv) = UFA + iadv - 1
     qpass_map(npassive + iadv) = QFA + iadv - 1
  enddo
  npassive = npassive + nadv

  if (QFS > -1) then
     do ispec = 1, nspecies+naux
        upass_map(npassive + ispec) = UFS + ispec - 1
        qpass_map(npassive + ispec) = QFS + ispec - 1
     enddo
     npassive = npassive + nspecies + naux
  endif

  !---------------------------------------------------------------------
  ! Particle state indices
  !---------------------------------------------------------------------
  PLOC  = 1 + pstate_loc
  PVEL  = 1 + pstate_vel
  PTEMP = 1 + pstate_T
  PDIA  = 1 + pstate_dia
  PRHO  = 1 + pstate_rho
  PSPC  = 1 + pstate_ys
  
  PFVEL  = 1 + pfld_vel
  PFRHO  = 1 + pfld_rho
  PFTEMP = 1 + pfld_T
  PFP    = 1 + pfld_p
  PFSPC  = 1 + pfld_ys

  !---------------------------------------------------------------------
  ! other initializations
  !---------------------------------------------------------------------

  ! This is a routine which links to the C++ ParallelDescriptor class

  call bl_pd_is_ioproc(ioproc)

  diffuse_cutoff_density       = diffuse_cutoff_density_in

  ! Note that the EOS may modify choice in namelist because of its
  ! internal limitations, so the small_dens and small_temp
  ! may be modified coming back out of this routine.

  call eos_init(small_dens=small_dens, small_temp=small_temp)

  ! Update device variables

  !$acc update &
  !$acc device(NTHERM, NVAR) &
  !$acc device(NQ) &
  !$acc device(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX) &
  !$acc device(USHK) &
  !$acc device(QTHERM, QVAR) &
  !$acc device(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME) &
  !$acc device(QFA, QFS, QFX) &
  !$acc device(NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE, QRSPEC) &
  !$acc device(small_dens, small_temp)

end subroutine set_method_params


subroutine clear_method_params() &
     bind(C, name="clear_method_params")

  use meth_params_module
  implicit none

  ! call to match parallel_initialize()?
  deallocate(qpass_map)
  deallocate(upass_map)
  ! call to match eos_init()?

end subroutine clear_method_params


subroutine init_godunov_indices() bind(C, name="init_godunov_indices")

  use meth_params_module, only : GDRHO, GDU, GDV, GDW, GDPRES, GDGAME, NGDNV, &
       QU, QV, QW

  implicit none

  NGDNV = 6
  GDRHO = 1
  GDU = 2
  GDV = 3
  GDW = 4
  GDPRES = 5
  GDGAME = 6

  ! sanity check
  if ((QU /= GDU) .or. (QV /= GDV) .or. (QW /= GDW)) then
     call bl_error("ERROR: velocity components for godunov and primitive state are not aligned")
  endif

end subroutine init_godunov_indices

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in,&
     Interior_in, UserBC_in, Inflow_in, Outflow_in, &
     Symmetry_in, SlipWall_in, NoSlipWall_in, &
     coord_type_in, &
     problo_in, probhi_in, center_in) &
     bind(C, name="set_problem_params")

  ! Passing data from C++ into f90

  use amrex_constants_module, only: ZERO
  use prob_params_module
  implicit none

  integer, intent(in) :: dm
  integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
  integer, intent(in) :: Interior_in, UserBC_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
  integer, intent(in) :: coord_type_in
  double precision, intent(in) :: problo_in(dm), probhi_in(dm), center_in(dm)

  dim = dm

  physbc_lo(1:dm) = physbc_lo_in(1:dm)
  physbc_hi(1:dm) = physbc_hi_in(1:dm)

  Interior   = Interior_in
  UserBC     = UserBC_in
  Inflow     = Inflow_in
  Outflow    = Outflow_in
  Symmetry   = Symmetry_in
  SlipWall   = SlipWall_in
  NoSlipWall = NoSlipWall_in

  coord_type = coord_type_in

  problo = ZERO
  probhi = ZERO
  center = ZERO

  problo(1:dm) = problo_in(1:dm)
  probhi(1:dm) = probhi_in(1:dm)
  center(1:dm) = center_in(1:dm)

  dg(:) = 1

  if (dim .lt. 2) then
     dg(2) = 0
  endif

  if (dim .lt. 3) then
     dg(3) = 0
  endif

end subroutine set_problem_params

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine set_grid_info(max_level_in, dx_level_in, domlo_in, domhi_in) &
     bind(C, name="set_grid_info")

  use prob_params_module, only: max_level, dx_level, domlo_level, domhi_level

  implicit none

  integer,          intent(in) :: max_level_in
  double precision, intent(in) :: dx_level_in(3*(max_level_in+1))
  integer,          intent(in) :: domlo_in(3*(max_level_in+1)), domhi_in(3*(max_level_in+1))

  integer :: lev, dir

  ! Sometimes this routine can get called multiple
  ! times upon initialization; in this case, just to
  ! be safe, we'll deallocate and start again.

  if (allocated(dx_level)) then
     deallocate(dx_level)
  endif
  if (allocated(domlo_level)) then
     deallocate(domlo_level)
  endif
  if (allocated(domhi_level)) then
     deallocate(domhi_level)
  endif

  max_level = max_level_in

  allocate(dx_level(1:3, 0:max_level))
  allocate(domlo_level(1:3, 0:max_level))
  allocate(domhi_level(1:3, 0:max_level))

  do lev = 0, max_level
     do dir = 1, 3
        dx_level(dir,lev) = dx_level_in(3*lev + dir)
        domlo_level(dir,lev) = domlo_in(3*lev + dir)
        domhi_level(dir,lev) = domhi_in(3*lev + dir)
     enddo
  enddo

end subroutine set_grid_info

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine clear_grid_info() &
     bind(C, name="clear_grid_info")

  use prob_params_module, only: max_level, dx_level, domlo_level, domhi_level

  implicit none

  deallocate(dx_level)
  deallocate(domlo_level)
  deallocate(domhi_level)

  max_level = 0

end subroutine clear_grid_info

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine pc_set_special_tagging_flag(dummy,flag) &
     bind(C, name="pc_set_special_tagging_flag")
  implicit none
  double precision :: dummy
  integer          :: flag
end subroutine pc_set_special_tagging_flag

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine get_tagging_params(name, namlen) &
     bind(C, name="get_tagging_params")

  use tagging_module

  implicit none
  ! Initialize the tagging parameters

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /tagging/ &
       denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
       enterr,     entgrad,   max_enterr_lev,   max_entgrad_lev, &
       velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
       vorterr,    max_vorterr_lev,                              &
       vfracerr,    max_vfracerr_lev,                              &
       presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
       temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev, &
       raderr,     radgrad,   max_raderr_lev,   max_radgrad_lev, &
       ftracerr, ftracgrad, max_ftracerr_lev,   max_ftracgrad_lev

  ! Set namelist defaults
  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = 10
  max_dengrad_lev = 10

  enterr = 1.d20
  entgrad = 1.d20
  max_enterr_lev = -1
  max_entgrad_lev = -1

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  velerr  = 1.d20
  velgrad = 1.d20
  max_velerr_lev = -1
  max_velgrad_lev = -1

  vorterr  = 1.d20
  max_vorterr_lev = -1

  vfracerr  = 1.d20
  max_vfracerr_lev = -1

  temperr  = 1.d20
  tempgrad = 1.d20
  max_temperr_lev = -1
  max_tempgrad_lev = -1

  raderr  = 1.d20
  radgrad = 1.d20
  max_raderr_lev = -1
  max_radgrad_lev = -1

  ftracerr  = 1.d20
  ftracgrad = 1.d20
  max_ftracerr_lev = -1
  max_ftracgrad_lev = -1

  ! create the filename
  if (namlen > maxlen) then
     call bl_error('probin file name too long')
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=tagging, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     call bl_error('ERROR: problem in the tagging namelist')
  endif

  close (unit=un)

end subroutine get_tagging_params

