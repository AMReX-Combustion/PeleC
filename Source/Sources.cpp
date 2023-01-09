#include "PeleC.H"
#include "IndexDefines.H"

void
PeleC::construct_old_source(
  int src,
  amrex::Real time,
  amrex::Real dt,
  int /*amr_iteration*/,
  int amr_ncycle,
  int sub_iteration,
  int sub_ncycle)
{
  AMREX_ASSERT(src >= 0 && src < num_src);

#ifndef PELEC_USE_SPRAY
  amrex::ignore_unused(amr_ncycle);
#endif
  switch (src) {

  case ext_src:
    construct_old_ext_source(time, dt);
    break;

  case les_src:
    construct_old_les_source(time, dt, sub_iteration, sub_ncycle);
    break;

  case forcing_src:
    construct_old_forcing_source(time, dt);
    break;

#ifdef PELEC_USE_SPRAY
  case spray_src:
    particleMKD(time, dt, sub_iteration, sub_ncycle, amr_ncycle);
    break;
#endif

#ifdef PELEC_USE_SOOT
  case soot_src:
    construct_old_soot_source(time, dt);
    break;
#endif

#ifdef PELEC_USE_MASA
  case mms_src:
    construct_old_mms_source(time);
    break;
#endif
  } // end switch
}

void
PeleC::construct_new_source(
  int src,
  amrex::Real time,
  amrex::Real dt,
  int /*amr_iteration*/,
  int amr_ncycle,
  int sub_iteration,
  int sub_ncycle)
{
  AMREX_ASSERT(src >= 0 && src < num_src);

#ifndef PELEC_USE_SPRAY
  amrex::ignore_unused(amr_ncycle);
#endif
  switch (src) {

  case ext_src:
    construct_new_ext_source(time, dt);
    break;

  case les_src:
    construct_new_les_source(time, dt, sub_iteration, sub_ncycle);
    break;

  case forcing_src:
    construct_new_forcing_source(time, dt);
    break;
#ifdef PELEC_USE_SPRAY
  case spray_src:
    particleMK(time, dt, sub_iteration, sub_ncycle, amr_ncycle);
    break;
#endif

#ifdef PELEC_USE_SOOT
  case soot_src:
    construct_new_soot_source(time, dt);
    break;
#endif

#ifdef PELEC_USE_MASA
  case mms_src:
    construct_new_mms_source(time);
    break;
#endif
  } // end switch
}

// Obtain the sum of all source terms.
void
PeleC::sum_of_sources(amrex::MultiFab& source)
{
  int ng = source.nGrow();

  source.setVal(0.0);

  for (int n = 0; n < src_list.size(); ++n) {
    amrex::MultiFab::Add(source, *old_sources[src_list[n]], 0, 0, NVAR, ng);
  }

  if (do_hydro) {
    amrex::MultiFab::Add(source, hydro_source, 0, 0, NVAR, ng);
  }

  for (int n = 0; n < src_list.size(); ++n) {
    amrex::MultiFab::Add(source, *new_sources[src_list[n]], 0, 0, NVAR, ng);
  }
}
