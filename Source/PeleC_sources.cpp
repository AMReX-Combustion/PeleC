#include "PeleC.H"
#include "PeleC_F.H"

using namespace amrex;

void
PeleC::construct_old_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

    case ext_src:
        construct_old_ext_source(time, dt);
        break;

    case les_src:
        construct_old_les_source(time, dt, sub_iteration, sub_ncycle);
	break;

    case forcing_src:
        construct_old_forcing_source(time, dt);
        break;

#ifdef USE_MASA
    case mms_src:
        construct_old_mms_source(time);
        break;
#endif

    } // end switch
}

void
PeleC::construct_new_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

    case ext_src:
        construct_new_ext_source(time, dt);
        break;

    case les_src:
        construct_new_les_source(time, dt, sub_iteration, sub_ncycle);
	break;

    case forcing_src:
        construct_new_forcing_source(time, dt);
        break;

#ifdef USE_MASA
    case mms_src:
        construct_new_mms_source(time);
        break;
#endif

    } // end switch
}

// Obtain the sum of all source terms.

void
PeleC::sum_of_sources(MultiFab& source)
{

    int ng = source.nGrow();

    source.setVal(0.0);

    for (int n = 0; n < src_list.size(); ++n)
    {
      MultiFab::Add(source, *old_sources[src_list[n]], 0, 0, NUM_STATE, ng);
    }

    if (do_hydro)
    {
      MultiFab::Add(source, hydro_source, 0, 0, NUM_STATE, ng);
    }

    for (int n = 0; n < src_list.size(); ++n)
    {
      MultiFab::Add(source, *new_sources[src_list[n]], 0, 0, NUM_STATE, ng);
    }
}
