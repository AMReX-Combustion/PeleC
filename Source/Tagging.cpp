#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "Tagging.H"

void
PeleC::read_tagging_params()
{
  amrex::ParmParse pp("tagging");

  pp.query("denerr", tagging_parm->denerr);
  pp.query("max_denerr_lev", tagging_parm->max_denerr_lev);
  pp.query("dengrad", tagging_parm->dengrad);
  pp.query("max_dengrad_lev", tagging_parm->max_dengrad_lev);

  pp.query("presserr", tagging_parm->presserr);
  pp.query("max_presserr_lev", tagging_parm->max_presserr_lev);
  pp.query("pressgrad", tagging_parm->pressgrad);
  pp.query("max_pressgrad_lev", tagging_parm->max_pressgrad_lev);

  pp.query("velerr", tagging_parm->velerr);
  pp.query("max_velerr_lev", tagging_parm->max_velerr_lev);
  pp.query("velgrad", tagging_parm->velgrad);
  pp.query("max_velgrad_lev", tagging_parm->max_velgrad_lev);

  pp.query("vorterr", tagging_parm->vorterr);
  pp.query("max_vorterr_lev", tagging_parm->max_vorterr_lev);

  pp.query("temperr", tagging_parm->temperr);
  pp.query("max_temperr_lev", tagging_parm->max_temperr_lev);
  pp.query("tempgrad", tagging_parm->tempgrad);
  pp.query("max_tempgrad_lev", tagging_parm->max_tempgrad_lev);

  pp.query("ftracerr", tagging_parm->ftracerr);
  pp.query("max_ftracerr_lev", tagging_parm->max_ftracerr_lev);
  pp.query("ftracgrad", tagging_parm->ftracgrad);
  pp.query("max_ftracgrad_lev", tagging_parm->max_ftracgrad_lev);

  pp.query("vfracerr", tagging_parm->vfracerr);
  pp.query("max_vfracerr_lev", tagging_parm->max_vfracerr_lev);
}
