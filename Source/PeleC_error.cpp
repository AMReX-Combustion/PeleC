
#include "PeleC.H"
#include "PeleC_error_F.H"

using std::string;
using namespace amrex;

void
PeleC::ErrorSetUp ()
{
  //
  // DEFINE ERROR ESTIMATION QUANTITIES
  //
  err_list.add("density",1,ErrorRec::Special,pc_denerror);
  err_list.add("Temp",1,ErrorRec::Special,pc_temperror);
  err_list.add("pressure",1,ErrorRec::Special,pc_presserror);
  err_list.add("x_velocity",1,ErrorRec::Special,pc_velerror);

#if (BL_SPACEDIM >= 2)
  err_list.add("y_velocity",1,ErrorRec::Special,pc_velerror);
#if (BL_SPACEDIM == 3)
  err_list.add("z_velocity",1,ErrorRec::Special,pc_velerror);
#endif
#endif

  err_list.add("magvort",1,ErrorRec::Special,pc_vorterror);

#ifdef PELE_USE_EB
  err_list.add("vfrac",1,ErrorRec::Special,pc_vfracerror);
#endif

  if (!flame_trac_name.empty())
  {
    int idx = -1;
    for (int i=0; i<spec_names.size(); ++i)
    {
      if (flame_trac_name == spec_names[i])
      {
        idx = i;
      }
    }

    if (idx >= 0)
    {
      const std::string name = "Y("+flame_trac_name+")";
      if (ParallelDescriptor::IOProcessor())
        std::cout << "Flame tracer will be " << name << '\n';
      err_list.add(name,1,ErrorRec::Special,pc_ftracerror);
    }
    else
    {
      amrex::Abort("Unknown species identified as flame_trac_name");
    }
  }
}
