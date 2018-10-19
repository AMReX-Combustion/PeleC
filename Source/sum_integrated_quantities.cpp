#include <iomanip>

#include <PeleC.H>
#include <PeleC_F.H>

using namespace amrex;

void
PeleC::sum_integrated_quantities ()
{
    BL_PROFILE("PeleC::sum_integrated_quantities()");

    if (verbose <= 0) return;

    bool local_flag = true;

    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real mom[3]      = { 0.0 };
    Real rho_e       = 0.0;
    Real rho_K       = 0.0;
    Real rho_E       = 0.0;
    Real enstr       = 0.0;
    Real fuel_prod   = 0;
    Real temp        = 0;
    int datwidth     = 14;
    int datprecision = 6;

    for (int lev = 0; lev <= finest_level; lev++)
    {
      PeleC& pc_lev = getLevel(lev);

      mass   += pc_lev.volWgtSum("density", time, local_flag);
      mom[0] += pc_lev.volWgtSum("xmom", time, local_flag);
      mom[1] += pc_lev.volWgtSum("ymom", time, local_flag);
      mom[2] += pc_lev.volWgtSum("zmom", time, local_flag);
      rho_e += pc_lev.volWgtSum("rho_e", time, local_flag);
      rho_K += pc_lev.volWgtSum("kineng", time, local_flag);
      rho_E += pc_lev.volWgtSum("rho_E", time, local_flag);
      enstr += pc_lev.volWgtSum("enstrophy", time, local_flag);

      if (fuel_name != "")
      {
	fuel_prod += pc_lev.volWgtSum("rho_omega_" + fuel_name, time, local_flag);
      }

      temp += pc_lev.volWgtSum("Temp", time, local_flag);
    }
 
    if (verbose > 0)
    {

        const int nfoo = 10;

	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], rho_e, rho_K, rho_E, enstr, fuel_prod, temp};
#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

	ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

	if (ParallelDescriptor::IOProcessor()) {

	    int i = 0;
	    mass       = foo[i++];
	    mom[0]     = foo[i++];
            mom[1]     = foo[i++];
            mom[2]     = foo[i++];
	    rho_e      = foo[i++];
	    rho_K      = foo[i++];
            rho_E      = foo[i++];
	    enstr      = foo[i++];
            fuel_prod  = foo[i++];
            temp       = foo[i++];

	    std::cout << '\n';
	    std::cout << "TIME= " << time << " MASS        = "   << mass      << '\n';
	    std::cout << "TIME= " << time << " XMOM        = "   << mom[0]    << '\n';
	    std::cout << "TIME= " << time << " YMOM        = "   << mom[1]    << '\n';
	    std::cout << "TIME= " << time << " ZMOM        = "   << mom[2]    << '\n';
	    std::cout << "TIME= " << time << " RHO*e       = "   << rho_e     << '\n';
	    std::cout << "TIME= " << time << " RHO*K       = "   << rho_K     << '\n';
	    std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
	    std::cout << "TIME= " << time << " ENSTROPHY   = "   << enstr     << '\n';
	    std::cout << "TIME= " << time << " FUEL PROD   = "   << fuel_prod << '\n';
	    std::cout << "TIME= " << time << " TEMP        = "   << temp      << '\n';

	    if (parent->NumDataLogs() > 0 ) {

	       std::ostream& data_log1 = parent->DataLog(0);

	       if (data_log1.good()) {

		  if (time == 0.0) {
		      data_log1 << std::setw(datwidth) <<  "          time";
		      data_log1 << std::setw(datwidth) <<  "          mass";
		      data_log1 << std::setw(datwidth) <<  "          xmom";
		      data_log1 << std::setw(datwidth) <<  "          ymom";
		      data_log1 << std::setw(datwidth) <<  "          zmom";
		      data_log1 << std::setw(datwidth) <<  "         rho_K";
		      data_log1 << std::setw(datwidth) <<  "         rho_e";
		      data_log1 << std::setw(datwidth) <<  "         rho_E";
		      data_log1 << std::setw(datwidth) <<  "         enstr";
		      data_log1 << std::setw(datwidth) <<  "     fuel_prod";
		      data_log1 << std::setw(datwidth) <<  "          temp";
		      data_log1 << std::endl;
		  }

		  // Write the quantities at this time
		  data_log1 << std::setw(datwidth) <<  time;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mass;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[0];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[1];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[2];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_K;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_e;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_E;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << enstr;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << fuel_prod;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << temp;
		  data_log1 << std::endl;

	       }

	    }
	}
#ifdef BL_LAZY
	});
#endif
    }
}
