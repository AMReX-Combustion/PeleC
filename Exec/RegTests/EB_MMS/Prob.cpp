#include <iomanip>

#include <PeleC.H>
#include <PeleC_F.H>
#include "Problem.H"
#include "Problem_F.H"

#ifdef USE_MASA
#include <masa.h>
using namespace MASA;
#endif

using namespace amrex;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
PeleC::problem_post_timestep(){

  if ((verbose <= 0) || (!do_mms)) return;

  bool local_flag = true;

  int finest_level   = parent->finestLevel();
  Real time          = state[State_Type].curTime();
  Real rho_mms_err = 0.0;
  Real u_mms_err = 0.0;
  Real v_mms_err = 0.0;
  Real w_mms_err = 0.0;
  Real p_mms_err = 0.0;
  Real rho_residual  = 0.0;
  Real rhou_residual = 0.0;
  Real rhov_residual = 0.0;
  Real rhow_residual = 0.0;
  Real rhoE_residual = 0.0;
  int datwidth     = 14;
  int datprecision = 6;

#ifdef USE_MASA
  if (level == 0)
  {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "... MMS problem post timestep" << std::endl;
    }

    // Calculate the errors and residuals
    for (int lev = 0; lev <= finest_level; lev++)
    {
      PeleC& pc_lev = getLevel(lev);

      rho_mms_err += pc_lev.volWgtSquaredSum("rhommserror", time, local_flag);
      u_mms_err += pc_lev.volWgtSquaredSum("ummserror", time, local_flag);
      v_mms_err += pc_lev.volWgtSquaredSum("vmmserror", time, local_flag);
      w_mms_err += pc_lev.volWgtSquaredSum("wmmserror", time, local_flag);
      p_mms_err += pc_lev.volWgtSquaredSum("pmmserror", time, local_flag);

      rho_residual += pc_lev.mms_volWgtSquaredSumDiff(Density, time, local_flag);
      rhou_residual += pc_lev.mms_volWgtSquaredSumDiff(Xmom, time, local_flag);
      rhov_residual += pc_lev.mms_volWgtSquaredSumDiff(Ymom, time, local_flag);
      rhow_residual += pc_lev.mms_volWgtSquaredSumDiff(Zmom, time, local_flag);
      rhoE_residual += pc_lev.mms_volWgtSquaredSumDiff(Eden, time, local_flag);
    }

    // Reductions
    ParallelDescriptor::ReduceRealSum(&rho_mms_err, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&u_mms_err, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&v_mms_err, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&w_mms_err, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&p_mms_err, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&rho_residual, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&rhou_residual, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&rhov_residual, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&rhow_residual, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&rhoE_residual, 1, ParallelDescriptor::IOProcessorNumber());

    // Get the norm and normalize it
    Real V = volume.sum(0, false);
    rho_mms_err = std::sqrt(rho_mms_err / V);
    u_mms_err = std::sqrt(u_mms_err / V);
    v_mms_err = std::sqrt(v_mms_err / V);
    w_mms_err = std::sqrt(w_mms_err / V);
    p_mms_err = std::sqrt(p_mms_err / V);
    rho_residual = std::sqrt(rho_residual / V);
    rhou_residual = std::sqrt(rhou_residual / V);
    rhov_residual = std::sqrt(rhov_residual / V);
    rhow_residual = std::sqrt(rhow_residual / V);
    rhoE_residual = std::sqrt(rhoE_residual/ V);

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "TIME= " << time << " RHO MMS ERROR  = " << rho_mms_err << '\n';
      std::cout << "TIME= " << time << " U MMS ERROR    = " << u_mms_err << '\n';
      std::cout << "TIME= " << time << " V MMS ERROR    = " << v_mms_err << '\n';
      std::cout << "TIME= " << time << " W MMS ERROR    = " << w_mms_err << '\n';
      std::cout << "TIME= " << time << " P MMS ERROR    = " << p_mms_err << '\n';
      std::cout << "TIME= " << time << " RHO RESIDUAL   = "   << rho_residual  << '\n';
      std::cout << "TIME= " << time << " RHO*U RESIDUAL = "   << rhou_residual << '\n';
      std::cout << "TIME= " << time << " RHO*V RESIDUAL = "   << rhov_residual << '\n';
      std::cout << "TIME= " << time << " RHO*W RESIDUAL = "   << rhow_residual << '\n';
      std::cout << "TIME= " << time << " RHO*E RESIDUAL = "   << rhoE_residual << '\n';

      if (parent->NumDataLogs() > 1 ) {

	std::ostream& data_log2 = parent->DataLog(1);

	// Write the quantities at this time
	data_log2 << std::setw(datwidth) <<  time;
	data_log2 << std::setw(datwidth) << std::setprecision(datprecision) << rho_mms_err;
	data_log2 << std::setw(datwidth) << std::setprecision(datprecision) << u_mms_err;
	data_log2 << std::setw(datwidth) << std::setprecision(datprecision) << v_mms_err;
	data_log2 << std::setw(datwidth) << std::setprecision(datprecision) << w_mms_err;
	data_log2 << std::setw(datwidth) << std::setprecision(datprecision) << p_mms_err;
	data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_residual;
	data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rhou_residual;
	data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rhov_residual;
	data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rhow_residual;
	data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rhoE_residual;
	data_log2 << std::endl;
      }
    }
  }

#else
    Error("MASA is not turned on. Turn on with USE_MASA=TRUE.");
#endif
}
#endif

#ifdef DO_PROBLEM_POST_INIT
void
PeleC::problem_post_init()
{

  if ((verbose <= 0) || (!do_mms)) return;

  Real time        = state[State_Type].curTime();
  int datwidth     = 14;
  int datprecision = 6;

  if (level == 0)
  {
    if (ParallelDescriptor::IOProcessor()) {

      if (parent->NumDataLogs() > 1 ) {

	std::ostream& data_log2 = parent->DataLog(1);
	if (time == 0.0) {
	  data_log2 << std::setw(datwidth) <<  "          time";
	  data_log2 << std::setw(datwidth) <<  "   rho_mms_err";
	  data_log2 << std::setw(datwidth) <<  "     u_mms_err";
	  data_log2 << std::setw(datwidth) <<  "     v_mms_err";
	  data_log2 << std::setw(datwidth) <<  "     w_mms_err";
	  data_log2 << std::setw(datwidth) <<  "     p_mms_err";
	  data_log2 << std::setw(datwidth) <<  "  rho_residual";
	  data_log2 << std::setw(datwidth) <<  " rhou_residual";
	  data_log2 << std::setw(datwidth) <<  " rhov_residual";
	  data_log2 << std::setw(datwidth) <<  " rhow_residual";
	  data_log2 << std::setw(datwidth) <<  " rhoE_residual";
	  data_log2 << std::endl;
	}
      }
    }
  }
}
#endif

#ifdef USE_MASA
Real
PeleC::mms_volWgtSquaredSumDiff (int comp,
				 Real time,
				 bool local)
{

  // Calculate volume weighted sum of the square of the difference
  // between the old and new quantity

  Real sum = 0.0;
  const Real* dx      = geom.CellSize();
  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);
  MultiFab diff(grids,dmap,1,0);

  // Calculate the difference between the states
  MultiFab::Copy(diff, S_old, comp, 0, 1, 0);
  MultiFab::Subtract(diff, S_new, comp, 0, 1, 0);

  if (level < parent->finestLevel())
  {
    const MultiFab& mask = getLevel(level+1).build_fine_mask();
    MultiFab::Multiply(diff, mask, 0, 0, 1, 0);
  }

#ifdef PELE_USE_EB
  MultiFab::Multiply(diff,vfrac,0,0,1,0);

  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
  for (MFIter mfi(diff,true); mfi.isValid(); ++mfi)
    {
      Real s = 0.0;
      const Box& box  = mfi.tilebox();
      const int* lo   = box.loVect();
      const int* hi   = box.hiVect();

#ifdef PELE_USE_EB
      const auto& flag_fab = flags[mfi];
      FabType typ = flag_fab.getType(box);
      if (typ == FabType::covered) {
        continue;
      }
#endif

      auto& fab = (diff)[mfi];

      pc_sumsquared(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
                    ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

      sum += s;
    }

  if (!local)
    ParallelDescriptor::ReduceRealSum(sum);

  return sum;

}
#endif
