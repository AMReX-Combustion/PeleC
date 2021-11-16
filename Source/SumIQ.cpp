#include <iomanip>

#include "PeleC.H"

void
PeleC::sum_integrated_quantities()
{
  BL_PROFILE("PeleC::sum_integrated_quantities()");

  if (verbose <= 0) {
    return;
  }

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real mass = 0.0;
  amrex::Real mom[3] = {0.0};
  amrex::Real rho_e = 0.0;
  amrex::Real rho_K = 0.0;
  amrex::Real rho_E = 0.0;
  amrex::Real enstr = 0.0;
  amrex::Real fuel_prod = 0;
  amrex::Real temp = 0;

  for (int lev = 0; lev <= finest_level; lev++) {
    PeleC& pc_lev = getLevel(lev);

    mass += pc_lev.volWgtSum("density", time, local_flag);
    mom[0] += pc_lev.volWgtSum("xmom", time, local_flag);
    mom[1] += pc_lev.volWgtSum("ymom", time, local_flag);
    mom[2] += pc_lev.volWgtSum("zmom", time, local_flag);
    rho_e += pc_lev.volWgtSum("rho_e", time, local_flag);
    rho_K += pc_lev.volWgtSum("kineng", time, local_flag);
    rho_E += pc_lev.volWgtSum("rho_E", time, local_flag);
    // enstr += pc_lev.volWgtSum("enstrophy", time, local_flag);

    if (!fuel_name.empty()) {
      fuel_prod += pc_lev.volWgtSum("rho_omega_" + fuel_name, time, local_flag);
    }

    temp += pc_lev.volWgtSum("Temp", time, local_flag);
  }

  if (verbose > 0) {
    const int nfoo = 10;
    amrex::Real foo[nfoo] = {mass,  mom[0], mom[1], mom[2],    rho_e,
                             rho_K, rho_E,  enstr,  fuel_prod, temp};
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealSum(
        foo, nfoo, amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor()) {
        int i = 0;
        mass = foo[i++];
        mom[0] = foo[i++];
        mom[1] = foo[i++];
        mom[2] = foo[i++];
        rho_e = foo[i++];
        rho_K = foo[i++];
        rho_E = foo[i++];
        enstr = foo[i++]; // NOLINT
        fuel_prod = foo[i++];
        temp = foo[i++];

        amrex::Print() << '\n';
        amrex::Print() << "TIME= " << time << " MASS        = " << mass << '\n';
        amrex::Print() << "TIME= " << time << " XMOM        = " << mom[0]
                       << '\n';
        amrex::Print() << "TIME= " << time << " YMOM        = " << mom[1]
                       << '\n';
        amrex::Print() << "TIME= " << time << " ZMOM        = " << mom[2]
                       << '\n';
        amrex::Print() << "TIME= " << time << " RHO*e       = " << rho_e
                       << '\n';
        amrex::Print() << "TIME= " << time << " RHO*K       = " << rho_K
                       << '\n';
        amrex::Print() << "TIME= " << time << " RHO*E       = " << rho_E
                       << '\n';
        // amrex::Print() << "TIME= " << time << " ENSTROPHY   = " << enstr
        //               << '\n';
        amrex::Print() << "TIME= " << time << " FUEL PROD   = " << fuel_prod
                       << '\n';
        amrex::Print() << "TIME= " << time << " TEMP        = " << temp << '\n';

        if (parent->NumDataLogs() > 0) {
          std::ostream& data_log1 = parent->DataLog(0);
          if (data_log1.good()) {
            const int datwidth = 14;
            if (time == 0.0) {
              data_log1 << std::setw(datwidth) << "          time";
              data_log1 << std::setw(datwidth) << "          mass";
              data_log1 << std::setw(datwidth) << "          xmom";
              data_log1 << std::setw(datwidth) << "          ymom";
              data_log1 << std::setw(datwidth) << "          zmom";
              data_log1 << std::setw(datwidth) << "         rho_K";
              data_log1 << std::setw(datwidth) << "         rho_e";
              data_log1 << std::setw(datwidth) << "         rho_E";
              // data_log1 << std::setw(datwidth) << "         enstr";
              data_log1 << std::setw(datwidth) << "     fuel_prod";
              data_log1 << std::setw(datwidth) << "          temp";
              data_log1 << std::endl;
            }

            // Write the quantities at this time
            const int datprecision = 6;
            data_log1 << std::setw(datwidth) << time;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << mass;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << mom[0];
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << mom[1];
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << mom[2];
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << rho_K;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << rho_e;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << rho_E;
            // data_log1 << std::setw(datwidth) <<
            // std::setprecision(datprecision)
            //          << enstr;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << fuel_prod;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << temp;
            data_log1 << std::endl;
          }
        }
      }
#ifdef AMREX_LAZY
    });
#endif
  }
}

void
PeleC::monitor_extrema()
{
  BL_PROFILE("PeleC::sum_integrated_quantities()");

  if (verbose <= 0) {
    return;
  }

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Vector<std::string> extrema_vars = {
    "density", "x_velocity", "y_velocity", "z_velocity",
    "eint_e",  "Temp",       "pressure",   "massfrac"};

  int nspec_extrema = 1;
  if (!fuel_name.empty()) {
    nspec_extrema += 1;
    extrema_vars.push_back(fuel_name);
  }
  if (!flame_trac_name.empty()) {
    nspec_extrema += 1;
    extrema_vars.push_back(flame_trac_name);
  }
  if (!extrema_spec_name.empty()) {
    nspec_extrema += 1;
    extrema_vars.push_back(extrema_spec_name);
  }
  auto nextrema = extrema_vars.size();
  constexpr amrex::Real neg_huge = std::numeric_limits<amrex::Real>::lowest();
  constexpr amrex::Real huge = std::numeric_limits<amrex::Real>::max();
  amrex::Vector<amrex::Real> extrema;
  extrema.resize(nextrema * 2);
  for (int ii = 0; ii < nextrema; ++ii) {
    extrema[ii] = neg_huge;
    extrema[nextrema + ii] = huge;
  }

  for (int lev = 0; lev <= finest_level; lev++) {
    PeleC& pc_lev = getLevel(lev);
    for (int ii = 0; ii < nextrema - nspec_extrema; ++ii) {
      extrema[ii] = amrex::max(
        extrema[ii], pc_lev.maxDerive(extrema_vars[ii], time, local_flag));
      extrema[nextrema + ii] = amrex::min(
        extrema[nextrema + ii],
        pc_lev.minDerive(extrema_vars[ii], time, local_flag));
    }

    // Handle species seperately
    auto mf = derive("massfrac", time, 0);
    BL_ASSERT(!(mf == nullptr));
    int idx_massfrac = nextrema - nspec_extrema;
    bool local = true;
    for (int ispec = 0; ispec < NUM_SPECIES; ispec++) {
      amrex::Real maxval = mf->max(ispec, 0, local);
      amrex::Real minval = mf->min(ispec, 0, local);
      // "massfrac" gets the extrema across all mass fractions
      extrema[idx_massfrac] = amrex::max(extrema[idx_massfrac], maxval);
      extrema[nextrema + idx_massfrac] =
        amrex::min(extrema[nextrema + idx_massfrac], minval);

      // Get values for any individual species if relevant
      for (int iext = 0; iext < nspec_extrema - 1; iext++) {
        int idx_spec = idx_massfrac + iext + 1;
        if (extrema_vars[idx_spec].compare(PeleC::spec_names[ispec]) == 0) {
          extrema[idx_spec] = amrex::max(extrema[idx_spec], maxval);
          extrema[nextrema + idx_spec] =
            amrex::min(extrema[nextrema + idx_spec], minval);
        }
      }
    }
  }

  if (verbose > 0) {
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(
        extrema.data(), nextrema * 2,
        amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor()) {

        amrex::Print() << std::endl;
        for (int ii = 0; ii < nextrema; ++ii) {
          const int datwidth = 22;
          const int datprecision = 14;
          amrex::Print() << "TIME= " << time;
          amrex::Print() << std::setw(datwidth) << extrema_vars[ii];
          amrex::Print() << " MIN= " << std::setw(datwidth)
                         << std::setprecision(datprecision)
                         << extrema[ii + nextrema];
          amrex::Print() << " MAX= " << std::setw(datwidth)
                         << std::setprecision(datprecision) << extrema[ii];
          amrex::Print() << std::endl;
        }

        if (parent->NumDataLogs() > 1) {
          std::ostream& data_log1 = parent->DataLog(1);
          if (data_log1.good()) {
            const int datwidth = 16;
            if (time == 0.0) {
              data_log1 << std::setw(datwidth) << "          time";
              for (int ii = 0; ii < nextrema; ++ii) {
                data_log1 << std::setw(datwidth - 4) << extrema_vars[ii]
                          << "-min";
                data_log1 << std::setw(datwidth - 4) << extrema_vars[ii]
                          << "-max";
              }
              data_log1 << std::endl;
            }

            // Write the quantities at this time
            const int datprecision = 8;
            data_log1 << std::setw(datwidth) << time;
            for (int ii = 0; ii < nextrema; ++ii) {
              data_log1 << std::setw(datwidth)
                        << std::setprecision(datprecision)
                        << extrema[ii + nextrema];
              data_log1 << std::setw(datwidth)
                        << std::setprecision(datprecision) << extrema[ii];
            }
            data_log1 << std::endl;
          }
        }
      }
#ifdef AMREX_LAZY
    });
#endif
  }
}
