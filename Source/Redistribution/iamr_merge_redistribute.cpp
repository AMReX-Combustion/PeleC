// clang-format off
#ifdef PELEC_USE_EB

#include <iamr_redistribution.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

void
Redistribution::MergeRedistribute ( Box const& bx, int ncomp,
                                    Array4<Real> const& dUdt_out,
                                    Array4<Real> const& dUdt_in,
                                    Array4<Real const> const& vfrac,
                                    Array4<int> const& itracker,
                                    Geometry const& lev_geom)
{
    bool debug_print = false;

    const Box domain = lev_geom.Domain();

    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    //
#if (AMREX_SPACEDIM == 2)
    Array<int,9> imap{0,-1, 0, 1,-1, 1,-1, 0, 1};
    Array<int,9> jmap{0,-1,-1,-1, 0, 0, 1, 1, 1};
    Array<int,9> kmap{0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
    Array<int,27>    imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    Array<int,27>    jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    Array<int,27>    kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

//  if (debug_print)
//      amrex::Print() << " IN MERGE_REDISTRIBUTE DOING BOX " << bx << " with ncomp " << ncomp << std::endl;

    const Real small_norm = 1.e-8;

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0)
        {
            for (int n = 0; n < ncomp; n++)
                dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
        } else {
            // We shouldn't need to do this but just in case ...
            for (int n = 0; n < ncomp; n++)
                dUdt_out(i,j,k,n) = 1.e100;
        }
    });

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       if (vfrac(i,j,k) > 0.0)
       {
           if (itracker(i,j,k,0) > 0)
           {
               for (int n = 0; n < ncomp; n++)
               {
                   Real sum_vol = vfrac(i,j,k);
                   Real sum_upd = vfrac(i,j,k) * dUdt_in(i,j,k,n);
                   for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
                   {
                       sum_upd +=   vfrac(i+imap[itracker(i,j,k,i_nbor)],
                                          j+jmap[itracker(i,j,k,i_nbor)],
                                          k+kmap[itracker(i,j,k,i_nbor)]) *
                                  dUdt_in(i+imap[itracker(i,j,k,i_nbor)],
                                          j+jmap[itracker(i,j,k,i_nbor)],
                                          k+kmap[itracker(i,j,k,i_nbor)],n);
                       sum_vol +=   vfrac(i+imap[itracker(i,j,k,i_nbor)],
                                          j+jmap[itracker(i,j,k,i_nbor)],
                                          k+kmap[itracker(i,j,k,i_nbor)]);
                   }

                   if (sum_vol < 0.5)
                   {
//                    amrex::Print() << "SUM_VOL STILL TOO SMALL at " <<
//                        IntVect(AMREX_D_DECL(i,j,k)) << " " << sum_vol << std::endl;
                      amrex::Abort();
                   }

                   Real avg_update = sum_upd / sum_vol;

                   dUdt_out(i,j,k,n) = avg_update;
               }
           }
       }
    });

    //
    // This tests whether the redistribution procedure was conservative -- it is only relevant if
    //      the box covers the entire domain
    //
#if 0
    { // STRT:SUM OF FINAL DUDT
        for (int n = 0; n < ncomp; n++)
        {
          Real sum1(0);
          Real sum2(0);
#if (AMREX_SPACEDIM == 2)
          int k = 0;
#else
          for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
#endif
           for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
           {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
            {
              if (vfrac(i,j,k) > 0.)
              {
                  sum1 += vfrac(i,j,k)*dUdt_in(i,j,k,n);
                  sum2 += vfrac(i,j,k)*dUdt_out(i,j,k,n);
              }
            }
          }

         if (std::abs(sum1-sum2) > 1.e-8 * sum1 && std::abs(sum1-sum2) > 1.e-8)
         {
//          amrex::Print() << " TESTING COMPONENT " << n << std::endl;
//          amrex::Print() << " SUMS DO NOT MATCH " << sum1 << " " << sum2 << std::endl;
            amrex::Abort(0);
         }
        }
    } //  END:SUM OF FINAL DUDT
#endif
}
#endif
// clang-format on
