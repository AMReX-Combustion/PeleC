// clang-format off
#ifdef PELEC_USE_EB

#include <iamr_redistribution.H>
#include <AMReX_EB_slopes_K.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)

void
Redistribution::MakeITracker ( Box const& bx,
                               Array4<Real const> const& apx,
                               Array4<Real const> const& apy,
                               Array4<Real const> const& vfrac,
                               Array4<int> const& itracker,
                               Geometry const& lev_geom,
                               std::string redist_type)
{
#if 0
    int debug_verbose = 0;
#endif

    const Box domain = lev_geom.Domain();

    // Note that itracker has 4 components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    // We identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->

    Array<int,9> imap{0,-1,0,1,-1,1,-1,0,1};
    Array<int,9> jmap{0,-1,-1,-1,0,0,1,1,1};
    Array<int,9> inv_map{0,8,7,6,5,4,3,2,1};

    // We will use small_norm as an off just to break the tie when at 45 degrees ...
    const Real small_norm = 1.e-8;

    const auto& is_periodic_x = lev_geom.isPeriodic(0);
    const auto& is_periodic_y = lev_geom.isPeriodic(1);

    int preferred_direction = 0; // x-direction is preferred
    // int preferred_direction = 1; // y-direction is preferred

//  if (debug_verbose > 0)
//      amrex::Print() << " IN MAKE_ITRACKER DOING BOX " << bx << std::endl;

    Box const& bxg4 = amrex::grow(bx,4);

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        itracker(i,j,k,0) = 0;
    });

    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       if (vfrac(i,j,k) > 0.0 && vfrac(i,j,k) < 0.5)
       {
           Real apnorm, apnorm_inv;
           const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
           const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
           apnorm = std::sqrt(dapx*dapx+dapy*dapy);
           apnorm_inv = 1.0/apnorm;
           Real nx = dapx * apnorm_inv;
           Real ny = dapy * apnorm_inv;

           // We use small_norm as an offset just to break the tie when at 45 degrees ...

           if (preferred_direction == 1)
           {
               // y-direction is preferred
               if (nx > 0)
                  nx -= small_norm;
               else
                  nx += small_norm;

           } else {
               // x-direction is preferred
               if (ny > 0)
                  ny -= small_norm;
               else
                  ny += small_norm;
           }

           // As a first pass, choose just based on the normal
           if (std::abs(nx) > std::abs(ny))
           {
               if (nx > 0)
                   itracker(i,j,k,1) = 5;
               else
                   itracker(i,j,k,1) = 4;

           } else {
               if (ny > 0)
                   itracker(i,j,k,1) = 7;
               else
                   itracker(i,j,k,1) = 2;
           }

           bool xdir_mns_ok = (is_periodic_x || (i != domain.smallEnd(0)));
           bool xdir_pls_ok = (is_periodic_x || (i != domain.bigEnd(0)  ));
           bool ydir_mns_ok = (is_periodic_y || (j != domain.smallEnd(1)));
           bool ydir_pls_ok = (is_periodic_y || (j != domain.bigEnd(1)  ));

           // Override above logic if trying to reach outside a domain boundary (and non-periodic)
           if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
                (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
           {
               itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
           }
           if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
                (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
           {
               itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
           }

           // (i,j) merges with at least one cell now
           itracker(i,j,k,0) += 1;

           // (i+ioff,j+joff) is in the nbhd of (i,j)
           int ioff = imap[itracker(i,j,k,1)];
           int joff = jmap[itracker(i,j,k,1)];

           // Sanity check
           if (vfrac(i+ioff,j+joff,k) == 0.)
               amrex::Abort(" Trying to merge with covered cell");

           Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k);

#if 0
           if (debug_verbose > 0)
               amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                 " trying to merge with " << IntVect(i+ioff,j+joff) <<
                                 " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                                 " to get new sum_vol " <<  sum_vol << std::endl;
#endif

           // If the merged cell isn't large enough, we try to merge in the other direction
           if (sum_vol < 0.5)
           {
               // Original offset was in y-direction, so we will add to the x-direction
               // Note that if we can't because it would go outside the domain, we don't
               if (ioff == 0) {
                   if (nx >= 0 && xdir_pls_ok)
                   {
                       itracker(i,j,k,2) = 5;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (nx <= 0 && xdir_mns_ok)
                   {
                       itracker(i,j,k,2) = 4;
                       itracker(i,j,k,0) += 1;
                   }

               // Original offset was in x-direction, so we will add to the y-direction
               // Note that if we can't because it would go outside the domain, we don't
               } else {
                   if (ny >= 0 && ydir_pls_ok)
                   {
                       itracker(i,j,k,2) = 7;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (ny <= 0 && ydir_mns_ok)
                   {
                       itracker(i,j,k,2) = 2;
                       itracker(i,j,k,0) += 1;
                   }
               }

               if (itracker(i,j,k,0) > 1)
               {

                   // (i+ioff,j+joff) is in the nbhd of (i,j)
                   int ioff_l = imap[itracker(i,j,k,2)];
                   int joff_l = jmap[itracker(i,j,k,2)];

                   sum_vol += vfrac(i+ioff_l,j+joff_l,k);
#if 0
                   if (debug_verbose > 0)
                       amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                         " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
                                         " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                                          " to get new sum_vol " <<  sum_vol << std::endl;
#endif
               }
           }

           // Now we merge in the corner direction if we have already claimed two
           if (itracker(i,j,k,0) == 2)
           {
               // We already have two offsets, and we know they are in different directions
               ioff = imap[itracker(i,j,k,1)] + imap[itracker(i,j,k,2)];
               joff = jmap[itracker(i,j,k,1)] + jmap[itracker(i,j,k,2)];

               if (ioff > 0 and joff > 0)
                   itracker(i,j,k,3) = 8;
               else if (ioff < 0 and joff > 0)
                   itracker(i,j,k,3) = 6;
               else if (ioff > 0 and joff < 0)
                   itracker(i,j,k,3) = 3;
               else
                   itracker(i,j,k,3) = 1;

               // (i,j) merges with at least three cells now
               itracker(i,j,k,0) += 1;

               sum_vol += vfrac(i+ioff,j+joff,k);
#if 0
               if (debug_verbose > 0)
                   amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                     " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
                                     " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                                     " to get new sum_vol " <<  sum_vol << std::endl;
#endif
           }
           if (sum_vol < 0.5)
           {
#if 0
             amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
                               IntVect(i,j) << " so stuck with sum_vol " << sum_vol << std::endl;
#endif
             amrex::Abort("Couldnt merge with enough cells to raise volume greater than 0.5");
           }
       }
    });

    if (redist_type == "State") return;

    // At this point every cell knows who it wants to merge with, but
    //   (1) not who wants to merge with it
    //   (2) not who its neighbor also wants to merge with
    // In this loop we only address (1)
    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       // Here we don't test on vfrac because some of the neighbors are full cells
       // We test on whether this cell has any neighbors
       if (itracker(i,j,k,0) > 0)
       {
          for (int ipair = 1; ipair <= itracker(i,j,k,0); ipair++)
          {
               // (i+ioff,j+joff) is in the nbhd of (i,j)
               int ioff = imap[itracker(i,j,k,ipair)];
               int joff = jmap[itracker(i,j,k,ipair)];

               if (bxg4.contains(IntVect(i+ioff,j+joff)))
               {
                   int n_of_nbor = itracker(i+ioff,j+joff,k,0);
                   bool found = false;
                   for (int ipair_nbor = 1; ipair_nbor <= n_of_nbor; ipair_nbor++)
                   {
                       if (imap[itracker(i+ioff,j+joff,k,ipair_nbor)] + ioff == 0 &&
                           jmap[itracker(i+ioff,j+joff,k,ipair_nbor)] + joff == 0)
                           found = true;
                   }

                   if  (!found)
                   {
                       // My neigbor didn't know about me so add me to my nbor's list of neighbors
                       itracker(i+ioff,j+joff,k,0) += 1;
                       itracker(i+ioff,j+joff,k,n_of_nbor+1) = inv_map[itracker(i,j,k,ipair)];
#if 0
                       if (debug_verbose > 1)
                           amrex::Print() << "Cell   " << IntVect(i,j) << " had nbor " << IntVect(i+ioff,j+joff)
                                          << " in its nbor list by taking inverse of " << itracker(i,j,k,ipair)
                                          << " which gave " << inv_map[itracker(i,j,k,ipair)] << std::endl;
                       if (debug_verbose > 1)
                           amrex::Print() << "Adding " << IntVect(i+ioff+imap[itracker(i+ioff,j+joff,k,n_of_nbor+1)],
                                                                  j+joff+jmap[itracker(i+ioff,j+joff,k,n_of_nbor+1)])
                                          << " to the nbor list of " << IntVect(i+ioff,j+joff) << std::endl;
#endif
                   } // found
               } // bxg4 contains
          } // ipair
       } // itracker
    });

    // Here we address (2), i.e. we want the neighbor of my neighbor to be my neighbor
    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       // Test on whether this cell has any neighbors and not already all three neighbors
       if (itracker(i,j,k,0) > 0 && itracker(i,j,k,0) < 3)
       {
          // Loop over my neighbors
          for (int ipair = 1; ipair <= itracker(i,j,k,0); ipair++)
          {
               // (i_n,j_n) is in the nbhd of (i,j)
               int i_n = i + imap[itracker(i,j,k,ipair)];
               int j_n = j + jmap[itracker(i,j,k,ipair)];

               if (bxg4.contains(IntVect(i_n,j_n)))
               {

                 // Loop over the neighbors of my neighbors
                 // If any of these aren't already my neighbor, make them my neighbor
                 int ipair_n = 1;
                 while (ipair_n <= itracker(i_n,j_n,k,0))
                 {
                    // (i_nn,j_nn) is in the nbhd of (i_n,j_n)
                    int i_nn = i_n + imap[itracker(i_n,j_n,k,ipair_n)];
                    int j_nn = j_n + jmap[itracker(i_n,j_n,k,ipair_n)];

                    bool found = false;

                    // Is this nbor of my nbor already my nbor (or me)??
                    for (int ipair_2 = 1; ipair_2 <= itracker(i,j,k,0); ipair_2++)
                    {
                        // (i_n2,j_n2) is in the nbhd of (i,j)
                        int i_n2 = i + imap[itracker(i,j,k,ipair_2)];
                        int j_n2 = j + jmap[itracker(i,j,k,ipair_2)];
                        if ( (i_nn == i_n2 && j_nn == j_n2) or (i_nn == i && j_nn == j) )
                            found = true;
                    }
#if 0
                    if (debug_verbose > 1)
                    {
                        if (!found)
                            amrex::Print() << "DOING CELL " << IntVect(i,j) << " who has nbor " << IntVect(i_n,j_n) <<
                                              " who has nbor " << IntVect(i_nn,j_nn) << " which was NOT found " << std::endl;
                        else
                            amrex::Print() << "DOING CELL " << IntVect(i,j) << " who has nbor " << IntVect(i_n,j_n) <<
                                              " who has nbor " << IntVect(i_nn,j_nn) << " which was found " << std::endl;
                    }
#endif

                    if (!found)
                    {
                        // My neighbor had a neighbor I didn't know so adding it here
                        itracker(i,j,k,0) += 1;
                        int n_nbor = itracker(i,j,k,0);
                        if (j_nn-j < 0)
                           itracker(i,j,k,n_nbor) = (i_nn-i)+2; // short-cut for mapping onto 1,2 or 3
                        else if (j_nn-j > 0)
                           itracker(i,j,k,n_nbor) = (i_nn-i)+7; // short-cut for mapping onto 6,7, or 8
                        else if (i_nn-i > 0)
                           itracker(i,j,k,n_nbor) = 5;
                        else
                           itracker(i,j,k,n_nbor) = 4;

#if 0
                        if (debug_verbose > 1)
                            amrex::Print() << "Adding " << IntVect(i_nn,j_nn) << " to the nbor list of " << IntVect(i,j) <<
                                              " in component " << n_nbor << std::endl;
                        if (debug_verbose > 1)
                            amrex::Print() << "Sanity check -- these should be the same: " << IntVect(i_nn,j_nn) << " " <<
                                IntVect(i+imap[itracker(i,j,k,n_nbor)],j+jmap[itracker(i,j,k,n_nbor)]) << std::endl;
#endif
                    }
                    ipair_n++;
                 } // while
               } // bxg4 contains
          } // ipair
       } // itracker
    });

}
#endif
#endif
// clang-format on
