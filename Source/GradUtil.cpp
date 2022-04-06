#include "GradUtil.H"

void
pc_compute_tangential_vel_derivs_eb(
  const amrex::Box& bx,
  const int dir,
  const amrex::Real dx1,
  const amrex::Real dx2,
  const EBBndryGeom* sv_ebg,
  const int Ncut,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<amrex::Real>& td)
{
  // dx1 and dx2 will be the trangential grid spacing
  const amrex::Real dx1inv = 1.0 / dx1;
  const amrex::Real dx2inv = 1.0 / dx2;
  const amrex::GpuArray<const amrex::Real, 3> weights{0.0, 1.0, 0.5};

  if (dir == 0) {
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      int i = 0;
      int j = 0;
      int k = 0;
      AMREX_D_TERM(i = sv_ebg[L].iv[0];, j = sv_ebg[L].iv[1];
                   , k = sv_ebg[L].iv[2];);
      const amrex::IntVect iv = amrex::IntVect{AMREX_D_DECL(i, j, k)};
      if (bx.contains(iv)) {
        const int jhip =
          j + static_cast<int>(flags(i, j, k).isConnected(0, 1, 0));
        const int jhim =
          j - static_cast<int>(flags(i, j, k).isConnected(0, -1, 0));
        const int jlop =
          j + static_cast<int>(flags(i - 1, j, k).isConnected(0, 1, 0));
        const int jlom =
          j - static_cast<int>(flags(i - 1, j, k).isConnected(0, -1, 0));
        const amrex::Real wjhi = weights[jhip - jhim];
        const amrex::Real wjlo = weights[jlop - jlom];
        td(i, j, k, 0) =
          0.5 * dx1inv *
          ((q(i, jhip, k, QU) - q(i, jhim, k, QU)) * wjhi +
           (q(i - 1, jlop, k, QU) - q(i - 1, jlom, k, QU)) * wjlo);
        td(i, j, k, 1) =
          0.5 * dx1inv *
          ((q(i, jhip, k, QV) - q(i, jhim, k, QV)) * wjhi +
           (q(i - 1, jlop, k, QV) - q(i - 1, jlom, k, QV)) * wjlo);
#if AMREX_SPACEDIM == 3
        td(i, j, k, 2) =
          0.5 * dx1inv *
          ((q(i, jhip, k, QW) - q(i, jhim, k, QW)) * wjhi +
           (q(i - 1, jlop, k, QW) - q(i - 1, jlom, k, QW)) * wjlo);

        const int khip =
          k + static_cast<int>(flags(i, j, k).isConnected(0, 0, 1));
        const int khim =
          k - static_cast<int>(flags(i, j, k).isConnected(0, 0, -1));
        const int klop =
          k + static_cast<int>(flags(i - 1, j, k).isConnected(0, 0, 1));
        const int klom =
          k - static_cast<int>(flags(i - 1, j, k).isConnected(0, 0, -1));
        const amrex::Real wkhi = weights[khip - khim];
        const amrex::Real wklo = weights[klop - klom];
        td(i, j, k, 3) =
          0.5 * dx2inv *
          ((q(i, j, khip, QU) - q(i, j, khim, QU)) * wkhi +
           (q(i - 1, j, klop, QU) - q(i - 1, j, klom, QU)) * wklo);
        td(i, j, k, 4) =
          0.5 * dx2inv *
          ((q(i, j, khip, QV) - q(i, j, khim, QV)) * wkhi +
           (q(i - 1, j, klop, QV) - q(i - 1, j, klom, QV)) * wklo);
        td(i, j, k, 5) =
          0.5 * dx2inv *
          ((q(i, j, khip, QW) - q(i, j, khim, QW)) * wkhi +
           (q(i - 1, j, klop, QW) - q(i - 1, j, klom, QW)) * wklo);
#endif
      }
    });
  } else if (dir == 1) {
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      int i = 0;
      int j = 0;
      int k = 0;
      AMREX_D_TERM(i = sv_ebg[L].iv[0];, j = sv_ebg[L].iv[1];
                   , k = sv_ebg[L].iv[2];);
      const amrex::IntVect iv = amrex::IntVect{AMREX_D_DECL(i, j, k)};
      if (bx.contains(iv)) {
        const int ihip =
          i + static_cast<int>(flags(i, j, k).isConnected(1, 0, 0));
        const int ihim =
          i - static_cast<int>(flags(i, j, k).isConnected(-1, 0, 0));
        const int ilop =
          i + static_cast<int>(flags(i, j - 1, k).isConnected(1, 0, 0));
        const int ilom =
          i - static_cast<int>(flags(i, j - 1, k).isConnected(-1, 0, 0));
        const amrex::Real wihi = weights[ihip - ihim];
        const amrex::Real wilo = weights[ilop - ilom];
        td(i, j, k, 0) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QU) - q(ihim, j, k, QU)) * wihi +
           (q(ilop, j - 1, k, QU) - q(ilom, j - 1, k, QU)) * wilo);
        td(i, j, k, 1) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QV) - q(ihim, j, k, QV)) * wihi +
           (q(ilop, j - 1, k, QV) - q(ilom, j - 1, k, QV)) * wilo);
#if AMREX_SPACEDIM == 3
        td(i, j, k, 2) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QW) - q(ihim, j, k, QW)) * wihi +
           (q(ilop, j - 1, k, QW) - q(ilom, j - 1, k, QW)) * wilo);

        const int khip =
          k + static_cast<int>(flags(i, j, k).isConnected(0, 0, 1));
        const int khim =
          k - static_cast<int>(flags(i, j, k).isConnected(0, 0, -1));
        const int klop =
          k + static_cast<int>(flags(i, j - 1, k).isConnected(0, 0, 1));
        const int klom =
          k - static_cast<int>(flags(i, j - 1, k).isConnected(0, 0, -1));
        const amrex::Real wkhi = weights[khip - khim];
        const amrex::Real wklo = weights[klop - klom];
        td(i, j, k, 3) =
          0.5 * dx2inv *
          ((q(i, j, khip, QU) - q(i, j, khim, QU)) * wkhi +
           (q(i, j - 1, klop, QU) - q(i, j - 1, klom, QU)) * wklo);
        td(i, j, k, 4) =
          0.5 * dx2inv *
          ((q(i, j, khip, QV) - q(i, j, khim, QV)) * wkhi +
           (q(i, j - 1, klop, QV) - q(i, j - 1, klom, QV)) * wklo);
        td(i, j, k, 5) =
          0.5 * dx2inv *
          ((q(i, j, khip, QW) - q(i, j, khim, QW)) * wkhi +
           (q(i, j - 1, klop, QW) - q(i, j - 1, klom, QW)) * wklo);
#endif
      }
    });
  } else if (dir == 2) {
    amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int L) {
      const int i = sv_ebg[L].iv[0];
      const int j = sv_ebg[L].iv[1];
      const int k = sv_ebg[L].iv[2];
      const amrex::IntVect iv = amrex::IntVect{AMREX_D_DECL(i, j, k)};
      if (bx.contains(iv)) {
        const int ihip =
          i + static_cast<int>(flags(i, j, k).isConnected(1, 0, 0));
        const int ihim =
          i - static_cast<int>(flags(i, j, k).isConnected(-1, 0, 0));
        const int ilop =
          i + static_cast<int>(flags(i, j, k - 1).isConnected(1, 0, 0));
        const int ilom =
          i - static_cast<int>(flags(i, j, k - 1).isConnected(-1, 0, 0));
        const amrex::Real wihi = weights[ihip - ihim];
        const amrex::Real wilo = weights[ilop - ilom];
        td(i, j, k, 0) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QU) - q(ihim, j, k, QU)) * wihi +
           (q(ilop, j, k - 1, QU) - q(ilom, j, k - 1, QU)) * wilo);
        td(i, j, k, 1) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QV) - q(ihim, j, k, QV)) * wihi +
           (q(ilop, j, k - 1, QV) - q(ilom, j, k - 1, QV)) * wilo);
        td(i, j, k, 2) =
          0.5 * dx1inv *
          ((q(ihip, j, k, QW) - q(ihim, j, k, QW)) * wihi +
           (q(ilop, j, k - 1, QW) - q(ilom, j, k - 1, QW)) * wilo);

        const int jhip =
          j + static_cast<int>(flags(i, j, k).isConnected(0, 1, 0));
        const int jhim =
          j - static_cast<int>(flags(i, j, k).isConnected(0, -1, 0));
        const int jlop =
          j + static_cast<int>(flags(i, j, k - 1).isConnected(0, 1, 0));
        const int jlom =
          j - static_cast<int>(flags(i, j, k - 1).isConnected(0, -1, 0));
        const amrex::Real wjhi = weights[jhip - jhim];
        const amrex::Real wjlo = weights[jlop - jlom];
        td(i, j, k, 3) =
          0.5 * dx2inv *
          ((q(i, jhip, k, QU) - q(i, jhim, k, QU)) * wjhi +
           (q(i, jlop, k - 1, QU) - q(i, jlom, k - 1, QU)) * wjlo);
        td(i, j, k, 4) =
          0.5 * dx2inv *
          ((q(i, jhip, k, QV) - q(i, jhim, k, QV)) * wjhi +
           (q(i, jlop, k - 1, QV) - q(i, jlom, k - 1, QV)) * wjlo);
        td(i, j, k, 5) =
          0.5 * dx2inv *
          ((q(i, jhip, k, QW) - q(i, jhim, k, QW)) * wjhi +
           (q(i, jlop, k - 1, QW) - q(i, jlom, k - 1, QW)) * wjlo);
      }
    });
  }
}
