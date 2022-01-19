#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

#include "PeleC.H"
#include "prob.H"

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;

  AMREX_GPU_HOST
  constexpr explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm)
    : lprobparm(d_prob_parm)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    // const amrex::Real* prob_hi = geom.ProbHi();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    amrex::Real s_int[NVAR] = {0.0};
    amrex::Real s_ext[NVAR] = {0.0};

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]-1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif
      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]+1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM > 1
    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]-1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif
      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]+1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {

      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domlo[idir], n);
      }

#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]-1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif

      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domhi[idir], n);
      }
#ifdef PELEC_USE_TURBINFLOW
      if(iv[idir] == domlo[idir]+1){
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = dest(iv[0], iv[1], iv[2], n);
        }
      }
#endif
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#endif
#endif
  }
};

struct PCReactFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

void
pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  ProbParmDevice* probparmDD = PeleC::d_prob_parm_device; // probparm data for device
  ProbParmDevice* probparmDH = PeleC::h_prob_parm_device; // host copy of probparm data for device
  ProbParmHost* probparmH = PeleC::prob_parm_host;        // probparm data for host
  constexpr int dim = AMREX_SPACEDIM;

#ifdef PELEC_USE_TURBINFLOW
  if (probparmH->do_turb) {

    // Copy problem parameter structs to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);
    for (int dir=0; dir<dim; ++dir) {

      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(),dir) & bx);
      // if(bc[idir] == amrex::BCType::ext_dir)
      if (bcr[1].lo()[dir]==EXT_DIR && bndryBoxLO.ok())
      {
        //Create box with ghost cells and set them to zero
        amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
        int Grow = PeleC::numGrow();
        for(int n=0;n<dim;n++)
          growVect[n] = Grow;
        growVect[dir] = 0;
        amrex::Box modDom = geom.Domain();
        modDom.grow(growVect);
        auto bndryBoxLO_ghost = amrex::Box(amrex::adjCellLo(modDom,dir) & bx);
        data.setVal<amrex::RunOn::Host>(0.0,bndryBoxLO_ghost,UMX,dim);
        
      	add_turb(bndryBoxLO, data, 0, geom, time, dir, amrex::Orientation::low, probparmDH->tp);
        probparmDH->turb_ok[dir] = true;
      }

      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(),dir) & bx);
      if (bcr[1].hi()[dir]==EXT_DIR && bndryBoxHI.ok())
      {
        //Create box with ghost cells and set them to zero
        amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
        int Grow = PeleC::numGrow();
        for(int n=0;n<dim;n++)
          growVect[n] = Grow;
        growVect[dir] = 0;
        amrex::Box modDom = geom.Domain();
        modDom.grow(growVect);
        auto bndryBoxHI_ghost = amrex::Box(amrex::adjCellHi(modDom,dir) & bx);
        data.setVal<amrex::RunOn::Host>(0.0,bndryBoxHI_ghost,UMX,dim);

        add_turb(bndryBoxHI, data, 0, geom, time, dir, amrex::Orientation::high, probparmDH->tp);
        probparmDH->turb_ok[dir+dim] = true;
      }
    }

    // Copy problem parameter structs back to device
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
  }
#endif

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(PCHypFillExtDir{probparmDD});
  hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);

#ifdef PELEC_USE_TURBINFLOW
  if (probparmH->do_turb) {

    // Copy problem parameter structs to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

    for (int dir=0; dir<dim; ++dir) {
      if (probparmDH->turb_ok[dir]) {
        // probparmH->turbfab[dir].clear();
        probparmDH->turb_ok[dir] = false;
      }
      if (probparmDH->turb_ok[dir+dim]) {
        // probparmH->turbfab[dir+dim].clear();
        probparmDH->turb_ok[dir+dim] = false;
      }
    }

    // Copy problem parameter structs back to device
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
    
  }
#endif
}

void
pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  amrex::GpuBndryFuncFab<PCReactFillExtDir> react_bndry_func(
    PCReactFillExtDir{});
  react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void
pc_nullfill(
  amrex::Box const& /*bx*/,
  amrex::FArrayBox& /*data*/,
  const int /*dcomp*/,
  const int /*numcomp*/,
  amrex::Geometry const& /*geom*/,
  const amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcr*/,
  const int /*bcomp*/,
  const int /*scomp*/)
{
}
