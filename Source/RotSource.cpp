#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"

void
PeleC::construct_old_rot_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0;

  old_sources[rot_src]->setVal(0.0);

  if (!do_mrf && !do_srf) {
    return;
  }

  fill_rot_source(S_old, S_old, *old_sources[rot_src], ng);

  old_sources[rot_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_rot_source(amrex::Real /*time*/, amrex::Real /*dt*/)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[rot_src]->setVal(0.0);

  if (!do_mrf && !do_srf) {
    return;
  }

  fill_rot_source(S_old, S_new, *new_sources[rot_src], ng);
}

void
PeleC::fill_rot_source(
  const amrex::MultiFab& state_old
  /*unused*/,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& rot_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  amrex::Real r_mrf=mrf_rad;
  int mrf_on=do_mrf;
  int srf_on=do_srf;
  amrex::Real omega = (srf_on==1)?srf_omega:mrf_omega;
  amrex::Real axis  = (srf_on==1)?srf_axis:mrf_axis;
  
  auto prob_lo = geom.ProbLoArray();
  auto prob_hi = geom.ProbHiArray();
  const auto dx = geom.CellSizeArray();

  //need a check to prevent mrf in 2 and 1d
  //read a vector instead for mrf_axis_x and y
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> axis_loc={mrf_axis_x,mrf_axis_y,mrf_axis_z};

  auto const& sarrs = state_new.const_arrays();
  auto const& srcs = mrf_src.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
      mrf_src, ngs, axis, r_mrf, mrf_on, srf_on, omega,
      [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
          
          if (!flagarrs[nbx](i, j, k).isCovered()) 
          {
              const auto& sarr = sarrs[nbx];
              const auto& src = srcs[nbx];

              RealVect r(0.0,0.0,0.0);
              RealVect w(0.0,0.0,0.0);
              RealVect v(0.0,0.0,0.0);

              r[0]=prob_lo[0]+(i+0.5)*dx[0]-axis_loc[0];
              r[1]=prob_lo[1]+(j+0.5)*dx[1]-axis_loc[1];
              r[2]=prob_lo[2]+(k+0.5)*dx[2]-axis_loc[2];

              w[axis]=omega;

              if(mrf_on)
              {
                  t1dir=(axis+1)%AMREX_SPACEDIM;
                  t2dir=(axis+2)%AMREX_SPACEDIM;
                  rmag=std::sqrt(std::pow(r[t1dir],2.0)+std::pow(r[t2dir],2.0));
                  if(rmag < r_mrf)
                  {
                      v[0]=sarr(i, j, k, UMX)/sarr(i, j, k, URHO);
                      v[1]=sarr(i, j, k, UMY)/sarr(i, j, k, URHO);
                      v[2]=sarr(i, j, k, UMZ)/sarr(i, j, k, URHO);

                      RealVect w_cross_v = w.crossProduct(v); 
                      //multi-reference frame case
                      src(i, j, k, UMX) = -sarr(i,j,k,URHO)*w_cross_v[0];
                      src(i, j, k, UMY) = -sarr(i,j,k,URHO)*w_cross_v[1];
                      src(i, j, k, UMZ) = -sarr(i,j,k,URHO)*w_cross_v[2];
                  }
              }
              if(srf_on)
              {
                  v[0]=sarr(i, j, k, UMX)/sarr(i, j, k, URHO);
                  v[1]=sarr(i, j, k, UMY)/sarr(i, j, k, URHO);
                  v[2]=sarr(i, j, k, UMZ)/sarr(i, j, k, URHO);

                  //non-inertial frame case
                  RealVect w_cross_v = w.crossProduct(v); 
                  RealVect w_cross_r = w.crossProduct(r);
                  RealVect w_cross_w_cross_r=w.crossProduct(w_cross_r);
                  
                  src(i, j, k, UMX) = -sarr(i,j,k,URHO)*(2.0*w_cross_v[0]+w_cross_w_cross_r[0]);
                  src(i, j, k, UMY) = -sarr(i,j,k,URHO)*(2.0*w_cross_v[1]+w_cross_w_cross_r[1]);
                  src(i, j, k, UMZ) = -sarr(i,j,k,URHO)*(2.0*w_cross_v[2]+w_cross_w_cross_r[2]);
              }
          }
      });

  amrex::Gpu::synchronize();
}
