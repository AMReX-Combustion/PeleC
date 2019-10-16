
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <Transport_F.H>
#include <drag_F.H>

using namespace amrex;


void
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{

}


void
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{

} 

void
SprayParticleContainer::InitParticlesUniform(const int& lev, const int& num_ppc)
{

  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {   
      const Box& tile_box  = mfi.tilebox();

      auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                                       mfi.LocalTileIndex())];
    
      for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
	for (int i_part=0; i_part<num_ppc;i_part++) {

	  Real r = (rand()%100)/99.;
	  Real x = plo[0] + (iv[0] + r)*dx[0];

	  r = (rand()%100)/99.;
	  Real y = plo[1] + (iv[1] + r)*dx[1];

	  ParticleType p;
	  p.id()  = ParticleType::NextID();
	  p.cpu() = ParallelDescriptor::MyProc();
	  p.pos(0) = x;
	  p.pos(1) = y;
	  if (AMREX_SPACEDIM>2) {
	    r = (rand()%100)/99.;
	    Real z = plo[2] + (iv[2] + r)*dx[2];
	    p.pos(2) = z;
	  }

	  //Now assign the real data carried by each particle
	  p.rdata(0) = 0.;//u-vel
	  p.rdata(1) = 0.;//v-vel
	  if (AMREX_SPACEDIM>2) {
	    p.rdata(2) = 0.;//w-vel
	  }    
	  p.rdata(AMREX_SPACEDIM) = 293.; // temperature
	  p.rdata(AMREX_SPACEDIM+1) = 0.0200; // diameter
	  p.rdata(AMREX_SPACEDIM+2) = 0.68141; // fuel density

	  particles.push_back(p);

	}   
      }

    }
    
}
