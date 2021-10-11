
#include "SprayParticles.H"
#ifdef SPRAY_PELE_LM
#include <pelelm_prob.H>
#else
#include <PeleC.H>
#include "prob.H"
#endif
using namespace amrex;

IntVect
unflatten_particles(const Long idx, const IntVect& max_parts)
{
  IntVect indx;
  Long cidx = idx;
#if AMREX_SPACEDIM > 1
#if AMREX_SPACEDIM > 2
  indx[2] = cidx / (max_parts[0] * max_parts[1]);
  cidx -= indx[2] * max_parts[0] * max_parts[1];
#endif
  indx[1] = cidx / max_parts[0];
#endif
  indx[0] = cidx % max_parts[0];
  return indx;
}

bool
SprayParticleContainer::injectParticles(
  Real /*time*/,
  Real /*dt*/,
  int /*nstep*/,
  int /*lev*/,
  int /*finest_level*/,
  ProbParmHost const& /*prob_parm*/,
  ProbParmDevice const& /*prob_parm_d*/)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(
  ProbParmHost const& prob_parm, ProbParmDevice const& /*prob_parm_d*/)
{
  const int lev = 0;
  const int MyProc = ParallelDescriptor::MyProc();
  const int NProcs = ParallelDescriptor::NProcs();
  int NRedist = prob_parm.numRedist; // Number of times to redistribute
  // TODO: This might be overkill but issues persisted at high Summit node
  // counts
  if (NRedist < 0) {
    NRedist = 1;
    if (NProcs <= 512) {
      NRedist = 2;
    } else if (NProcs <= 1024) {
      NRedist = 4;
    } else if (NProcs <= 2048) {
      NRedist = 32;
    } else if (NProcs <= 4096) {
      NRedist = 48;
    }
  }
  Real part_dia = prob_parm.partDia;
  Real T_ref = prob_parm.partTemp;
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateY = m_sprayIndx.pstateY;
  const IntVect num_part = prob_parm.partNum;
  const RealVect part_vel = prob_parm.partVel;
  // Reference values for the particles
  Real part_vals[NAR_SPR + NSR_SPR];
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    part_vals[pstateVel + dir] = part_vel[dir];
  }
  part_vals[pstateT] = T_ref;
  part_vals[pstateDia] = part_dia;
  for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
    part_vals[pstateY + sp] = 0.;
  }
  part_vals[pstateY] = 1.; // Only use the first fuel species
  const Long total_part_num =
    AMREX_D_TERM(num_part[0], *num_part[1], *num_part[2]);
  const RealVect dx_part(AMREX_D_DECL(
    Geom(lev).ProbLength(0) / Real(num_part[0]),
    Geom(lev).ProbLength(1) / Real(num_part[1]),
    Geom(lev).ProbLength(2) / Real(num_part[2])));
  Long parts_pp = total_part_num / NProcs;
  // Number of particles per processor to be initialized
  Long cur_parts_pp = parts_pp;
  // Give any remaining particles to the last processor
  if (MyProc == NProcs - 1) {
    cur_parts_pp += (total_part_num % NProcs);
  }
  // Starting particle for this processor
  const Long first_part = MyProc * parts_pp;
  Gpu::HostVector<ParticleType> nparticles;
  Vector<Gpu::HostVector<Real>> nreals;
  if (NAR_SPR > 0) {
    nreals.resize(NAR_SPR);
  }
  for (int prc = 0; prc < cur_parts_pp; ++prc) {
    Long cur_part = first_part + prc;
    IntVect indx = unflatten_particles(cur_part, num_part);
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) = (Real(indx[dir]) + 0.5) * dx_part[dir];
    }
    for (int n = 0; n < NSR_SPR; ++n) {
      p.rdata(n) = part_vals[n];
    }
    for (int n = 0; n < NAR_SPR; ++n) {
      nreals[n].push_back(part_vals[n]);
    }
    nparticles.push_back(p);
  }
  ParticleLocData pld;
  // Only copy particle data for certain processors at a time
  int NRchunk = NProcs / NRedist;
  for (int nr = 0; nr < NRedist; ++nr) {
    std::map<std::pair<int, int>, Gpu::HostVector<ParticleType>> host_particles;
    std::map<std::pair<int, int>, std::array<Gpu::HostVector<Real>, NAR_SPR>>
      host_real_attribs;
    if (m_verbose > 0) {
      amrex::Print() << "Redistributing from processor " << nr * NRchunk
                     << " to " << (nr + 1) * NRchunk - 1 << '\n';
    }
    for (int which = nr * NRchunk; which < (nr + 1) * NRchunk; ++which) {
      if (which == MyProc) {
        while (!nparticles.empty()) {
          // Retrieve the last particle entry and add it to host_particles
          ParticleType& p = nparticles.back();
          Where(p, pld);
          std::pair<int, int> ind(pld.m_grid, pld.m_tile);
          host_particles[ind].push_back(p);
          for (int n = 0; n < NAR_SPR; ++n) {
            host_real_attribs[ind][n].push_back(nreals[n].back());
          }
          // Remove the particle just read
          nparticles.pop_back();
          for (int n = 0; n < NAR_SPR; ++n) {
            nreals[n].pop_back();
          }
        }
      } // if (which == MyProc)
    }   // for (int which ...
    for (auto& kv : host_particles) {
      auto grid = kv.first.first;
      auto tile = kv.first.second;
      const auto& src_tile = kv.second;
      auto& dst_tile = GetParticles(lev)[std::make_pair(grid, tile)];
      auto old_size = dst_tile.GetArrayOfStructs().size();
      auto new_size = old_size + src_tile.size();
      dst_tile.resize(new_size);

      // Copy the AoS part of the host particles to the GPU
      Gpu::copy(
        Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
        dst_tile.GetArrayOfStructs().begin() + old_size);
      for (int i = 0; i < NAR_SPR; ++i) {
        Gpu::copy(
          Gpu::hostToDevice,
          host_real_attribs[std::make_pair(grid, tile)][i].begin(),
          host_real_attribs[std::make_pair(grid, tile)][i].end(),
          dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
      }
    }
    Redistribute();
  } // for (int nr ...
  // Now copy over any remaining processors
  for (int which = NRedist * NRchunk; which < NProcs; ++which) {
    std::map<std::pair<int, int>, Gpu::HostVector<ParticleType>> host_particles;
    std::map<std::pair<int, int>, std::array<Gpu::HostVector<Real>, NAR_SPR>>
      host_real_attribs;
    if (m_verbose > 0) {
      amrex::Print() << "Redistributing from processor " << NRedist * NRchunk
                     << " to " << NProcs << '\n';
    }
    if (which == MyProc) {
      while (!nparticles.empty()) {
        // Retrieve the last particle entry and add it to host_particles
        ParticleType& p = nparticles.back();
        Where(p, pld);
        std::pair<int, int> ind(pld.m_grid, pld.m_tile);
        host_particles[ind].push_back(p);
        for (int n = 0; n < NAR_SPR; ++n) {
          host_real_attribs[ind][n].push_back(nreals[n].back());
        }
        // Remove the particle just read
        nparticles.pop_back();
        for (int n = 0; n < NAR_SPR; ++n) {
          nreals[n].pop_back();
        }
      }
    } // if (which == MyProc)
    for (auto& kv : host_particles) {
      auto grid = kv.first.first;
      auto tile = kv.first.second;
      const auto& src_tile = kv.second;
      auto& dst_tile = GetParticles(lev)[std::make_pair(grid, tile)];
      auto old_size = dst_tile.GetArrayOfStructs().size();
      auto new_size = old_size + src_tile.size();
      dst_tile.resize(new_size);

      // Copy the AoS part of the host particles to the GPU
      Gpu::copy(
        Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
        dst_tile.GetArrayOfStructs().begin() + old_size);
      for (int i = 0; i < NAR_SPR; ++i) {
        Gpu::copy(
          Gpu::hostToDevice,
          host_real_attribs[std::make_pair(grid, tile)][i].begin(),
          host_real_attribs[std::make_pair(grid, tile)][i].end(),
          dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
      }
    }
    Redistribute();
  } // for (int which ...
  Gpu::streamSynchronize();
}
