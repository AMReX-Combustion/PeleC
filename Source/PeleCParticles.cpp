#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "PeleC.H"
#include "PeleC_F.H"
#include "Spray_F.H"
#include <chemistry_file.H>

using namespace amrex;

#ifdef AMREX_PARTICLES

namespace
{
    bool virtual_particles_set = false;

    std::string ascii_particle_file;
    std::string binary_particle_file;

    const std::string chk_particle_file("Spray");

    //
    // We want to call this routine when on exit to clean up particles.
    //

    //
    // Containers for the real "active" Particles
    //
    SprayParticleContainer* SprayPC = 0;
    //
    // Container for temporary, virtual Particles
    //
    SprayParticleContainer* VirtPC  = 0;
    //
    // Container for temporary, ghost Particles
    //
    SprayParticleContainer* GhostPC  = 0;

    void RemoveParticlesOnExit ()
    {
        delete SprayPC;
        delete GhostPC;
        delete VirtPC;
    }
}

int PeleC::do_spray_particles          =  0;
int PeleC::particle_verbose            =  0;
Real PeleC::particle_cfl               = 0.05;
int PeleC::particle_mass_tran = 0;
int PeleC::particle_heat_tran = 0;
int PeleC::particle_mom_tran = 0;
int PeleC::particle_nfuel_species = 1;
Vector<Real> PeleC::partf_mf;
Vector<Real> PeleC::partf_crit_T;
Vector<Real> PeleC::partf_boil_T;
Vector<Real> PeleC::partf_cp;
Vector<Real> PeleC::partf_latent;
Vector<Real> PeleC::partf_molwt;
Vector<std::string> PeleC::partf_fuel_names;
Real PeleC::partf_ref_T;
Vector<int> PeleC::partf_indx;

namespace {
    std::string       particle_init_file;
    int               particle_init_uniform = 0;
    std::string       timestamp_dir;
    std::vector<int>  timestamp_indices;
}

SprayParticleContainer*
PeleC::theSprayPC()
{
    return SprayPC;
}

SprayParticleContainer*
PeleC::theVirtPC ()
{
    return VirtPC;
}

SprayParticleContainer*
PeleC::theGhostPC ()
{
    return GhostPC;
}

void
PeleC::particle_est_time_step (Real& est_dt)
{
    BL_PROFILE("PeleC::particle_est_time_step()");
    Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

    if (est_dt_particle > 0)
        est_dt = std::min(est_dt, est_dt_particle);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        if (est_dt_particle > 0)
        {
            std::cout << "...estdt from particles at level "
                      << level << ": " << est_dt_particle << '\n';
        }
        else
        {
            std::cout << "...there are no particles at level "
                      << level << '\n';
        }
    }
}

void 
PeleC::read_particle_params ()
{
  if (verbose)
    {
      Print() << "reading particle parameters" << std::endl;
    }
    ParmParse pp("pelec");
    
    pp.query("do_spray_particles",do_spray_particles);

    ParmParse ppp("particles");
    //
    // Control the verbosity of the Particle class
    //
    ppp.query("v",particle_verbose);

    ppp.get("mass_transfer",particle_mass_tran);
    ppp.get("heat_transfer",particle_heat_tran);
    ppp.get("mom_transfer",particle_mom_tran);
    //
    // Used in import_fuel_properties()
    //
    particle_nfuel_species = ppp.countval("fuel_species");
    const int nfuel = particle_nfuel_species;
    partf_mf.resize(nfuel);
    partf_crit_T.resize(nfuel);
    partf_boil_T.resize(nfuel);
    partf_latent.resize(nfuel);
    partf_cp.resize(nfuel);
    partf_molwt.resize(nfuel);
    partf_indx.assign(nfuel, -1);
    partf_fuel_names.assign(nfuel, "");
    // Read in the mass fractions of the particles from the input file and ensure they sum to 1
    Real total_fuel_mf = 0.;
    for (int i = 0; i != nfuel; ++i)
      {
	ppp.getkth("fuel_species", i, partf_fuel_names[i], 0);
	ppp.getkth("fuel_mass_frac", i, partf_mf[i], 0);
	total_fuel_mf += partf_mf[i];
      }
    if (std::abs(total_fuel_mf - 1.) > 1.E-8)
      {
	Abort("Particle fuel mass fractions do not sum to 1");
      }
    // TODO: Would be nice if these were tabulated or computed and not input directly
    ppp.getarr("fuel_crit_temp", partf_crit_T, 0, nfuel);
    ppp.getarr("fuel_boil_temp", partf_boil_T, 0, nfuel);
    ppp.getarr("fuel_latent", partf_latent, 0, nfuel);
    ppp.getarr("fuel_cp", partf_cp, 0, nfuel);
    // Must use same reference temperature for all fuels
    // TODO: This means the reference temperature must be the same for all fuel species
    ppp.get("fuel_ref_temp", partf_ref_T);
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in initData() on startup to set a uniform particle field
    //
    ppp.query("particle_init_uniform", particle_init_uniform);
    //
    // Used in post_restart() to read in a file of particles.
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    // ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    //ppp.query("particle_output_file", particle_output_file);
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(timestamp_dir, 0755))
            amrex::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
}

void
PeleC::define_particles()
{
  // Fluid species mass fractions
  Vector<Real> species_mw(NumSpec);
  get_mw(species_mw.dataPtr());
  for (int i = 0; i != particle_nfuel_species; ++i)
    {
      for (int ns = 0; ns != NumSpec; ++ns)
	{
	  std::string gas_spec = spec_names[ns];
	  if (gas_spec == partf_fuel_names[i])
	    {
	      partf_indx[i] = ns;
	      partf_molwt[i] = species_mw[ns];
	    }
	}
      if (partf_indx[i] < 0)
	{
	  Abort("Particle fuel species not found in fluid species list");
	}
    }
}

void
PeleC::setup_virtual_particles()
{
    BL_PROFILE("PeleC::setup_virtual_particles()");

    if(PeleC::theSprayPC() != 0 && !virtual_particles_set)
    {
        SprayParticleContainer::AoS virts;
        if (level < parent->finestLevel())
        {
            ((PeleC*) &parent->getLevel(level+1))->setup_virtual_particles();
            PeleC::theVirtPC()->CreateVirtualParticles(level+1, virts);
            PeleC::theVirtPC()->AddParticlesAtLevel(virts, level);

            PeleC::theSprayPC()->CreateVirtualParticles(level+1, virts);
            PeleC::theVirtPC()->AddParticlesAtLevel(virts, level);
        }
        virtual_particles_set = true;
    }
}

void
PeleC::remove_virtual_particles()
{
    BL_PROFILE("PeleC::remove_virtual_particles()");
    if (VirtPC != 0)
        VirtPC->RemoveParticlesAtLevel(level);
    virtual_particles_set = false;
}

void
PeleC::setup_ghost_particles(int ngrow)
{
    BL_PROFILE("PeleC::setup_ghost_particles()");
    BL_ASSERT(level < parent->finestLevel());
    if(PeleC::theSprayPC() != 0)
    {
        SprayParticleContainer::AoS ghosts;
        PeleC::theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
        PeleC::theGhostPC()->AddParticlesAtLevel(ghosts, level+1, ngrow);
    }
}

void
PeleC::remove_ghost_particles()
{
    BL_PROFILE("PeleC::remove_ghost_particles()");
    if (GhostPC != 0)
        GhostPC->RemoveParticlesAtLevel(level);
}


/**
 * Initialize the particles on the grid at level 0
 **/
void
PeleC::init_particles ()
{
    BL_PROFILE("PeleC::init_particles()");

    if (level > 0)
        return;

    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    amrex::ExecOnFinalize(RemoveParticlesOnExit);

    if (do_spray_particles)
    {
  	BL_ASSERT(theSprayPC() == 0);
	
	SprayPC = new SprayParticleContainer(parent,&phys_bc);
	theSprayPC()->SetVerbose(particle_verbose);

        if (parent->subCycle())
        {
            VirtPC = new SprayParticleContainer(parent,&phys_bc);
            GhostPC = new SprayParticleContainer(parent,&phys_bc);
        }
	
        // fuel properties from user input file
        import_fuel_properties(&particle_nfuel_species, partf_mf.data(), 
			       partf_crit_T.data(), partf_latent.data(), partf_boil_T.data(),
                               partf_cp.data(), partf_molwt.data(), partf_indx.data());

        import_control_parameters(&particle_heat_tran, &particle_mass_tran, &particle_mom_tran);

	if (! particle_init_file.empty())
	{
	    theSprayPC()->InitFromAsciiFile(particle_init_file,NSR_SPR);
	} 
        else if (particle_init_uniform > 0) 
        {
            // Initialize uniform particle distribution
            //  {vel(DIM), temperature, diameter, density}
            theSprayPC()->InitParticlesUniform(level, particle_init_uniform);
        }
    }
}


void
PeleC::particle_post_restart (const std::string& restart_file, bool is_checkpoint)
{
    if (level > 0)
       return;

    if (do_spray_particles)
    {
        BL_ASSERT(SprayPC == 0);

        SprayPC = new SprayParticleContainer(parent,&phys_bc);
        theSprayPC()->SetVerbose(particle_verbose);

        if (parent->subCycle())
        {
            VirtPC = new SprayParticleContainer(parent,&phys_bc);
            GhostPC = new SprayParticleContainer(parent,&phys_bc);
        }

        // fuel properties from user input file
        import_fuel_properties(&particle_nfuel_species, partf_mf.data(),
                               partf_crit_T.data(), partf_latent.data(), partf_boil_T.data(),
                               partf_cp.data(), partf_molwt.data(),partf_indx.data());

        import_control_parameters(&particle_heat_tran, &particle_mass_tran, &particle_mom_tran);

        // Make sure to call RemoveParticlesOnExit() on exit.
        //
        amrex::ExecOnFinalize(RemoveParticlesOnExit);

        theSprayPC()->Restart(parent->theRestartFile(), chk_particle_file, is_checkpoint);
    }
}

std::unique_ptr<MultiFab>
PeleC::particle_derive(const std::string& name,
		       Real               time,
		       int                ngrow)
{
    BL_PROFILE("PeleC::particle_derive()");

    if (theSprayPC() && name == "particle_count")
    {
        std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
	MultiFab    temp_dat(grids,dmap,1,0);
	temp_dat.setVal(0);
	theSprayPC()->Increment(temp_dat,level);
	MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
	return derive_dat;
    }
    else if (theSprayPC() && name == "total_particle_count")
    {
	//
	// We want the total particle count at this level or higher.
	//
	std::unique_ptr<MultiFab> derive_dat = particle_derive("particle_count",time,ngrow);

	IntVect trr(D_DECL(1,1,1));

	for (int lev = level+1; lev <= parent->finestLevel(); lev++)
	{
            auto ba = parent->boxArray(lev);
            const auto& dm = parent->DistributionMap(lev);
            MultiFab temp_dat(ba, dm, 1, 0);

            trr *= parent->refRatio(lev - 1);

            ba.coarsen(trr);
            MultiFab ctemp_dat(ba, dm, 1, 0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            theSprayPC()->Increment(temp_dat, lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab = temp_dat[mfi];
                FArrayBox& cfab = ctemp_dat[mfi];
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
                    const Real val = ffab(p);
                    if (val > 0)
                        cfab(amrex::coarsen(p, trr)) += val;
                }
            }

            temp_dat.clear();

            MultiFab dat(grids, dmap, 1, 0);
            dat.setVal(0);
            dat.copy(ctemp_dat);

            MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
	}

	return derive_dat;
    }
    else
    {
	return AmrLevel::derive(name,time,ngrow);
    }
}

void
PeleC::particle_redistribute (int lbase, bool init_part)
{
    BL_PROFILE("PeleC::particle_redistribute()");
    if (theSprayPC())
    {
        //  
        // If we are calling with init_part = true, then we want to force the redistribute
        //    without checking whether the grids have changed.
        //  
        if (init_part)
        {
            theSprayPC()->Redistribute(lbase);
            return;
        }

        //
        // These are usually the BoxArray and DMap from the last regridding.
        //
        static Vector<BoxArray>            ba;
        static Vector<DistributionMapping> dm;

        bool changed = false;

        int flev = parent->finestLevel();
	
        while ( parent->getAmrLevels()[flev] == nullptr ) {
            flev--;
	}
 
        if (ba.size() != flev+1)
        {
            ba.resize(flev+1);
            dm.resize(flev+1);
            changed = true;
        }
        else
        {
            for (int i = 0; i <= flev && !changed; i++)
            {
                if (ba[i] != parent->boxArray(i))
                    //
                    // The BoxArrays have changed in the regridding.
                    //
                    changed = true;

                if ( ! changed)
                {
                    if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
                        //
                        // The DistributionMaps have changed in the regridding.
                        //
                        changed = true;
                }
            }
        }

        if (changed)
        {
            //
            // We only need to call Redistribute if the BoxArrays or DistMaps have changed.
	    // We also only call it for particles >= lbase. This is
	    // because of we called redistribute during a subcycle, there may be particles not in
	    // the proper position on coarser levels.
            //
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "Calling redistribute because changed " << '\n';

            theSprayPC()->Redistribute(lbase);
            //
            // Use the new BoxArray and DistMap to define ba and dm for next time.
            //
            for (int i = 0; i <= flev; i++)
            {
                ba[i] = parent->boxArray(i);
                dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
            }
        }
        else
        {
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "NOT calling redistribute because NOT changed " << '\n';
        }
    }
}


void
PeleC::TimestampParticles (int ngrow)
{
#if 1
    return;
    amrex::Abort("Need to fill in PeleC::TimestampParticles");
#else
    static bool first = true;
    static int imax = -1;
    if (first)
    {
	first = false;

	ParmParse ppp("particles");

	// have to do it here, not in read_particle_params, because Density, ..., are set after
	// read_particle_params is called.

	int timestamp_density = 1;
	ppp.query("timestamp_density", timestamp_density);
	if (timestamp_density) {
	    timestamp_indices.push_back(Density);
	    std::cout << "Density = " << Density << std::endl;
	}

	int timestamp_temperature = 0;
	ppp.query("timestamp_temperature", timestamp_temperature);
	if (timestamp_temperature) {
	    timestamp_indices.push_back(Temp);
	    std::cout << "Temp = " << Temp << std::endl;
	}
	
	if (!timestamp_indices.empty()) {
	    imax = *(std::max_element(timestamp_indices.begin(), timestamp_indices.end()));
	}
    }

    if ( SprayPC && !timestamp_dir.empty())
    {
	std::string basename = timestamp_dir;
		
	if (basename[basename.length()-1] != '/') basename += '/';
	
	basename += "Timestamp";
	
	int finest_level = parent->finestLevel();
	Real time        = state[State_Type].curTime();

	for (int lev = level; lev <= finest_level; lev++)
	{
	    if (SprayPC->NumberOfParticlesAtLevel(lev) <= 0) continue;
	    
	    MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

	    if (imax >= 0) {  // FillPatchIterator will fail otherwise
		int ng = (lev == level) ? ngrow : 1;
		FillPatchIterator fpi(parent->getLevel(lev), S_new, 
				      ng, time, State_Type, 0, imax+1);
		const MultiFab& S = fpi.get_mf();
		SprayPC->Timestamp(basename, S    , lev, time, timestamp_indices);
	    } else {
		SprayPC->Timestamp(basename, S_new, lev, time, timestamp_indices);
	    }
	}
    }	
#endif
}

#endif
