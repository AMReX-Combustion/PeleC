
#include "AMReX_LevelBld.H"
#include "PeleC.H"

using namespace amrex;

class PeleCBld
    :
    public LevelBld
{
    virtual void variableSetUp ();
    virtual void variableCleanUp ();
    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                       const DistributionMapping& dm,
                                  Real            time);
};

PeleCBld PeleC_bld;

LevelBld*
getLevelBld ()
{
    return &PeleC_bld;
}

void
PeleCBld::variableSetUp ()
{
    PeleC::variableSetUp();
}

void
PeleCBld::variableCleanUp ()
{
    PeleC::variableCleanUp();
}

AmrLevel*
PeleCBld::operator() ()
{
    return new PeleC;
}

AmrLevel*
PeleCBld::operator() (Amr&            papa,
                       int             lev,
                       const Geometry& level_geom,
                       const BoxArray& ba,
                       const DistributionMapping& dm,
                       Real            time)
{
    return new PeleC(papa, lev, level_geom, ba, dm, time);
}
