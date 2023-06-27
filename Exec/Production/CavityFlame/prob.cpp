#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
    void
    amrex_probinit(
        const int* /*init*/,
        const int* /*name*/,
        const int* /*namelen*/,
        const amrex::Real* problo,
        const amrex::Real* probhi)
    {
        amrex::ParmParse pp("prob");

        //pp.query("T_ox", PeleC::h_prob_parm_device->T_ox);

    }

    void
    PeleC::problem_post_timestep()
    {
    }

    void
    PeleC::problem_post_init()
    {
    }

    void
    PeleC::problem_post_restart()
    {
    }
}
