#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{

    amrex::Real pamb      = 101325.0;
    amrex::Real v_in      = 18.407;
    amrex::Real T_in      = 300.0;
    amrex::Real inflow_diam = 6.604e-3;
    amrex::Real thermocouple_datum = 5.31e-3;
    amrex::Real z_wall1    = 25e-3;
    amrex::Real T_wall1    = 300.0;
    amrex::Real z_wall2    = 50e-3;
    amrex::Real T_wall2    = 300.0; 
    amrex::Real z_wall3    = 75e-3;
    amrex::Real T_wall3    = 300.0;   
    amrex::Real z_wall4    = 125e-3;
    amrex::Real T_wall4    = 300.0;
    amrex::Real z_wall5    = 150e-3;
    amrex::Real T_wall5    = 300.0;

};

struct ProbParmHost
{


};

#endif
