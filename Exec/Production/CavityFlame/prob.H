#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"

#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Utilities.H"
#include "turbinflow.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real Cs;
  massfrac[O2_ID]=prob_parm.Yox_amb;
  massfrac[N2_ID]=1.0-prob_parm.Yox_amb;


  amrex::Real rho = 0.0, eint = 0;
  eos.PYT2R(prob_parm.Pres_amb, massfrac, prob_parm.Temp_amb, rho);
  eos.RTY2E(rho, prob_parm.Temp_amb, massfrac, eint);
  eos.RTY2Cs(rho, prob_parm.Temp_amb, massfrac, Cs);
  
  amrex::Real x_velocity = Cs*prob_parm.Mach_channel;
  amrex::Real y_velocity = 0.0;
  amrex::Real z_velocity = 0.0;

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * x_velocity;
  state(i, j, k, UMY) = rho * y_velocity;
  state(i, j, k, UMZ) = rho * z_velocity;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) =
    rho * (eint + 0.5 * (x_velocity * x_velocity + y_velocity * y_velocity +
                         z_velocity * z_velocity));
  state(i, j, k, UTEMP) = prob_parm.Temp_amb;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::Real s_int[NVAR],
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real time,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
{

    if (idir == 1) 
    {

        // default wall condition
        s_ext[URHO] = s_int[URHO];
        s_ext[UMX] = -s_int[UMX];
        s_ext[UMY] = -s_int[UMY];
        s_ext[UMZ] = -s_int[UMZ];
        s_ext[UEINT] = s_int[UEINT];
        s_ext[UEDEN] = s_int[UEDEN];
        s_ext[UTEMP] = s_int[UTEMP];
        for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[UFS + n] = s_int[UFS + n];
        }

        for(int nh=0;nh<prob_parm.nholes;nh++)
        {
            amrex::Real dphi_angle=2.0*M_PI/amrex::Real(prob_parm.nholes);
            amrex::Real phi_angle = nh*dphi_angle;
            amrex::Real hole_cx = prob_parm.centx + prob_parm.r_circ*cos(phi_angle);
            amrex::Real hole_cz = prob_parm.centz + prob_parm.r_circ*sin(phi_angle);

            const amrex::Real radius = std::sqrt(
                (x[0] - hole_cx) * (x[0] - hole_cx) +
                (x[2] - hole_cz) * (x[2] - hole_cz));

            if (radius < prob_parm.r_hole) 
            {
                amrex::Real massfrac[NUM_SPECIES] = {0.0};
                auto eos = pele::physics::PhysicsType::eos();

                // interior pressure
                amrex::Real p_int;
                amrex::Real massfrac_int[NUM_SPECIES] = {0.0};
                amrex::Real massfrac_ext[NUM_SPECIES] = {0.0};
                for (int n = 0; n < NUM_SPECIES; n++) {
                    massfrac_int[n] = s_int[UFS + n] / s_int[URHO];
                }
                eos.RTY2P(s_int[URHO], s_int[UTEMP], massfrac_int, p_int);

                amrex::Real decayfac=0.0;
                amrex::Real eta = (1.0 - std::pow(std::tanh(decayfac*radius/prob_parm.r_hole),2.0));
                massfrac_ext[H2_ID]=prob_parm.Yfuel_jet*eta;
                massfrac_ext[O2_ID]=0.0;
                massfrac_ext[N2_ID]=1.0-massfrac_ext[H2_ID]-massfrac_ext[O2_ID];
                amrex::Real vjet_pipeflow = prob_parm.vel_jet*(1.0 - std::pow(radius/prob_parm.r_hole,2.0));
                amrex::Real vy_in = vjet_pipeflow*cos(prob_parm.cone_angle*M_PI/180.0);
                amrex::Real vr_in = vjet_pipeflow*sin(prob_parm.cone_angle*M_PI/180.0);
                amrex::Real vx_in = vr_in*cos(phi_angle);
                amrex::Real vz_in = vr_in*sin(phi_angle);

                s_ext[UTEMP] = prob_parm.temp_jet; 
                amrex::Real rho, eint;
                eos.PYT2RE(p_int, massfrac, s_ext[UTEMP],rho,eint);

                s_ext[URHO] = rho;
                s_ext[UEINT] = rho * eint;

                for (int n = 0; n < NUM_SPECIES; n++) {
                    s_ext[UFS + n] = rho * massfrac[n];
                }
                // Velocity
                s_ext[UMX] = rho * vx_in;
                s_ext[UMY] = rho * vy_in;
                s_ext[UMZ] = rho * vz_in;
                s_ext[UEDEN] = s_ext[UEINT] + 0.5 * rho *
                (vx_in * vx_in + vy_in * vy_in + vz_in * vz_in);
            }
        }
    }
    if(idir == 0)
    {

        //supersonic inflow        
        amrex::Real rho,eint,Cs;
        s_ext[UTEMP]=prob_parm.Temp_amb;
        auto eos = pele::physics::PhysicsType::eos();
        amrex::Real massfrac_ext[NUM_SPECIES] = {0.0};
        massfrac_ext[O2_ID]=prob_parm.Yox_channel;
        massfrac_ext[H2_ID]=0.0;
        massfrac_ext[N2_ID]=1-prob_parm.Yox_channel-massfrac_ext[H2_ID];
        eos.PYT2RE(prob_parm.Pres_amb, massfrac_ext, s_ext[UTEMP],rho,eint);
        eos.RTY2Cs(rho, prob_parm.Temp_amb, massfrac_ext, Cs);

        s_ext[URHO]   = rho;
        s_ext[UMX]    = rho * Cs * prob_parm.Mach_channel;
        s_ext[UMY]    = 0.0;
        s_ext[UMZ]    = 0.0;
        s_ext[UEINT]  = rho*eint;
        s_ext[UEDEN]  = s_ext[UEINT] + 0.5*s_ext[UMX]*s_ext[UMX]/rho;

        for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[UFS + n] = rho * massfrac_ext[n];
        }
    }
}

struct MyProbTagStruct
{
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    static void set_problem_tags(
        const int /*i*/,
        const int /*j*/,
        const int /*k*/,
        amrex::Array4<char> const& /*tag*/,
        amrex::Array4<amrex::Real const> const& /*field*/,
        char /*tagval*/,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
        const amrex::Real /*time*/,
        const int /*level*/,
        ProbParmDevice const& /*d_prob_parm_device*/) noexcept
    {
        // could do problem specific tagging here
    }
};

using ProblemTags = MyProbTagStruct;

struct MyProbDeriveStruct
{
    static void
    add(amrex::DeriveList& /*derive_lst*/, amrex::DescriptorList& /*desc_lst*/)
    {
        // Add derives as follows and define the derive function below:
        // derive_lst.add(
        //  "varname", amrex::IndexType::TheCellType(), 1, pc_varname,
        //  the_same_box);
        // derive_lst.addComponent("varname", desc_lst, State_Type, 0, NVAR);
    }

    static void pc_varname(
        const amrex::Box& /*bx*/,
        amrex::FArrayBox& /*derfab*/,
        int /*dcomp*/,
        int /*ncomp*/,
        const amrex::FArrayBox& /*datfab*/,
        const amrex::Geometry& /*geomdata*/,
        amrex::Real /*time*/,
        const int* /*bcrec*/,
        int /*level*/)
    {
        // auto const dat = datfab.array();
        // auto arr = derfab.array();
        // amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        // { do something with arr
        // });
    }
};

void pc_prob_close();

using ProblemDerives = MyProbDeriveStruct;

#endif
