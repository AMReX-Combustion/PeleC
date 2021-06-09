#include "EB.H"

void
pc_fill_bndry_grad_stencil_ls(
  const amrex::Box& bx,
  const amrex::Real dx,
  const int /*Nebg*/,
  const EBBndryGeom* ebg,
  const int Nsten,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  EBBndrySten* grad_stencil)
{
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  const amrex::Real area = std::pow(dx, AMREX_SPACEDIM - 1);
  const amrex::Real fac = area / dx;

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) 
  {
    if (is_inside(ebg[L].iv, lo, hi)) 
    {
        const amrex::Real n[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};

        // From ivs, move to center of stencil, then move to lower-left of that
        const int baseiv[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ivs[0] + (int)amrex::Math::copysign(1.0, n[0]) - 1,
        ivs[1] + (int)amrex::Math::copysign(1.0, n[1]) - 1,
        ivs[2] + (int)amrex::Math::copysign(1.0, n[2]) - 1)};

        amrex::Real delta_x_i[27],delta_y_i[27],delta_z_i[27];
        amrex::Real qmat[6],wvec[27][3];

        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) 
        {
          grad_stencil[L].iv[dir] = ebg[L].iv[dir];
          grad_stencil[L].iv_base[dir] = baseiv[dir];
        }

        get_qmat(delta_x_i,delta_y_i,delta_z_i,N,qmat);
        get_weightvec(delta_x_i,delta_y_i,delta_z_i,N,qmat,wvec);


      for (int ii = 0; ii < 3; ii++) 
      {
        for (int jj = 0; jj < 3; jj++) 
        {
#if AMREX_SPACEDIM > 2
          for (int kk = 0; kk < 3; kk++) 
          {
#endif
            grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
              fac * ebg[L].eb_area * tsten AMREX_D_TERM([ii], [jj], [kk]);
#if AMREX_SPACEDIM > 2
          }
#endif
        }
      }
      grad_stencil[L].bcval_sten = fac * ebg[L].eb_area * bcs;
    }
  });
}

AMREX_GPU_DEVICE
void get_qmat(amrex::Real delta_x_i[27],amrex::Real delta_y_i[27],
        amrex::Real delta_z_i[27],int N,amrex::Real qmat[6])
{
        amrex::Real r11=0.0;
        amrex::Real r12=0.0;
        amrex::Real r22=0.0;
        amrex::Real r13=0.0;
        amrex::Real r23=0.0;
        amrex::Real r33=0.0;

        for(int i=0;i<N;i++)
        {
                r11 += pow(delta_x_i[i],2.0);
                r12 += delta_x_i[i]*delta_y_i[i];
                r22 += pow(delta_y_i[i],2.0);
                r13 += delta_x_i[i]*delta_z_i[i];
                r23 += delta_y_i[i]*delta_z_i[i];
                r33 += pow(delta_z_i[i],2.0);
        }

        r11=sqrt(r11);
        r12=r12/r11;
        r22=sqrt(r22-r12*r12);
        r13=r13/r11;
        r23=r23/r22-r12*r13/r22;
        r33=sqrt(r33-r13*r13-r23*r23);

        qmat[0]=r11;
        qmat[1]=r12;
        qmat[2]=r22;
        qmat[3]=r13;
        qmat[4]=r23;
        qmat[5]=r33;
}

AMREX_GPU_DEVICE
void get_weightvec(amrex::Real delta_x_i[27],amrex::Real delta_y_i[27],
        amrex::Real delta_z_i[27],int N,amrex::Real qmat[6],amrex::Real wvec[27,3])
{
    amrex::Real r11=qmat[0];
    amrex::Real r12=qmat[1];
    amrex::Real r22=qmat[2];
    amrex::Real r13=qmat[3];
    amrex::Real r23=qmat[4];
    amrex::Real r33=qmat[5];

    amrex::Real alpha_ij_1,alpha_ij_2,alpha_ij_3,beta;

    beta=(r12*r23-r13*r22)/(r11*r22);
 
    for(int i=0;i<N;i++)
   {
       alpha_ij_1=delta_x_i[i]/(r11*r11);
       alpha_ij_2=(delta_y_i[i]-r12/r11*delta_x_i[i])/(r22*r22);
       alpha_ij_3=(delta_z_i[i]-r23/r22*delta_y_i[i]+beta*delta_x_i[i]);

       wvec[i][0]=alpha_ij_1 - r12/r11*alpha_ij_2+beta*alpha_ij_3;
       wvec[i][1]=alpha_ij_2 - r23/r22*alpha_ij_3;
       wvec[i][2]=alpha_ij_3;
   }
        
}
