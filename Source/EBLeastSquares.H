#include "EB.H"

AMREX_GPU_DEVICE
void get_qmat(amrex::Real delta_x_i[NLSPTS],amrex::Real delta_y_i[NLSPTS],
        amrex::Real delta_z_i[NLSPTS],int N,amrex::Real qmat[NEL_TRIMAT])
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
void get_weightvec(amrex::Real delta_x_i[NLSPTS],amrex::Real delta_y_i[NLSPTS],
        amrex::Real delta_z_i[NLSPTS],int N,amrex::Real qmat[NEL_TRIMAT],amrex::Real wvec[NLSPTS][AMREX_SPACEDIM])
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
        alpha_ij_3=(delta_z_i[i]-r23/r22*delta_y_i[i]+beta*delta_x_i[i])/(r33*r33);

        wvec[i][0]=alpha_ij_1 - r12/r11*alpha_ij_2 + beta*alpha_ij_3;
        wvec[i][1]=alpha_ij_2 - r23/r22*alpha_ij_3;
        wvec[i][2]=alpha_ij_3;
    }

}

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
  
  AMREX_ASSERT(AMREX_SPACEDIM > 1);

  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  const amrex::Real area = std::pow(dx, AMREX_SPACEDIM - 1);
  const amrex::Real fac = area / dx;

  //amrex::Print()<<"box lo:"<<lo.x<<"\t"<<lo.y<<"\t"<<lo.z<<"\n";
  //amrex::Print()<<"box hi:"<<hi.x<<"\t"<<hi.y<<"\t"<<hi.z<<"\n";

  amrex::ParallelFor(Nsten, [=] AMREX_GPU_DEVICE(int L) 
  {
    if (is_inside(ebg[L].iv, lo, hi)) 
    {
        const amrex::Real n[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_normal[0], ebg[L].eb_normal[1], ebg[L].eb_normal[2])};

        amrex::Real centcoord[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ebg[L].eb_centroid[0], ebg[L].eb_centroid[1],
        ebg[L].eb_centroid[2])};
      
        const int ivs[AMREX_SPACEDIM] = {
        AMREX_D_DECL(ebg[L].iv[0], ebg[L].iv[1], ebg[L].iv[2])};

        // From ivs, move to center of stencil, then move to lower-left of that
        const int baseiv[AMREX_SPACEDIM] = {AMREX_D_DECL(
        ivs[0] + (int)amrex::Math::copysign(1.0, n[0]) - 1,
        ivs[1] + (int)amrex::Math::copysign(1.0, n[1]) - 1,
        ivs[2] + (int)amrex::Math::copysign(1.0, n[2]) - 1)};

        //amrex::Print()<<"L:"<<L<<"\n";
        //amrex::Print()<<"ivs:"<<ivs[0]<<"\t"<<ivs[1]<<"\t"<<ivs[2]<<"\n";
        //amrex::Print()<<"normals:"<<n[0]<<"\t"<<n[1]<<"\t"<<n[2]<<"\n";
        //amrex::Print()<<"centcoord:"<<centcoord[0]<<"\t"<<centcoord[1]<<"\t"<<centcoord[2]<<"\n";
        //amrex::Print()<<"baseiv:"<<baseiv[0]<<"\t"<<baseiv[1]<<"\t"<<baseiv[2]<<"\n";
        //amrex::Print()<<"eb_area/vfrac:"<<ebg[L].eb_area<<"\t"<<ebg[L].eb_vfrac<<"\n";

        //set iv and iv_base
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) 
        {
          grad_stencil[L].iv[dir] = ebg[L].iv[dir];
          grad_stencil[L].iv_base[dir] = baseiv[dir];
        }

        //do not include self cell
        amrex::Real delta_x_i[NLSPTS],delta_y_i[NLSPTS],delta_z_i[NLSPTS];
        amrex::Real qmat[NEL_TRIMAT],wvec[NLSPTS][AMREX_SPACEDIM];
        amrex::Real xi[AMREX_SPACEDIM];
      
        int iter=0;

        for (int ii = 0; ii < 3; ii++) 
        {
           for (int jj = 0; jj < 3; jj++) 
           {
             for (int kk = 0; kk < 3; kk++) 
             {
                int global_i=baseiv[0]+ii;
                int global_j=baseiv[1]+jj;
                int global_k=baseiv[2]+kk;

                if(!(global_i==ivs[0] && global_j==ivs[1] && global_k==ivs[2]))
                        //!flags(global_i,global_j,global_k).isCovered())
                {
                    xi[0]=global_i-ivs[0];
                    xi[1]=global_j-ivs[1];
                    xi[2]=global_k-ivs[2];

                    delta_x_i[iter]  = (xi[0]-centcoord[0]);
                    delta_y_i[iter]  = (xi[1]-centcoord[1]);
                    delta_z_i[iter]  = (xi[2]-centcoord[2]);
                    iter++;
                }
            }
          }
        }
        //amrex::Print()<<"iter:"<<iter<<"\n";

        for(int npt=0;npt<NLSPTS;npt++)
        {
            //amrex::Print()<<delta_x_i[npt]<<"\t"<<delta_y_i[npt]<<"\t"<<delta_z_i[npt]<<"\n";
        }

        if(iter<3)
        {
                //amrex::Print()<<"iter:"<<iter<<"\n";
        }

        get_qmat(delta_x_i,delta_y_i,delta_z_i,iter,qmat);
        for(int npt=0;npt<NEL_TRIMAT;npt++)
        {
            //amrex::Print()<<qmat[npt]<<"\t";
        }
        //amrex::Print()<<"\n";
        get_weightvec(delta_x_i,delta_y_i,delta_z_i,iter,qmat,wvec);
        for(int npt=0;npt<NLSPTS;npt++)
        {
            //amrex::Print()<<wvec[npt][0]<<"\t"<<wvec[npt][1]<<"\t"<<wvec[npt][2]<<"\n";
        }


        iter=0;
        amrex::Real selfweight=0.0;
        for (int ii = 0; ii < 3; ii++) 
        {
            for (int jj = 0; jj < 3; jj++) 
            {
                for (int kk = 0; kk < 3; kk++) 
                {
                    int global_i=baseiv[0]+ii;
                    int global_j=baseiv[1]+jj;
                    int global_k=baseiv[2]+kk;
                
                    if(!(global_i==ivs[0] && global_j==ivs[1] && global_k==ivs[2]))
                            //!flags(global_i,global_j,global_k).isCovered())
                    {
                        grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) =
                        fac * ebg[L].eb_area * (wvec[iter][0]*n[0]
                                               +wvec[iter][1]*n[1]
                                               +wvec[iter][2]*n[2]);
                        selfweight += grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]);
                        iter++;
                    }
                    else
                    {
                        grad_stencil[L].val PELEC_D_TERM_REVERSE([kk], [jj], [ii]) = 0.0;
                    }
                }
            }
        }
        grad_stencil[L].bcval_sten = -selfweight;
    }
  });
}

