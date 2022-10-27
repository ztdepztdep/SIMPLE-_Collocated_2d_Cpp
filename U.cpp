/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2018  <copyright holder> <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "U.h"



U::U(const U & other):AGrid(other.AGrid),AllVar(other.AllVar)
{ }

U::~U()
{  }


void U::Setup()
{
    double Diff,Flow;
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

     for(int icol=1;icol<L1;icol++)
      {
        Flow=AGrid.Xcv[icol]*AllVar.MIV[0][icol]*AllVar.Rho[0][icol];
        Diff=AGrid.Xcv[icol]*AllVar.Mu[0][icol]/AGrid.Ydif[1];
        AllVar.ArowM[1][icol]=Diff+max(0.0,Flow);
      }


     for(int jrow=Jstart;jrow<Jend;jrow++)
      {
        Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][0]*AllVar.Rho[jrow][0];
        Diff=AGrid.Ycv[jrow]*AllVar.Mu[jrow][0]/AGrid.Xdif[1] ;
        AllVar.AcolM[jrow][1]=Diff+max(0.0,Flow);  // the most left column

       for(int icol=1;icol<L1;icol++)   ///the first inner and the last inner node
        {
            if(icol==L1-1)
            {                                        //icol=L1 values,
                Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][icol]*AllVar.Rho[jrow][L1];
                Diff=AGrid.Ycv[jrow]*AllVar.Mu[jrow][L1]/AGrid.Xdif[L1] ;
            }
            else
            {
                Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][icol];
                Flow=Flow*(AGrid.Fx[icol+1]*AllVar.Rho[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Rho[jrow][icol]);
                Diff=2.0*AGrid.Ycv[jrow]*AllVar.Mu[jrow][icol]*AllVar.Mu[jrow][icol+1];
                Diff=Diff/(AGrid.Xcv[icol+1]*AllVar.Mu[jrow][icol]+AGrid.Xcv[icol]*AllVar.Mu[jrow][icol+1]);
            }
            AllVar.AcolM[jrow][icol+1]=Diff+max(0.0,Flow);
            AllVar.AcolP[jrow][icol]=AllVar.AcolM[jrow][icol+1]-Flow;

         //   if(icol==1) cout <<" "<<setw(13)<<fixed<<setprecision(8)<<AllVar.AcolP[jrow][icol]<<"  "<<Flow<<endl;


            //    cout<<jrow<<"  "<<icol<<"  "<<AllVar.AcolP[jrow][icol]<<endl;


        if(jrow==M1-1)
        {
             Flow=AGrid.Xcv[icol]*AllVar.MIV[jrow][icol]*AllVar.Rho[M1][icol];
              Diff=AGrid.Xcv[icol]*AllVar.Mu[M1][icol]/AGrid.Ydif[M1];
        }
        else{
                Flow=AGrid.Xcv[icol]*AllVar.MIV[jrow][icol];
                Flow=Flow*(AGrid.Fy[jrow+1]*AllVar.Rho[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Rho[jrow][icol] );
                Diff=AGrid.Xcv[icol]*2.0*AllVar.Mu[jrow][icol]*AllVar.Mu[jrow+1][icol];
                Diff=Diff/(AGrid.Ycv[jrow+1]*AllVar.Mu[jrow][icol]+AGrid.Ycv[jrow]*AllVar.Mu[jrow+1][icol]);
           }

              AllVar.ArowM[jrow+1][icol]=Diff+max(0.0,Flow);
              AllVar.ArowP[jrow][icol]=AllVar.ArowM[jrow+1][icol]-Flow;

               double VOL=AGrid.Ycv[jrow]*AGrid.Xcv[icol];

               double APT=AllVar.Rho[jrow][icol] /AllVar.Dt;

               AllVar.Ap[jrow][icol]= 0.0-APT;

               double SC=0.0;

               AllVar.Con[jrow][icol]=SC+APT*AllVar.Uold1[jrow][icol]-AllVar.Rho[jrow][icol]*AllVar.Gx*AllVar.Beta*(AllVar.T[jrow][icol]-AllVar.Treference);
                AllVar.Ap[jrow][icol]=-AllVar.Ap[jrow][icol]*VOL+AllVar.AcolM[jrow][icol]+AllVar.AcolP[jrow][icol]
                                           +AllVar.ArowM[jrow][icol]+AllVar.ArowP[jrow][icol];

                AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]*VOL;
                AllVar.UCon[jrow][icol]=AllVar.Con[jrow][icol];


                AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]+AllVar.Ap[jrow][icol]*AllVar.U[jrow][icol]/AllVar.EUV;
                AllVar.Ap[jrow][icol]=AllVar.Ap[jrow][icol]*( 1.0+1.0/AllVar.EUV );

                AllVar.UArowM[jrow][icol]=AllVar.ArowM[jrow][icol];
                AllVar.UArowP[jrow][icol]=AllVar.ArowP[jrow][icol];
                AllVar.UAcolM[jrow][icol]=AllVar.AcolM[jrow][icol];
                AllVar.UAcolP[jrow][icol]=AllVar.AcolP[jrow][icol];
                AllVar.UAp[jrow][icol]=AllVar.Ap[jrow][icol];

                double PE=AGrid.Fx[icol+1]*AllVar.P[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.P[jrow][icol];
                double PW=AGrid.Fx[icol]*AllVar.P[jrow][icol]+AGrid.Fxm[icol]*AllVar.P[jrow][icol-1];
                AllVar.Con[jrow][icol]=AllVar.Con[jrow][icol]-(PE-PW)*AGrid.Ycv[jrow];

            //cout<<jrow<<"  "<<icol<<"  "<<AllVar.Ap[jrow][icol]<<"  "<<AllVar.Con[jrow][icol]<<endl;
        }
    }
    AllVar.SUD(AllVar.U);
  //   AllVar.Show(AllVar.UAp);
}




void U::GenerateUhat()
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        for(int icol=1;icol<L1;icol++)
        {

            AllVar.Uhat[jrow][icol]=AllVar.UAcolM[jrow][icol]*AllVar.U[jrow][icol-1]
            +AllVar.UAcolP[jrow][icol]*AllVar.U[jrow][icol+1]
            +AllVar.UArowM[jrow][icol]*AllVar.U[jrow-1][icol]
            +AllVar.UArowP[jrow][icol]*AllVar.U[jrow+1][icol]
            +AllVar.UCon[jrow][icol];

            AllVar.Uhat[jrow][icol]=AllVar.Uhat[jrow][icol]/AllVar.UAp[jrow][icol];
            AllVar.Du[jrow][icol]=AGrid.Ycv[jrow]/AllVar.UAp[jrow][icol];

            AllVar.Tu[jrow][icol]=-0.5*AllVar.Rho[jrow][icol]*AllVar.Gx*AllVar.Beta/AllVar.UAp[jrow][icol];
        }

    }


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        for(int icol=1;icol<L1-1;icol++)
        {
            double Uhate=AGrid.Fx[icol+1]*AllVar.Uhat[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Uhat[jrow][icol];
            double Due=AGrid.Fx[icol+1]*AllVar.Du[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Du[jrow][icol];
            AllVar.MIU[jrow][icol]=Uhate+Due*(AllVar.P[jrow][icol]-AllVar.P[jrow][icol+1]);
            AllVar.MIU[jrow][icol]=AllVar.MIU[jrow][icol]*( 1.0+1.0/AllVar.EUV );
        }

    }

}
