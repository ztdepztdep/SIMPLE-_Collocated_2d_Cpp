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

#include "Temperature.h"

void Temperature::Setup()
{

    double Diff,Flow;
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int icol=1;icol<L1;icol++)
    {
        Flow=AGrid.Xcv[icol]*AllVar.MIV[0][icol]*AllVar.RhoCp[0][icol];
        Diff=AGrid.Xcv[icol]*AllVar.Lamda[0][icol]/AGrid.Ydif[1];
        AllVar.ArowM[1][icol]=Diff+max(0.0,Flow);
    }


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][0]*AllVar.RhoCp[jrow][0];
        Diff=AGrid.Ycv[jrow]*AllVar.Lamda[jrow][0]/AGrid.Xdif[1] ;
        AllVar.AcolM[jrow][1]=Diff+max(0.0,Flow);  // the most left column

        for(int icol=1;icol<L1;icol++)   ///the first inner and the last inner node
        {
            if(icol==L1-1)
            {                                        //icol=L1 values,
                Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][icol]*AllVar.RhoCp[jrow][L1];
                Diff=AGrid.Ycv[jrow]*AllVar.Lamda[jrow][L1]/AGrid.Xdif[L1] ;
            }
            else
            {
                Flow=AGrid.Ycv[jrow]*AllVar.MIU[jrow][icol];
                Flow=Flow*(AGrid.Fx[icol+1]*AllVar.RhoCp[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.RhoCp[jrow][icol]);
                Diff=2.0*AGrid.Ycv[jrow]*AllVar.Lamda[jrow][icol]*AllVar.Lamda[jrow][icol+1];
                Diff=Diff/(AGrid.Xcv[icol+1]*AllVar.Lamda[jrow][icol]+AGrid.Xcv[icol]*AllVar.Lamda[jrow][icol+1]);
            }
            AllVar.AcolM[jrow][icol+1]=Diff+max(0.0,Flow);
            AllVar.AcolP[jrow][icol]=AllVar.AcolM[jrow][icol+1]-Flow;

            //   if(icol==1) cout <<" "<<setw(13)<<fixed<<setprecision(8)<<AllVar.AcolP[jrow][icol]<<"  "<<Flow<<endl;


            //    cout<<jrow<<"  "<<icol<<"  "<<AllVar.AcolP[jrow][icol]<<endl;


            if(jrow==M1-1)
            {
                Flow=AGrid.Xcv[icol]*AllVar.MIV[jrow][icol]*AllVar.RhoCp[M1][icol];
                Diff=AGrid.Xcv[icol]*AllVar.Lamda[M1][icol]/AGrid.Ydif[M1];
            }
            else{
                Flow=AGrid.Xcv[icol]*AllVar.MIV[jrow][icol];
                Flow=Flow*(AGrid.Fy[jrow+1]*AllVar.RhoCp[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.RhoCp[jrow][icol] );
                Diff=AGrid.Xcv[icol]*2.0*AllVar.Lamda[jrow][icol]*AllVar.Lamda[jrow+1][icol];
                Diff=Diff/(AGrid.Ycv[jrow+1]*AllVar.Lamda[jrow][icol]+AGrid.Ycv[jrow]*AllVar.Lamda[jrow+1][icol]);
            }

            AllVar.ArowM[jrow+1][icol]=Diff+max(0.0,Flow);
            AllVar.ArowP[jrow][icol]=AllVar.ArowM[jrow+1][icol]-Flow;

            double VOL=AGrid.Ycv[jrow]*AGrid.Xcv[icol];

            double APT=AllVar.RhoCp[jrow][icol] /AllVar.Dt;

            AllVar.Ap[jrow][icol]= 0.0-APT;

            double SC=0.0;

            AllVar.Con[jrow][icol]=SC+APT*AllVar.Told1[jrow][icol];
            AllVar.Ap[jrow][icol]=-AllVar.Ap[jrow][icol]*VOL+AllVar.AcolM[jrow][icol]+AllVar.AcolP[jrow][icol]
                                      +AllVar.ArowM[jrow][icol]+AllVar.ArowP[jrow][icol];

            AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]*VOL;


            AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]+AllVar.Ap[jrow][icol]*AllVar.T[jrow][icol]/AllVar.EOTHER;
            AllVar.Ap[jrow][icol]=AllVar.Ap[jrow][icol]*( 1.0+1.0/AllVar.EOTHER );


     //cout<<jrow<<"  "<<icol<<"  "<<AllVar.Ap[jrow][icol]<<"  "<<AllVar.Con[jrow][icol]<<endl;
        }
    }
    AllVar.SUD(AllVar.T);
    //   AllVar.Show(AllVar.UAp);
}
