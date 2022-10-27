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

#include "Solver.h"
#include <time.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Solver::SimpleAlgorithm(const int Nstep)
{
      //  system ("rm ./output/*.*");
        cout<<" The SIMPLE Solver is initlized ........"<<endl;

        for (int n = 1; n <=Nstep; n++)
        {
            AllVar.Uold1=AllVar.U;
            AllVar.Vold1=AllVar.V;

            AllVar.Rhoold1=AllVar.Rho;

            AllVar.UResidual=0.0;
            AllVar.VResidual=0.0;
            AllVar.MResidual=0.0;
            AllVar.TResidual=0.0;
        
            int counter=0;
            bool condition=true;

            while(condition)
            {
               UEq.Setup();
               AllVar.EvaluateResidual(AllVar.U,AllVar.UResidual);
               MSIP(AllVar.U,2);
               // Jacobi(AllVar.U,5);

                 VEq.Setup();
                AllVar.EvaluateResidual(AllVar.V,AllVar.VResidual);
                MSIP(AllVar.V,2);
              //  Jacobi(AllVar.V,5);

                UEq.GenerateUhat();
                VEq.GenerateVhat();

                PressureEq.Setup();
                MSIP(AllVar.Pc,AllVar.SolveNtime);
               // Jacobi(AllVar.Pc,5);

                PressureEq.UpdatePressure();
                PressureEq.UpdateMassFlux();

                TEq.Setup();
                AllVar.EvaluateResidual(AllVar.T,AllVar.TResidual);
                SIP(AllVar.T,AllVar.SolveNtime);


                AllVar.UResVec.push_back(AllVar.UResidual);
                AllVar.VResVec.push_back(AllVar.VResidual);
                AllVar.MResVec.push_back(AllVar.MResidual);
                AllVar.TResVec.push_back(AllVar.TResidual);
                cout << setprecision(4)<<scientific<<" U :  " << AllVar.UResidual <<"  V l: "<<AllVar.VResidual
                <<"  Mass: "<<AllVar.MResidual<<"  T : "<<AllVar.TResidual<<"  Counter "<< counter<<" "<<n<<"  with time "<<"  at time "<<AllVar.Time<<" Dt "<<AllVar.Dt<<endl;

                 if(counter>=AllVar.BreakStep) break;
                    counter++;

                if(counter>0)  condition=(AllVar.UResidual>AllVar.UVtol || AllVar.VResidual>AllVar.UVtol || AllVar.MResidual>AllVar.Mtol || AllVar.TResidual>1E-6);
            }
            AllVar.Time=AllVar.Time+AllVar.Dt;          
        }
     
      
        string FileName ="ResSIMPLE.dat";
        FileName="./output/"+FileName;
        ofstream outfile;
        outfile.open(FileName);

       for(int i=3;i<AllVar.UResVec.size();i++)
                 outfile<< setprecision(15)<<scientific<<i<<" "<<AllVar.UResVec[i]<<"  "<<AllVar.VResVec[i]<<"  "<<AllVar.MResVec[i]<<"  "<<AllVar.TResVec[i]<<endl;
       outfile.close();
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Solver::Jacobi(vector<vector<double>>& output,const int& Ntimes)
{
    int IST=1;
    int JST=1;
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    for(int k=0;k<Ntimes;k++)
    {
    for(int icol=1;icol<L1;icol++)
        for(int jrow=1; jrow<M1-1;jrow++)
       {  output[jrow][icol]=AllVar.Con[jrow][icol]+AllVar.ArowM[jrow][icol]*output[jrow-1][icol]
            +AllVar.ArowP[jrow][icol]*output[jrow+1][icol]+AllVar.AcolM[jrow][icol]*output[jrow][icol-1]
            +AllVar.AcolP[jrow][icol]*output[jrow][icol+1];
            output[jrow][icol]=output[jrow][icol]/AllVar.Ap[jrow][icol];

    }
       }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Solver::SIP(vector<vector<double>>& output,const int& Ntimes)
{
    double P1,P2,ALFA=0.9;
    int IST=1;
    int JST=1;
    int L1=AGrid.L1;
    int M1=AGrid.M1;

    vector<vector<double>> LW;
    vector<vector<double>> LS;
    vector<vector<double>> LPR;
    vector<vector<double>> UUN;
    vector<vector<double>> UUE;
    vector<vector<double>> RES;

    LW.resize(AGrid.LocalRowSize);
    for(auto&&i:LW) i.assign(AGrid.L1+1,0.0);

    LS.resize(AGrid.LocalRowSize);
    for(auto&&i:LS) i.assign(AGrid.L1+1,0.0);

    LPR.resize(AGrid.LocalRowSize);
    for(auto&&i:LPR) i.assign(AGrid.L1+1,0.0);

    UUN.resize(AGrid.LocalRowSize);
    for(auto&&i:UUN) i.assign(AGrid.L1+1,0.0);

    UUE.resize(AGrid.LocalRowSize);
    for(auto&&i:UUE) i.assign(AGrid.L1+1,0.0);

    RES.resize(AGrid.LocalRowSize);
    for(auto&&i:RES) i.assign(AGrid.L1+1,0.0);

    for(int icol=1;icol<AGrid.L1;icol++) {
        for(int jrow=1;jrow<AGrid.LocalRowSize-1;jrow++) {
            LW[jrow][icol]=-AllVar.AcolM[jrow][icol]/(1.0+ALFA*UUN[jrow][icol-1]);
            LS[jrow][icol]=-AllVar.ArowM[jrow][icol]/(1.0+ALFA*UUE[jrow-1][icol]);
            P1=ALFA*LW[jrow][icol]*UUN[jrow][icol-1];
            P2=ALFA*LS[jrow][icol]*UUE[jrow-1][icol];
            LPR[jrow][icol]=1.0/(AllVar.Ap[jrow][icol]+P1+P2-LW[jrow][icol]*UUE[jrow][icol-1]-LS[jrow][icol]*UUN[jrow-1][icol]+1.E-30);
            UUN[jrow][icol]=(-AllVar.ArowP[jrow][icol]-P1)*LPR[jrow][icol];
            UUE[jrow][icol]=(-AllVar.AcolP[jrow][icol]-P2)*LPR[jrow][icol];
        }
    }


    for(int k=0;k<Ntimes;k++)
    {
        for(int icol=1;icol<AGrid.L1;icol++)
            for(int jrow=1; jrow<AGrid.LocalRowSize-1;jrow++)
            {
                RES[jrow][icol]=AllVar.Con[jrow][icol]-AllVar.Ap[jrow][icol]*output[jrow][icol]+AllVar.ArowM[jrow][icol]*output[jrow-1][icol]
                +AllVar.ArowP[jrow][icol]*output[jrow+1][icol]+AllVar.AcolM[jrow][icol]*output[jrow][icol-1]
                +AllVar.AcolP[jrow][icol]*output[jrow][icol+1];
                RES[jrow][icol]=(RES[jrow][icol]-LS[jrow][icol]*RES[jrow-1][icol]-LW[jrow][icol]*RES[jrow][icol-1])*LPR[jrow][icol];
            }


            for(int icol=AGrid.L1-1;icol>0;icol--)
                for(int jrow=AGrid.LocalRowSize-2; jrow>0;jrow--)
                {
                    RES[jrow][icol]=RES[jrow][icol]-UUN[jrow][icol]*RES[jrow+1][icol]-UUE[jrow][icol]*RES[jrow][icol+1];
                    output[jrow][icol]=output[jrow][icol]+RES[jrow][icol];
                }
    }

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Solver::MSIP(vector<vector<double>>& output,const int& Ntimes)
{

    double P1,P2,ALFA=0.85;
    int IST=1;
    int JST=1;
    int L1=AGrid.L1;
    int M1=AGrid.M1;

    vector<vector<double>> LW;
    vector<vector<double>> LS;
    vector<vector<double>> LPR;
    vector<vector<double>> UUN;
    vector<vector<double>> UUE;
    vector<vector<double>> RES;
    vector<vector<double>> C;
    vector<vector<double>> G;

    LW.resize(AGrid.LocalRowSize);
    for(auto&&i:LW) i.assign(AGrid.L1+1,0.0);

    LS.resize(AGrid.LocalRowSize);
    for(auto&&i:LS) i.assign(AGrid.L1+1,0.0);

    LPR.resize(AGrid.LocalRowSize);
    for(auto&&i:LPR) i.assign(AGrid.L1+1,0.0);

    UUN.resize(AGrid.LocalRowSize);
    for(auto&&i:UUN) i.assign(AGrid.L1+1,0.0);

    UUE.resize(AGrid.LocalRowSize);
    for(auto&&i:UUE) i.assign(AGrid.L1+1,0.0);

    RES.resize(AGrid.LocalRowSize);
    for(auto&&i:RES) i.assign(AGrid.L1+1,0.0);

    C.resize(AGrid.LocalRowSize);
    for(auto&&i:C) i.assign(AGrid.L1+1,0.0);

    G.resize(AGrid.LocalRowSize);
    for(auto&&i:G) i.assign(AGrid.L1+1,0.0);


    for(int icol=1;icol<AGrid.L1;icol++) {
        for(int jrow=1;jrow<AGrid.LocalRowSize-1;jrow++) {
            double fim=UUN[jrow][icol-1];
            double fim1=UUN[jrow+1][icol-1];
            LW[jrow][icol]=-AllVar.AcolM[jrow][icol]/(1.0-ALFA*fim*fim1);

            double bi=LW[jrow][icol];
            C[jrow][icol]=-bi*fim;

            LS[jrow][icol]=(-AllVar.ArowM[jrow][icol]-bi*G[jrow][icol-1])/(1.0+2.0*ALFA*G[jrow-1][icol]);
            double di=LS[jrow][icol];

            P1=2.0*ALFA*fim1-G[jrow+1][icol-1];
            P2=2.0*ALFA*G[jrow-1][icol]-UUN[jrow-1][icol];
            LPR[jrow][icol]=1.0/(AllVar.Ap[jrow][icol]-bi*UUE[jrow][icol-1]+C[jrow][icol]*P1+P2*di+1.E-30);

            double F= UUE[jrow+1][icol-1]+2.0*ALFA*fim1;
            UUN[jrow][icol]=(-AllVar.ArowP[jrow][icol]-C[jrow][icol]*F)*LPR[jrow][icol];

            G[jrow][icol]=-di*UUE[jrow-1][icol]*LPR[jrow][icol];

            UUE[jrow][icol]=(-AllVar.AcolP[jrow][icol]-ALFA*di*G[jrow-1][icol])*LPR[jrow][icol];
        }
    }


    for(int k=0;k<Ntimes;k++)
    {
        for(int icol=1;icol<AGrid.L1;icol++)
            for(int jrow=1; jrow<AGrid.LocalRowSize-1;jrow++)
            {
                RES[jrow][icol]=AllVar.Con[jrow][icol]-AllVar.Ap[jrow][icol]*output[jrow][icol]+AllVar.ArowM[jrow][icol]*output[jrow-1][icol]
                +AllVar.ArowP[jrow][icol]*output[jrow+1][icol]+AllVar.AcolM[jrow][icol]*output[jrow][icol-1]
                +AllVar.AcolP[jrow][icol]*output[jrow][icol+1];
                RES[jrow][icol]=(RES[jrow][icol]-C[jrow][icol]*RES[jrow+1][icol-1]-LS[jrow][icol]*RES[jrow-1][icol]-LW[jrow][icol]*RES[jrow][icol-1])*LPR[jrow][icol];
            }


            for(int icol=AGrid.L1-1;icol>0;icol--)
                for(int jrow=AGrid.LocalRowSize-2; jrow>0;jrow--)
                {
                    RES[jrow][icol]=RES[jrow][icol]-G[jrow][icol]*RES[jrow-1][icol+1]-UUN[jrow][icol]*RES[jrow+1][icol]-UUE[jrow][icol]*RES[jrow][icol+1];
                    output[jrow][icol]=output[jrow][icol]+RES[jrow][icol];
                }
    }


}
