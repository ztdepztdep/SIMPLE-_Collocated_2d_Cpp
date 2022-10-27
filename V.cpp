
#include "V.h"

V::V(const V & other):AGrid(other.AGrid),AllVar(other.AllVar)
{ }

V::~V()
{  }


void V::Setup()
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

            AllVar.Con[jrow][icol]=SC+APT*AllVar.Vold1[jrow][icol]-AllVar.Rho[jrow][icol]*AllVar.Gy*AllVar.Beta*(AllVar.T[jrow][icol]-AllVar.Treference);
            AllVar.Ap[jrow][icol]=-AllVar.Ap[jrow][icol]*VOL+AllVar.AcolM[jrow][icol]+AllVar.AcolP[jrow][icol]
                         +AllVar.ArowM[jrow][icol]+AllVar.ArowP[jrow][icol];

            AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]*VOL;
            AllVar.VCon[jrow][icol]=AllVar.Con[jrow][icol];


            AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]+AllVar.Ap[jrow][icol]*AllVar.V[jrow][icol]/AllVar.EUV;
            AllVar.Ap[jrow][icol]=AllVar.Ap[jrow][icol]*( 1.0+1.0/AllVar.EUV );

            AllVar.VArowM[jrow][icol]=AllVar.ArowM[jrow][icol];
            AllVar.VArowP[jrow][icol]=AllVar.ArowP[jrow][icol];
            AllVar.VAcolM[jrow][icol]=AllVar.AcolM[jrow][icol];
            AllVar.VAcolP[jrow][icol]=AllVar.AcolP[jrow][icol];
            AllVar.VAp[jrow][icol]=AllVar.Ap[jrow][icol];

            double PN=AGrid.Fy[jrow+1]*AllVar.P[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.P[jrow][icol] ;
            double PS=AGrid.Fy[jrow]*AllVar.P[jrow][icol]+AGrid.Fym[jrow]*AllVar.P[jrow-1][icol];
            AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]-(PN-PS)*AGrid.Xcv[icol];

            //  cout<<jrow<<"  "<<icol<<"  "<<AllVar.Ap[jrow][icol]<<"  "<<AllVar.Con[jrow][icol]<<endl;
        }
    }
    AllVar.SUD(AllVar.V);
}



void V::GenerateVhat()
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

            AllVar.Vhat[jrow][icol]=AllVar.VAcolM[jrow][icol]*AllVar.V[jrow][icol-1]
            +AllVar.VAcolP[jrow][icol]*AllVar.V[jrow][icol+1]
            +AllVar.VArowM[jrow][icol]*AllVar.V[jrow-1][icol]
            +AllVar.VArowP[jrow][icol]*AllVar.V[jrow+1][icol]
            +AllVar.VCon[jrow][icol];

            AllVar.Vhat[jrow][icol]=AllVar.Vhat[jrow][icol]/AllVar.VAp[jrow][icol];
            AllVar.Dv[jrow][icol]=AGrid.Xcv[icol]/AllVar.VAp[jrow][icol];
            AllVar.Tv[jrow][icol]=-0.5*AllVar.Rho[jrow][icol]*AllVar.Gy*AllVar.Beta/AllVar.VAp[jrow][icol];
        }
    }


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        for(int icol=1;icol<L1;icol++)
        {
            if(jrow<M1-1) {
                double Vhate=AGrid.Fy[jrow+1]*AllVar.Vhat[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Vhat[jrow][icol] ;
                double Dve=AGrid.Fy[jrow+1]*AllVar.Dv[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Dv[jrow][icol];
                AllVar.MIV[jrow][icol]=Vhate+Dve*(AllVar.P[jrow][icol]-AllVar.P[jrow+1][icol]);
                AllVar.MIV[jrow][icol]=AllVar.MIV[jrow][icol]*( 1.0+1.0/AllVar.EUV );
            }

        }
    }

}
