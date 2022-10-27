#include "Pressure.h"

Pressure::Pressure(Grid& grid, VarBox& allvar):AGrid(grid),AllVar(allvar)
{}

Pressure::Pressure(const Pressure& orig):AGrid(orig.AGrid),AllVar(orig.AllVar)
{}

Pressure::~Pressure()
{}

void Pressure::Setup()
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;
    double ARho,Flow;
    double Csum=0.0;

    for(int jrow=Jstart;jrow<Jend;jrow++)
        for(int icol=1;icol<L1;icol++)
         AllVar.Con[jrow][icol]=0.0;

        for(int icol=1;icol<L1;icol++)
        {
            Flow=AGrid.Xcv[icol]*AllVar.MIV[0][icol]*AllVar.Rho[0][icol];
            AllVar.Con[1][icol]= AllVar.Con[1][icol]+Flow;
            AllVar.ArowM[1][icol]=0.0;
        }

        for(int jrow=Jstart;jrow<Jend;jrow++)
        {
            ARho= AGrid.Ycv[jrow]*AllVar.Rho[jrow][0];
            Flow=ARho*AllVar.MIU[jrow][0];

            AllVar.Con[jrow][1]= AllVar.Con[jrow][1]+Flow;
            AllVar.AcolM[jrow][1]=0.0;

            for(int icol=1;icol<L1;icol++)   ///the first inner and the last inner node
            {
                if(icol==L1-1)
                {
                    ARho=AGrid.Ycv[jrow]*AllVar.Rho[jrow][L1];
                    Flow=ARho*AllVar.MIU[jrow][icol];
                    AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]-Flow;
                    AllVar.AcolP[jrow][icol]=0.0;
                }
                else
                {
                    ARho=AGrid.Ycv[jrow];
                    ARho=ARho*(AGrid.Fx[icol+1]*AllVar.Rho[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Rho[jrow][icol]);

                    Flow=ARho*AllVar.MIU[jrow][icol];
                    AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]-Flow;
                    AllVar.Con[jrow][icol+1]= AllVar.Con[jrow][icol+1]+Flow;
                    double Due=AGrid.Fx[icol+1]*AllVar.Du[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Du[jrow][icol];
                    AllVar.AcolP[jrow][icol]=ARho*Due;;
                    AllVar.AcolM[jrow][icol+1]=AllVar.AcolP[jrow][icol];
                }

                if(jrow==M1-1)
                {
                        ARho=AGrid.Xcv[icol]*AllVar.Rho[M1][icol];
                        Flow=ARho*AllVar.MIV[jrow][icol];
                        AllVar.ArowP[jrow][icol]=0.0;
                        AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]-Flow;
                }
                else
                {

                    ARho=AGrid.Xcv[icol];
                    ARho=ARho*(AGrid.Fy[jrow+1]*AllVar.Rho[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Rho[jrow][icol] );
                    Flow=ARho*AllVar.MIV[jrow][icol];
                    AllVar.Con[jrow][icol]= AllVar.Con[jrow][icol]-Flow;
                    AllVar.Con[jrow+1][icol]= AllVar.Con[jrow+1][icol]+Flow;
                    double Dve=AGrid.Fy[jrow+1]*AllVar.Dv[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Dv[jrow][icol];
                    AllVar.ArowP[jrow][icol]=ARho*Dve;
                    AllVar.ArowM[jrow+1][icol]=AllVar.ArowP[jrow][icol];
                }

                AllVar.Ap[jrow][icol]=AllVar.AcolM[jrow][icol]+AllVar.AcolP[jrow][icol]
                                            +AllVar.ArowM[jrow][icol]+AllVar.ArowP[jrow][icol];
                AllVar.Pc[jrow][icol]=0.0;

                Csum=max(Csum, abs(AllVar.Con[jrow][icol]));
                AllVar.MResidual=Csum;
                //    cout<<jrow<<"  "<<icol<<"  "<<AllVar.AcolP[jrow][icol]<<endl;
            }
        }

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Pressure::UpdatePressure()
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;
    for(int jrow=Jstart;jrow<Jend;jrow++)
        for(int icol=1;icol<L1;icol++)
        {
            AllVar.P[jrow][icol]= AllVar.P[jrow][icol]+0.3*AllVar.Pc[jrow][icol];

            AllVar.P[jrow][0]=AllVar.P[jrow][1];
            AllVar.P[jrow][L1]=AllVar.P[jrow][L1-1];
            AllVar.P[0][icol]= AllVar.P[1][icol];
            AllVar.P[M1][icol]= AllVar.P[M1-1][icol];
        }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Pressure::UpdateMassFlux()
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int jrow=Jstart;jrow<Jend;jrow++)
        for(int icol=1;icol<L1-1;icol++)
        {
            double PE=AllVar.Pc[jrow][icol+1];
            double PW=AllVar.Pc[jrow][icol];
            double Due=AGrid.Fx[icol+1]*AllVar.Du[jrow][icol+1]+AGrid.Fxm[icol+1]*AllVar.Du[jrow][icol];
            AllVar.MIU[jrow][icol]=AllVar.MIU[jrow][icol]-(PE-PW)*Due;

        }

        for(int jrow=Jstart;jrow<Jend-1;jrow++)
            for(int icol=1;icol<L1;icol++)
            {
                      double PN=AllVar.Pc[jrow+1][icol] ;
                    double PS=AllVar.Pc[jrow][icol];
                    double Dve=AGrid.Fy[jrow+1]*AllVar.Dv[jrow+1][icol]+AGrid.Fym[jrow+1]*AllVar.Dv[jrow][icol];
                    AllVar.MIV[jrow][icol]= AllVar.MIV[jrow][icol]-(PN-PS)*Dve;
            }



}
