#include <iostream>
#include "Grid.h"
#include "VarBox.h"
#include "Solver.h"
#include "U.h"
#include "V.h"
#include "Pressure.h"
int main(int argc, char **argv)
{
    Grid  FvmMesh(100,100,0.0,1.0,0.0,1.0) ;

    FvmMesh.Generate();
    VarBox  AllVar(FvmMesh);

    AllVar.SetRho(1.0);
    AllVar.SetMu(0.001);
    AllVar.SetCp(1);
    AllVar.SetLamda(1);


    AllVar.Gx=0.0;
    AllVar.Gy=-0.0;
    AllVar.Beta=0.00333;
    AllVar.Treference=27.0;

    AllVar.SolveNtime=6;
    AllVar.BreakStep=6000;
    AllVar.UVtol=1E-6;

    double uvrelax=0.85;
    AllVar.EUV=uvrelax/(1.0-uvrelax);
    double trelax=0.999;
    AllVar.EOTHER=trelax/(1.0-trelax);

    AllVar.InitUVField_LidDriven();

    U U(FvmMesh,AllVar);
    V V(FvmMesh,AllVar);
    Pressure Pressure(FvmMesh,AllVar);

    Temperature Temperature(FvmMesh,AllVar);

    Solver S(FvmMesh,AllVar,U,V,Pressure,Temperature);
    S.SimpleAlgorithm(1);
    return 0;
}
