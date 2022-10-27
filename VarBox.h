
#ifndef VARBOX_H
#define VARBOX_H
#include "Grid.h"
#include "silo.h"
#include<string>
#include<random>
#include <sstream>
using namespace std;

class VarBox {
public:
    VarBox(Grid& grid);

    VarBox(const VarBox& orig);

    virtual ~VarBox();

    void SetLamda(const double& a );



    void SetMu(const double& a );

    void SetRho(const double& a );

    void SetCp(const double& a );

    void UpdateBoudaryInsulation(vector<vector<double>> & in);

    void EvaluateResidual(vector<vector<double>> & in, double& out );

    void InitUVField_LidDriven();

    void GenerateVar();

    void SUD(vector<vector<double>> & in );

    void outputsilo(const int& step);

    string d2string(const double a);

    int SolveNtime=5;

    double Dt=5E30;
    int MaxTimeStep=1;
    double EOTHER=2.0;
    double EUV=2.0;
    int EveryOutPut=25;
    int BreakStep=122200;

    double Ttol=1e-7;
    double UVtol=1e-7;
    double Mtol=1e-7;

    double Gx=0.0;
    double Gy=-9.8;
    double Beta=0.0;
    double Treference=0.0;

    void Show(vector<vector<double>> & in);
 
private:

    Grid& AGrid;
    vector<vector<double>> ArowM;
    vector<vector<double>> ArowP;
    vector<vector<double>> AcolM;
    vector<vector<double>> AcolP;
    vector<vector<double>> Ap;
    vector<vector<double>> Con;

    vector<vector<double>> UArowM;
    vector<vector<double>> UArowP;
    vector<vector<double>> UAcolM;
    vector<vector<double>> UAcolP;
    vector<vector<double>> UAp;
    vector<vector<double>> UCon;


    vector<vector<double>> VArowM;
    vector<vector<double>> VArowP;
    vector<vector<double>> VAcolM;
    vector<vector<double>> VAcolP;
    vector<vector<double>> VAp;
    vector<vector<double>> VCon;


    vector<vector<double>> MIU;   // from 1 to L1
    vector<vector<double>> MIV;


    vector<vector<double>> U;
    vector<vector<double>> V;
    vector<vector<double>> T;
    vector<vector<double>> P;
    vector<vector<double>> Pc;

    vector<vector<double>> Uold1;
    vector<vector<double>> Vold1;
    vector<vector<double>> Told1;



    vector<vector<double>> Rho;
    vector<vector<double>> Rhoold1;
    vector<vector<double>> Mu;
    vector<vector<double>> Lamda;
    vector<vector<double>> Cp;
    vector<vector<double>> RhoCp;

    vector<vector<double>> Uhat;
    vector<vector<double>> Vhat;
    vector<vector<double>> Du;
    vector<vector<double>> Dv;

    vector<vector<double>> Tu;
    vector<vector<double>> Tv;
    vector<vector<double>> Tc;

    vector<double>UResVec;
    vector<double>VResVec;
    vector<double>MResVec;
    vector<double>TResVec;
    int  NTIMES=3;

    double Time=0.0;
    double TResidual=0.0;
    double CResidual=0.0;
    double PhiResidual=0.0;
    double UResidual=0.0;
    double VResidual=0.0;
    double MResidual=0.0;

    friend class Temperature;
    friend class Solver;
    friend class U;
    friend class V;
    friend class Pressure;
};

#endif /* VARBOX_H */

