/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   VarBox.cpp
 * Author: ztdep
 *
 * Created on January 14, 2017, 10:13 AM
 */

#include "VarBox.h"

VarBox::VarBox(Grid& grid):AGrid(grid)
{
    GenerateVar();
}

VarBox::VarBox(const VarBox& orig):AGrid(orig.AGrid)
{}


VarBox::~VarBox()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox:: GenerateVar()
{
    ArowM.resize(AGrid.LocalRowSize);
    for(auto&&i:ArowM) i.assign(AGrid.L1+1,0.0);

    ArowP.resize(AGrid.LocalRowSize);
    for(auto&&i:ArowP) i.assign(AGrid.L1+1,0.0);

    AcolM.resize(AGrid.LocalRowSize);
    for(auto&&i:AcolM) i.assign(AGrid.L1+1,0.0);

    AcolP.resize(AGrid.LocalRowSize);
    for(auto&&i:AcolP) i.assign(AGrid.L1+1,0.0);

    Con.resize(AGrid.LocalRowSize);
    for(auto&&i:Con) i.assign(AGrid.L1+1,0.0);

    Ap.resize(AGrid.LocalRowSize);
    for(auto&&i:Ap) i.assign(AGrid.L1+1,0.0);

    UArowM.resize(AGrid.LocalRowSize);
    for(auto&&i:UArowM) i.assign(AGrid.L1+1,0.0);

    UArowP.resize(AGrid.LocalRowSize);
    for(auto&&i:UArowP) i.assign(AGrid.L1+1,0.0);

    UAcolM.resize(AGrid.LocalRowSize);
    for(auto&&i:UAcolM) i.assign(AGrid.L1+1,0.0);

    UAcolP.resize(AGrid.LocalRowSize);
    for(auto&&i:UAcolP) i.assign(AGrid.L1+1,0.0);

    UCon.resize(AGrid.LocalRowSize);
    for(auto&&i:UCon) i.assign(AGrid.L1+1,0.0);

    UAp.resize(AGrid.LocalRowSize);
    for(auto&&i:UAp) i.assign(AGrid.L1+1,0.0);

    VArowM.resize(AGrid.LocalRowSize);
    for(auto&&i:VArowM) i.assign(AGrid.L1+1,0.0);

    VArowP.resize(AGrid.LocalRowSize);
    for(auto&&i:VArowP) i.assign(AGrid.L1+1,0.0);

    VAcolM.resize(AGrid.LocalRowSize);
    for(auto&&i:VAcolM) i.assign(AGrid.L1+1,0.0);

    VAcolP.resize(AGrid.LocalRowSize);
    for(auto&&i:VAcolP) i.assign(AGrid.L1+1,0.0);

    VCon.resize(AGrid.LocalRowSize);
    for(auto&&i:VCon) i.assign(AGrid.L1+1,0.0);

    VAp.resize(AGrid.LocalRowSize);
    for(auto&&i:VAp) i.assign(AGrid.L1+1,0.0);

    U.resize(AGrid.LocalRowSize);
    for(auto&&i:U) i.assign(AGrid.L1+1,0.0);

    V.resize(AGrid.LocalRowSize);
    for(auto&&i:V) i.assign(AGrid.L1+1,0.0);

    Uold1.resize(AGrid.LocalRowSize);
    for(auto&&i:Uold1) i.assign(AGrid.L1+1,0.0);

    Vold1.resize(AGrid.LocalRowSize);
    for(auto&&i:Vold1) i.assign(AGrid.L1+1,0.0);

    T.resize(AGrid.LocalRowSize);
    for(auto&&i:T) i.assign(AGrid.L1+1,0.0);

    Told1.resize(AGrid.LocalRowSize);
    for(auto&&i:Told1) i.assign(AGrid.L1+1,0.0);

    P.resize(AGrid.LocalRowSize);
    for(auto&&i:P) i.assign(AGrid.L1+1,0.0);

    Pc.resize(AGrid.LocalRowSize);
    for(auto&&i:Pc) i.assign(AGrid.L1+1,0.0);

    Tc.resize(AGrid.LocalRowSize);
    for(auto&&i:Tc) i.assign(AGrid.L1+1,0.0);

    Rho.resize(AGrid.LocalRowSize);
    for(auto&&i:Rho) i.assign(AGrid.L1+1,1.0);

    RhoCp.resize(AGrid.LocalRowSize);
    for(auto&&i:RhoCp) i.assign(AGrid.L1+1,1.0);

    Rhoold1.resize(AGrid.LocalRowSize);
    for(auto&&i:Rhoold1) i.assign(AGrid.L1+1,1.0);

    Mu.resize(AGrid.LocalRowSize);
    for(auto&&i:Mu) i.assign(AGrid.L1+1,0.0);

    Lamda.resize(AGrid.LocalRowSize);
    for(auto&&i:Lamda) i.assign(AGrid.L1+1,0.0);

    Cp.resize(AGrid.LocalRowSize);
    for(auto&&i:Cp) i.assign(AGrid.L1+1,1.0);

    MIU.resize(AGrid.LocalRowSize);
    for(auto&&i:MIU) i.assign(AGrid.L1+1,0.0);

    MIV.resize(AGrid.LocalRowSize);
    for(auto&&i:MIV) i.assign(AGrid.L1+1,0.0);

    Uhat.resize(AGrid.LocalRowSize);
    for(auto&&i:Uhat) i.assign(AGrid.L1+1,0.0);

    Vhat.resize(AGrid.LocalRowSize);
    for(auto&&i:Vhat) i.assign(AGrid.L1+1,0.0);

    Du.resize(AGrid.LocalRowSize);
    for(auto&&i:Du) i.assign(AGrid.L1+1,0.0);

    Dv.resize(AGrid.LocalRowSize);
    for(auto&&i:Dv) i.assign(AGrid.L1+1,0.0);

    Tu.resize(AGrid.LocalRowSize);
    for(auto&&i:Tu) i.assign(AGrid.L1+1,0.0);

    Tv.resize(AGrid.LocalRowSize);
    for(auto&&i:Tv) i.assign(AGrid.L1+1,0.0);

    cout<<" All Var has been allocated  ";
    system ("hostname");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::InitUVField_LidDriven()
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int jrow=0;jrow<M1+1;jrow++)
    {
        U[jrow][0]=0.0;
        V[jrow][0]=0.0;
        U[jrow][L1]=0.0;
        V[jrow][L1]=0.0;
    }

       for(int icol=0;icol<AGrid.L1+1;icol++)
        {
            U[M1][icol]=1.0;
            V[M1][icol]=0.0;
            U[0][icol]=0.0;
            V[0][icol]=0.0;

        }


    for(int icol=1;icol<L1;icol++)
    {
            MIV[0][icol]=V[0][icol];
            MIV[M1][icol]=V[M1][icol];
            MIV[M1-1][icol]=V[M1][icol];
    }


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        int icol=0;
        MIU[jrow][icol]=U[jrow][icol];
        T[jrow][icol]=Treference+10.0;

        icol=L1;
        MIU[jrow][icol]=U[jrow][icol];
        MIU[jrow][icol-1]=U[jrow][icol];

        T[jrow][icol]=Treference;
    }


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
string VarBox::d2string(const double a) {

    std::ostringstream stm ;
    stm <<std::fixed<<std::setprecision(std::numeric_limits<double>::digits10) <<a ;
    return stm.str() ;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::SetLamda(const double& a )
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int jrow=0;jrow<M1+1;jrow++)
        for(int icol=0;icol<AGrid.L1+1;icol++)
       {

    Lamda[jrow][icol]=a;

     Lamda[0][icol]=0.0;
     Lamda[M1][icol]=0.0;}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::SetMu(const double& a )
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;
    for(int jrow=0;jrow<M1+1;jrow++)
     for(int icol=0;icol<AGrid.L1+1;icol++)
       Mu[jrow][icol]=a;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void VarBox::SetRho(const double& a )
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int jrow=0;jrow<M1+1;jrow++)
        for(int icol=0;icol<AGrid.L1+1;icol++)
            Rho[jrow][icol]=a;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::SetCp(const double& a )
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    for(int jrow=0;jrow<M1+1;jrow++)
        for(int icol=0;icol<AGrid.L1+1;icol++)
        {  Cp[jrow][icol]=a;
           RhoCp[jrow][icol]=a*Rho[jrow][icol];
        }


}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::UpdateBoudaryInsulation(vector<vector<double>> & in)
{
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    if(Jstart==1){
        for(int i=1;i<AGrid.L1;i++)
            in[0][i]=in[1][i];
    }

    if(Jend==AGrid.M1)
    {
        for(int i=1;i<AGrid.L1;i++)
            in[AGrid.M1][i]=in[AGrid.M1-1][i];
    }


    for(int i=0;i<AGrid.LocalRowSize;i++)
    {
        in[i][0]=in[i][1];
        in[i][AGrid.L1]=in[i][AGrid.L1-1];
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::Show(vector<vector<double>> & in)
{
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;

    for(int jrow=Jstart;jrow<Jend;jrow++)
    {   cout<<"jrow: "<<jrow;
        for(int icol=1;icol<L1;icol++)
        {
           cout <<" "<<setw(13)<<fixed<<setprecision(8)<<in[jrow][icol]<<"  ";

        }
         cout<<endl;
    }




}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::EvaluateResidual(vector<vector<double>> & in, double& out )
{
    double temp=0.0;
    double temp1=1E-15;
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;


    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        for(int icol=1;icol<L1;icol++)
        {

          temp=temp+pow(Ap[jrow][icol]*in[jrow][icol]-AcolM[jrow][icol]*in[jrow][icol-1]
            -AcolP[jrow][icol]*in[jrow][icol+1]-ArowM[jrow][icol]*in[jrow-1][icol]
            -ArowP[jrow][icol]*in[jrow+1][icol]-Con[jrow][icol],2.0);
            temp1=temp1+pow(Ap[jrow][icol]*in[jrow][icol],2.0);
        }
    }

    out=sqrt(temp/temp1);
   // cout<<" temp temp1 "<<temp<<"  "<<temp1<<endl;
}






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void VarBox::SUD(vector<vector<double>> & in )
{
    double Diff,Flowe,Floww,Flown,Flows,deffe,deffw,deffs,deffn;
    int L1=AGrid.L1;
    int M1=AGrid.M1;
    int Jstart=AGrid.Jrow0;
    int Jend=AGrid.Jrow1;
    int Jzero=AGrid.Jzero;

    // deffered corecction for the EW directions

    for(int jrow=Jstart;jrow<Jend;jrow++)
    {
        for(int icol=1;icol<L1;icol++)   ///the first inner and the last inner node
        {
            if(icol==1)
            {
                Floww=AGrid.Ycv[jrow]*MIU[jrow-Jzero][0]*Rho[jrow-Jzero][0];
                Flowe=AGrid.Ycv[jrow]*MIU[jrow-Jzero][icol];
                Flowe=Flowe*(AGrid.Fx[icol+1]*Rho[jrow-Jzero][icol+1]+AGrid.Fxm[icol+1]*Rho[jrow-Jzero][icol]);

                deffw=0.0;

                if(Flowe>0.0)
                    deffe=0.0;
                else
                    deffe=-Flowe*0.5*(in[jrow-Jzero][icol+1]-in[jrow-Jzero][icol+2]);

                         Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffw+deffe;
            }

            else if(icol==L1-1)
            {  Flowe=AGrid.Ycv[jrow]*MIU[jrow-Jzero][icol]*Rho[jrow-Jzero][icol];
                Floww=AGrid.Ycv[jrow]*MIU[jrow-Jzero][icol-1];
                Floww=Floww*(AGrid.Fx[icol]*Rho[jrow-Jzero][icol]+AGrid.Fxm[icol]*Rho[jrow-Jzero][icol-1]);

                if(Floww>0.0)
                    deffw=Floww*0.5*(in[jrow-Jzero][icol-1]-in[jrow-Jzero][icol-2]);
                else
                    deffw=0.0;

                deffe=0.0;


                Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffw+deffe;

            }
            else
            {
                Flowe=AGrid.Ycv[jrow]*MIU[jrow-Jzero][icol];
                Flowe=Flowe*(AGrid.Fx[icol+1]*Rho[jrow-Jzero][icol+1]+AGrid.Fxm[icol+1]*Rho[jrow-Jzero][icol]);

                Floww=AGrid.Ycv[jrow]*MIU[jrow-Jzero][icol-1];
                Floww=Floww*(AGrid.Fx[icol]*Rho[jrow-Jzero][icol]+AGrid.Fxm[icol]*Rho[jrow-Jzero][icol-1]);

                if(Floww>0.0)
                    deffw=Floww*0.5*(in[jrow-Jzero][icol-1]-in[jrow-Jzero][icol-2]);
                else
                    deffw=Floww*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero][icol+1]);

                if(Flowe>0.0)
                    deffe=-Flowe*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero][icol-1]);
                else
                    deffe=-Flowe*0.5*(in[jrow-Jzero][icol+1]-in[jrow-Jzero][icol+2]);


                Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffw+deffe;
            }
        }
    }



    // coe for the NS

    // USE CD for the NS direction.

    for(int jrow=Jstart;jrow<Jend;jrow++){
        for(int icol=1;icol<L1;icol++)  {
            if(jrow==1)   {
                Flows=AGrid.Xcv[icol]*MIV[0][icol]*Rho[0][icol] ;
                if(Flows>0.0)
                    deffs=0.0;
                else
                    deffs=0.0;

                Flown=AGrid.Xcv[icol]*MIV[jrow-Jzero][icol];
                Flown=Flown*(AGrid.Fy[jrow+1]*Rho[jrow-Jzero+1][icol]+AGrid.Fym[jrow+1]*Rho[jrow-Jzero][icol] );

                if(Flown>0.0)
                    deffn=-Flown*0.5*(in[jrow-Jzero+1][icol]-in[jrow-Jzero][icol]);
                else
                    deffn=-Flown*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero+1][icol]);


                Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffn+deffs;
            }

            else if(jrow==M1-1)  {
                Flown=AGrid.Xcv[icol]*MIV[jrow-Jzero][icol]*Rho[M1-Jzero][icol];
                if(Flown>0.0)
                    deffn=0.0;
                else
                    deffn=0.0;

                Flows=AGrid.Xcv[icol]*MIV[jrow-Jzero-1][icol];
                Flows=Flows*(AGrid.Fy[jrow]*Rho[jrow-Jzero][icol]+AGrid.Fym[jrow]*Rho[jrow-Jzero-1][icol] );

                if(Flows>0.0)
                    deffs=Flows*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero-1][icol]);
                else
                    deffs=Flows*0.5*(in[jrow-Jzero-1][icol]-in[jrow-Jzero][icol]);



                Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffs+deffn;
            }

            else{

                Flown=AGrid.Xcv[icol]*MIV[jrow-Jzero][icol];
                Flown=Flown*(AGrid.Fy[jrow+1]*Rho[jrow-Jzero+1][icol]+AGrid.Fym[jrow+1]*Rho[jrow-Jzero][icol] );

                if(Flown>0.0)
                    deffn=-Flown*0.5*(in[jrow-Jzero+1][icol]-in[jrow-Jzero][icol]);
                else
                    deffn=-Flown*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero+1][icol]);

                Flows=AGrid.Xcv[icol]*MIV[jrow-Jzero-1][icol];
                Flows=Flows*(AGrid.Fy[jrow]*Rho[jrow-Jzero][icol]+AGrid.Fym[jrow]*Rho[jrow-Jzero-1][icol] );

                if(Flows>0.0)
                    deffs=Flows*0.5*(in[jrow-Jzero][icol]-in[jrow-Jzero-1][icol]);
                else
                    deffs=Flows*0.5*(in[jrow-Jzero-1][icol]-in[jrow-Jzero][icol]);



                Con[jrow-Jzero][icol]= Con[jrow-Jzero][icol]+deffs+deffn;

                //    cout<<jrow-Jzero<<"  "<<icol<<"  "<<ArowM[jrow-Jzero][icol]<<endl;
                //     cout<<jrow-Jzero<<"  "<<icol<<"  "<<ArowP[jrow-Jzero][icol]<<endl;
            }

        }  }

}



