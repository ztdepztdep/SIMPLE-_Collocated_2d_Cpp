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

#include "Grid.h"

Grid::Grid()
{

}


Grid::Grid(const int nx,const int ny,const double& x0, const double& x1,const double& y0, const double& y1)
       :NX(nx),NY(ny),X0(x0),X1(x1),Y0(y0),Y1(y1)
{


}


Grid::Grid(const Grid & other)
{

}

Grid::~Grid()
{

}

void Grid::Generate()
{
    L1=NX+1;
    M1=NY+1;

    XU.resize(NX+2);
    X.resize(NX+2);
    Xdif.resize(NX+2);
    Xcv.resize(NX+2);
    Fx.resize(NX+2);
    Fxm.resize(NX+2);

    YV.resize(NY+2);
    Y.resize(NY+2);
    Ydif.resize(NY+2);
    Ycv.resize(NY+2);
    Fym.resize(NY+2);
    Fy.resize(NY+2);

    XU[1]=X0;
    double DX=(X1-X0)/NX;
    for(int i=2;i<L1+1;i++) XU[i]=XU[i-1]+DX;

    YV[1]=Y0;
    double DY=(Y1-Y0)/NY;
    for(int i=2;i<M1+1;i++) YV[i]=YV[i-1]+DY;

    X[0]=X0;
    for(int i=1;i<L1;i++) X[i]=0.5*(XU[i]+XU[i+1]);
    X[L1]=X1;

    Y[0]=Y0;
    for(int i=1;i<M1;i++) Y[i]=0.5*(YV[i]+YV[i+1]);
    Y[M1]=Y1;


    for(int i=1;i<L1+1;i++) Xdif[i]=X[i]-X[i-1];
    for(int i=1;i<M1+1;i++) Ydif[i]=Y[i]-Y[i-1];

    for(int i=1;i<L1;i++) Xcv[i]=XU[i+1]-XU[i];
    for(int i=1;i<M1;i++) Ycv[i]=YV[i+1]-YV[i];

    for(int i=2;i<L1;i++){
        Fx[i]=0.5*Xcv[i-1]/Xdif[i];
        Fxm[i]=1.0-Fx[i];
    }

    Fx[1]=0.0;
    Fxm[1]=1.0;

    Fx[L1]=1.0;
    Fxm[L1]=0.0;


    for(int i=2;i<M1;i++){
        Fy[i]=0.5*Ycv[i-1]/Ydif[i];
        Fym[i]=1.0-Fy[i];
    }
    Fy[1]=0.0;
    Fym[1]=1.0;

    Fy[M1]=1.0;
    Fym[M1]=0.0;

    // gives each cpu the start of row and col index in the global coordinates in natural ordering.
    Jrow0=1;    // in FVM, +1 means the zero volume control volume in the start and end of computatinoal region
    int RowSize= NY;
    Jrow1=Jrow0+RowSize;  // the last pointJrow1-1==the last inner point.
    Jzero=Jrow0-1;
    LocalRowSize=RowSize+2; // to include two boundary face

    cout<<" Jrow0  Jrow1  LocalRowSize  "<<" "<<Jrow0<<"   "<<Jrow1<<"  "<<LocalRowSize<<endl;


  //  for(auto i: X)    cout<<i<<endl;

}
