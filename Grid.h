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

#ifndef GRID_H
#define GRID_H
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include <algorithm>
#include <iterator>
using namespace std;
class Grid
{
public:
    Grid();
    Grid(const int nx,const int ny,const double& x0, const double& x1,const double& y0, const double& y1);
    Grid ( const Grid& other );
   ~Grid();

   void Generate();

private:

    int L1;

    int M1;

    void setX0X1(const double& x0, const double& x1){X0=x0;X1=x1;}

    void setY0Y1(const double& y0, const double& y1){Y0=y0;Y1=y1;}

    int NX=10;
    int NY=10;

    double X0=0.0;
    double X1=1.0;
    double Y0=0.0;
    double Y1=1.0;

    int Jrow0;
    int Jrow1;
    int Jzero;
    int LocalRowSize;

    vector<double> XU; // start from 0 to L1
    vector<double> YV;

    vector<double> X;
    vector<double> Y;

    vector<double> Xdif;
    vector<double> Ydif;

    vector<double> Xcv;
    vector<double> Ycv;

    vector<double> Fxm;
    vector<double> Fx;

    vector<double> Fym;
    vector<double> Fy;
    friend class Temperature;
    friend class VarBox;
    friend class Solver;
    friend class U;
    friend class V;
    friend class Pressure;
};

#endif // GRID_H
