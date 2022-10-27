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

#ifndef SOLVER_H
#define SOLVER_H
#include "Grid.h"
#include "VarBox.h"
#include "U.h"
#include "V.h"
#include "Pressure.h"
#include "Temperature.h"
class Solver
{
public:

    Solver(Grid& grid, VarBox& allvar,U& aU,V& aV,Pressure& aPressure,Temperature& aTem)
    :AGrid(grid),AllVar(allvar),UEq(aU),VEq(aV),PressureEq(aPressure),TEq(aTem)
    {}

    void SimpleAlgorithm(const int Nstep);
private:

    U& UEq;

    V& VEq;

    Pressure& PressureEq;

    Temperature& TEq;

    double t = 0.0;

    Grid& AGrid;

    VarBox& AllVar;

    void SIP(vector<vector<double>>& output,const int& Ntimes);
    void Jacobi(vector<vector<double>>& output,const int& Ntimes);
    void MSIP(vector<vector<double>>& output,const int& Ntimes);
};

#endif // SOLVER_H
