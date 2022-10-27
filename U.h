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

#ifndef U_H
#define U_H
#include "Grid.h"

#include "VarBox.h"

#include <cmath>
class U
{
public:
    U(Grid& grid, VarBox& allvar):AGrid(grid),AllVar(allvar)
    {}
    U ( const U& other );
    ~U();
  void Setup();

 void GenerateUhat();

private:

    Grid& AGrid;

    VarBox& AllVar;
};

#endif // U_H
