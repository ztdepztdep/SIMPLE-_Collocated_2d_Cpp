
#ifndef PRESSURE_H
#define PRESSURE_H
#include "Grid.h"

#include "VarBox.h"
class Pressure {
public:
    Pressure(Grid& grid, VarBox& allvar);

    Pressure(const Pressure& orig);

    virtual ~Pressure();

    void Setup();

    void UpdatePressure();

    void UpdateMassFlux();

private:

    Grid& AGrid;

    VarBox& AllVar;
};
#endif /* PRESSURE_H */

