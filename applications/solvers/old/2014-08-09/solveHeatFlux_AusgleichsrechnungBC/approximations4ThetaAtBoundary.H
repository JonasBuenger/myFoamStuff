#ifndef APPROXIMATIONS4THETAATBOUNDARY_H
#define APPROXIMATIONS4THETAATBOUNDARY_H

#include "fvCFD.H"

namespace Foam{

scalarField LeastSquareApproximation(const int type, const fvPatch& patch, const objectRegistry& db, scalarField boundaryTheta);

}

#endif // APPROXIMATIONS4THETAATBOUNDARY_H
