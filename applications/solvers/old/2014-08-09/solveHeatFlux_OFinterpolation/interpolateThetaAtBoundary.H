#ifndef APPROXIMATIONSTHETAATBOUNDARY_H
#define APPROXIMATIONSTHETAATBOUNDARY_H

#include "fvCFD.H"

namespace Foam{

scalarField OFinterpolate(const fvPatch& patch, const objectRegistry& db, scalarField boundaryTheta, const vectorField& normals);

}

#endif // APPROXIMATIONSTHETAATBOUNDARY_H
