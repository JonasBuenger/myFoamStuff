#ifndef INTERPOLATIONWEIGHTS
#define INTERPOLATIONWEIGHTS

#include "fvCFD.H"

template<class type>
const tmp<surfaceScalarField> interpolWeights(GeometricField<type, fvPatchField, volMesh>& vf, const word& name);

#include "interpolationWeights.C"

#endif
