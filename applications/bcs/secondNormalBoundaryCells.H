#ifndef SECONDNORMALBOUNDARYCELLS
#define SECONDNORMALBOUNDARYCELLS

#include"fvCFD.H"

namespace Foam{

//template<class type>
class secondNormalBoundaryCells{

private:
    List<bool> secondNormalFaceIsOwner;
    List<scalar> secondDeltas;
    List<label> secondFaces;
    List<scalar> secondInternalCoeffs;

public:
   secondNormalBoundaryCells(const fvPatch& p);

   //const Field<type>& valueSecondInternalCoeffs();

};

}

#include "secondNormalBoundaryCells.C"

#endif

