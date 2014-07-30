#ifndef MYIMPLICITSCHEME_H
#define MYIMPLICITSCHEME_H

#include "fvMesh.H"
#include "GeometricFields.H"
#include "fvMatrices.H"

namespace Foam{


class myImplicitScheme
{
public:
    const volVectorField& vf_;
    const fvMesh& mesh_;
public:
    myImplicitScheme(const volVectorField);
};


}
#endif // MYIMPLICITSCHEME_H
