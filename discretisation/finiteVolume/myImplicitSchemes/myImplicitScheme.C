#include "myImplicitScheme.H"

#include "fvCFD.H"

Foam::myImplicitScheme::myImplicitScheme(const volVectorField vf):
    vf_(vf),
    mesh_(vf.mesh())
{
}
