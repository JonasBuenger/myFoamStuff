// JONAS

#include "fvmGrad.H"
#include "gradScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

    template<class Type>
    tmp<fvMatrix<typename outerProduct<vector, Type>::type> > grad
    (
        GeometricField<Type, fvPatchField, volMesh>& vf
    )
    {
        return fv::gradScheme<Type>::New
        (
            vf.mesh(),
            vf.mesh().schemesDict().gradScheme(vf.name())
        )().fvmGrad(vf);
    }


} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
