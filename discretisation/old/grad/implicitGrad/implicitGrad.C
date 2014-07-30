// JONAS

#include "implicitGrad.H"
#include "fvMatrices.H"
#include "fvMatrix.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
implicitGrad<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    Info << "Diese Funktion ist nicht sinnvoll!!!" << endl;

    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tfGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    return tfGrad;
}


// JONAS
template<class Type>
tmp<fvMatrix<typename outerProduct<vector, Type>::type> >
implicitGrad<Type>::fvmGrad
(
     GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    Info << "implicitGrad::fvmGrad" << endl;

    //typedef typename outerProduct<Type, vector>::type GradType;

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()/dimLength
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    return tfvm;

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
