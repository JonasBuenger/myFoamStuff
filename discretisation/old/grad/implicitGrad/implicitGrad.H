// JONAS

#ifndef implicitGrad_H
#define implicitGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class implicitGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class implicitGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        implicitGrad(const implicitGrad&);

        //- Disallow default bitwise assignment
        void operator=(const implicitGrad&);


public:

    //- Runtime type information
    TypeName("implicit");


    // Constructors

        //- Construct from Istream
        implicitGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


	// Member Functions

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > grad
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;



        // JONAS
        tmp<fvMatrix<typename outerProduct<vector, Type>::type> >
        fvmGrad
        (
            GeometricField<Type, fvPatchField, volMesh>&
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "implicitGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
