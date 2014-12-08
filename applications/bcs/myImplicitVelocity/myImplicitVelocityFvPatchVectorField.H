/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::myImplicitVelocityFvPatchVectorField

Description
    Foam::myImplicitVelocityFvPatchVectorField

SourceFiles
    myImplicitVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myImplicitVelocityFvPatchVectorField_H
#define myImplicitVelocityFvPatchVectorField_H

#include "implicitExtrapolationFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class implicitExtrapolationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class myImplicitVelocityFvPatchVectorField
:
    public implicitExtrapolationFvPatchVectorField
{

private:

    Field<tensor> internalCoeffsDiagTensor;
    Field<tensor> internalCoeffsLowerTensor;
    Field<tensor> internalCoeffsUpperTensor;

public:

    //- Runtime vector information
    TypeName("myImplicitVelocity")
    // Constructors

        //- Construct from patch and internal field
        myImplicitVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myImplicitVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given myImplicitVelocityFvPatchVectorField<vector>
        //  onto a new patch
        myImplicitVelocityFvPatchVectorField
        (
            const myImplicitVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myImplicitVelocityFvPatchVectorField
        (
            const myImplicitVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new myImplicitVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myImplicitVelocityFvPatchVectorField
        (
            const myImplicitVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new myImplicitVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            void extrapolateUsingOldValues();

            virtual void updateCoeffs();

            virtual tmp<Field<tensor> > valueInternalCoeffsTensor
            (
                const tmp<scalarField>&
            ) const;

            virtual tmp<Field<tensor> > valueInternalCoeffsUpperTensor
            (
                const tmp<scalarField>&
            ) const;

            virtual tmp<Field<tensor> > valueInternalCoeffsLowerTensor
            (
                const tmp<scalarField>&
            ) const;

            tmp<Field<vector> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                notImplemented(type() + "::valueInternalCoeffs(const tmp<Field<scalar> >&)");

                return tmp<Field<vector> > (new Field<vector>(this->size(), pTraits<vector>::zero));
            }

            tmp<Field<vector> > valueInternalCoeffsUpper
            (
                const tmp<scalarField>&
            ) const
            {
                notImplemented(type() + "::valueInternalCoeffsUpper(const tmp<Field<scalar> >&)");

                return tmp<Field<vector> > (new Field<vector>(this->size(), pTraits<vector>::zero));
            }

            tmp<Field<vector> > valueInternalCoeffsLower
            (
                const tmp<scalarField>&
            ) const
            {
                notImplemented(type() + "::valueInternalCoeffsLower(const tmp<Field<scalar> >&)");

                return tmp<Field<vector> > (new Field<vector>(this->size(), pTraits<vector>::zero));
            }

            tmp<Field<vector> > gradientInternalCoeffs() const;

            tmp<Field<vector> > gradientBoundaryCoeffs() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
