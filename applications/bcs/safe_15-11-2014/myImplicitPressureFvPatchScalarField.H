/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::myImplicitPressureFvPatchScalarField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    by its magnitude.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        alpha     	| reference value         | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            myHeatFlux;
        alpha        	-10;           // Faktor aus Gleichung
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing Scalars

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    myImplicitPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myImplicitPressureFvPatchScalarField_H
#define myImplicitPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "secondNormalBoundaryCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class myImplicitPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class myImplicitPressureFvPatchScalarField
:
    public zeroGradientFvPatchScalarField,
    public secondNormalBoundaryCells
{

    // Private data


public:

    //- Runtime type information
    TypeName("myPressure")


    // Constructors

        //- Construct from patch and internal field
        myImplicitPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myImplicitPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myImplicitPressureFvPatchScalarField
        //  onto a new patch
        myImplicitPressureFvPatchScalarField
        (
            const myImplicitPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myImplicitPressureFvPatchScalarField
        (
            const myImplicitPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myImplicitPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myImplicitPressureFvPatchScalarField
        (
            const myImplicitPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myImplicitPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
