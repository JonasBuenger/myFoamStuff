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
    Foam::mySigmaFvPatchSymmTensorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal symmTensor boundary condition
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
        type            myVelocity;
        alpha        	-10;           // Faktor aus Gleichung
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing Scalars

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    mySigmaFvPatchSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef mySigmaFvPatchSymmTensorField_H
#define mySigmaFvPatchSymmTensorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class mySigmaFvPatchSymmTensorField Declaration
\*---------------------------------------------------------------------------*/

class mySigmaFvPatchSymmTensorField
:
    public fixedValueFvPatchSymmTensorField
{
    // Private data

public:

    //- Runtime type information
    TypeName("mySigma");


    // Constructors

        //- Construct from patch and internal field
        mySigmaFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mySigmaFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mySigmaFvPatchSymmTensorField
        //  onto a new patch
        mySigmaFvPatchSymmTensorField
        (
            const mySigmaFvPatchSymmTensorField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mySigmaFvPatchSymmTensorField
        (
            const mySigmaFvPatchSymmTensorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new mySigmaFvPatchSymmTensorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mySigmaFvPatchSymmTensorField
        (
            const mySigmaFvPatchSymmTensorField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new mySigmaFvPatchSymmTensorField
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
                const fvPatchSymmTensorField&,
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
