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
    Foam::myHeatFluxFvPatchVectorField

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
    myHeatFluxFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myHeatFluxFvPatchVectorField_H
#define myHeatFluxFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "patchNormalNeighbours.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class myHeatFluxFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class myHeatFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

    patchNormalNeighbours _normalNeighbours;

        scalarField alpha;
        scalarField Theta_wall;
        scalarField gamma;
        vectorField g0;

public:

    //- Runtime type information
    TypeName("myHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        myHeatFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myHeatFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myHeatFluxFvPatchVectorField
        //  onto a new patch
        myHeatFluxFvPatchVectorField
        (
            const myHeatFluxFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myHeatFluxFvPatchVectorField
        (
            const myHeatFluxFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myHeatFluxFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myHeatFluxFvPatchVectorField
        (
            const myHeatFluxFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new myHeatFluxFvPatchVectorField
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
                const fvPatchVectorField&,
                const labelList&
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


	// My functions

        fvPatchField<scalar> interpolateBoundaryTheta
		(
		);
		
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
