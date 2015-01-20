/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "myPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myPressureFvPatchScalarField::
myPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)//,
    //secondNormalBoundaryCells(p)
//    alpha(p.size())
{
//    Info << "myPressure-Constructor 1" << endl;
}

Foam::myPressureFvPatchScalarField::
myPressureFvPatchScalarField
(
    const myPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(p, iF)
 //   alpha(ptf.alpha, mapper)
{
    //Info << "myPressure-Constructor 2" << endl;
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
/*
    fixedValueFvPatchScalarField::operator=
    (
        scalarField
        (
            ptf.alpha*ptf.patch().nf(),
            mapper
        )
    );
*/
}


Foam::myPressureFvPatchScalarField::
myPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //fixedValueFvPatchScalarField(p, iF, dict),
    fixedValueFvPatchScalarField(p, iF)//,
//    alpha("alpha", dict, p.size())
{
//    Info << "myPressure-Constructor 3" << endl;
    //fvPatchScalarField::operator=(alpha*patch().nf());
    if (dict.found("value"))
    {
        Info << "found" << endl;
        fixedValueFvPatchField<scalar>::operator==
        (
             Field<scalar>("value", dict, p.size())
        );
    }
}


Foam::myPressureFvPatchScalarField::
myPressureFvPatchScalarField
(
    const myPressureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf)
//    alpha(pivpvf.alpha)
{
//    Info << "myPressure-Constructor 4" << endl;
}


Foam::myPressureFvPatchScalarField::
myPressureFvPatchScalarField
(
    const myPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF)
//    alpha(pivpvf.alpha)
{
//    Info << "myPressure-Constructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
//    alpha.autoMap(m);
}


void Foam::myPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const myPressureFvPatchScalarField& tiptf =
        refCast<const myPressureFvPatchScalarField>(ptf);

//    alpha.rmap(tiptf.alpha, addr);
}

void Foam::myPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField& boundaryField = *this;

    const volScalarField& p = db().lookupObject<volScalarField>("p");
    const vectorField pGradient = fvc::grad(p);
    tmp<vectorField> tpatchDeltas = patch().delta();
    const vectorField& patchDeltas = tpatchDeltas();

    for(int i=0; i<patch().size(); i++){

        const label faceCell = patch().faceCells()[i];

        boundaryField[i] = p[faceCell]
                           + ( pGradient[faceCell] & patchDeltas[i] );

    }

    fixedValueFvPatchScalarField::updateCoeffs();

}
//---------------------------------------------------------------------------//

void Foam::myPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
//    alpha.writeEntry("alpha", os);
 	 this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myPressureFvPatchScalarField
    );
}

// ************************************************************************* //
