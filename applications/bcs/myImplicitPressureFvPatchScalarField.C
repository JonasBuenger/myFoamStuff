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
#include "myImplicitPressureFvPatchScalarField.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myImplicitPressureFvPatchScalarField::
myImplicitPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    secondNormalBoundaryCells(p)
//    alpha(p.size())
{
//	Info << "myPressure-Constructor 1" << endl;
}

Foam::myImplicitPressureFvPatchScalarField::
myImplicitPressureFvPatchScalarField
(
    const myImplicitPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(p, iF),
    secondNormalBoundaryCells(p)
 //   alpha(ptf.alpha, mapper)
{
//	Info << "myPressure-Constructor 2" << endl;
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


Foam::myImplicitPressureFvPatchScalarField::
myImplicitPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    secondNormalBoundaryCells(p)
    //fvPatchField<scalar>(p, iF)//,
//    alpha("alpha", dict, p.size())
{
//	Info << "myPressure-Constructor 3" << endl;
    //fvPatchScalarField::operator=(alpha*patch().nf());
}


Foam::myImplicitPressureFvPatchScalarField::
myImplicitPressureFvPatchScalarField
(
    const myImplicitPressureFvPatchScalarField& pivpvf
)
:
    zeroGradientFvPatchScalarField(pivpvf),
    secondNormalBoundaryCells(pivpvf.patch())
//    alpha(pivpvf.alpha)
{
//	Info << "myPressure-Constructor 4" << endl;
}


Foam::myImplicitPressureFvPatchScalarField::
myImplicitPressureFvPatchScalarField
(
    const myImplicitPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(pivpvf, iF),
    secondNormalBoundaryCells(pivpvf.patch())
//    alpha(pivpvf.alpha)
{
//	Info << "myPressure-Constructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myImplicitPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchScalarField::autoMap(m);
//    alpha.autoMap(m);
}


void Foam::myImplicitPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    //const myImplicitPressureFvPatchScalarField& tiptf =
    //    refCast<const myImplicitPressureFvPatchScalarField>(ptf);

//    alpha.rmap(tiptf.alpha, addr);
}

void Foam::myImplicitPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField& boundaryField = *this;

    const volScalarField& p = db().lookupObject<volScalarField>("p");
    const vectorField pGradient = fvc::grad(p);          // klären, wie gradient-methode festgelegt werden kann
    tmp<vectorField> tpatchDeltas = patch().delta();
    const vectorField& patchDeltas = tpatchDeltas();

    for(int i=0; i<patch().size(); i++){

        const label faceCell = patch().faceCells()[i];

        boundaryField[i] = p[faceCell]
                           + ( pGradient[faceCell] & patchDeltas[i] );

    }

    fvPatchField<scalar>::updateCoeffs();

}
//---------------------------------------------------------------------------//

void Foam::myImplicitPressureFvPatchScalarField::write(Ostream& os) const
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
        myImplicitPressureFvPatchScalarField
    );
}

// ************************************************************************* //
