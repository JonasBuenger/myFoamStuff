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
#include "surfaceInterpolate.H"
#include "myTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

#include "primitivePatchInterpolation.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTemperatureFvPatchScalarField::
myTemperatureFvPatchScalarField
(
    const fvPatch& Theta,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(Theta, iF),
    _normalNeighbours(Theta)
//    alpha(Theta.size())
{
//	Info << "myTemperature-Constructor 1" << endl;
}


Foam::myTemperatureFvPatchScalarField::
myTemperatureFvPatchScalarField
(
    const myTemperatureFvPatchScalarField& ptf,
    const fvPatch& Theta,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(Theta, iF),
    _normalNeighbours(Theta)
 //   alpha(ptf.alpha, mapper)
{
//	Info << "myTemperature-Constructor 2" << endl;
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


Foam::myTemperatureFvPatchScalarField::
myTemperatureFvPatchScalarField
(
    const fvPatch& Theta,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(Theta, iF, dict),
    _normalNeighbours(Theta)//,
    //fixedValueFvPatchScalarField(Theta, iF)//,
//    alpha("alpha", dict, Theta.size())
{
//	Info << "myTemperature-Constructor 3" << endl;
    //fvPatchScalarField::operator=(alpha*patch().nf());
}


Foam::myTemperatureFvPatchScalarField::
myTemperatureFvPatchScalarField
(
    const myTemperatureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    _normalNeighbours(pivpvf._normalNeighbours)
//    alpha(pivpvf.alpha)
{
//	Info << "myTemperature-Constructor 4" << endl;
}


Foam::myTemperatureFvPatchScalarField::
myTemperatureFvPatchScalarField
(
    const myTemperatureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    _normalNeighbours(pivpvf._normalNeighbours)
//    alpha(pivpvf.alpha)
{
//	Info << "myTemperature-Constructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
//    alpha.autoMap(m);
}


void Foam::myTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    //const myTemperatureFvPatchScalarField& tiptf =
    //    refCast<const myTemperatureFvPatchScalarField>(ptf);

//    alpha.rmap(tiptf.alpha, addr);
}

// *** Eigene Member-Function *** //

Foam::fvPatchField<Foam::scalar> Foam::myTemperatureFvPatchScalarField::interpolateBoundaryTheta
(
)
{

    const int order = 1;
    const volScalarField& internalTheta = db().lookupObject<volScalarField>("Theta");
    const label patchID = patch().boundaryMesh().findPatchID(patch().name());
    fvPatchField<scalar> boundaryTheta = internalTheta.boundaryField()[patchID];

    for(int i=0; i<boundaryTheta.size(); i++){
        simpleMatrix<scalar> interpolationMatrix(_normalNeighbours.get_interpolationMatrix(i, order));
        for(int j=0; j<order; j++){
            interpolationMatrix.source()[j] = internalTheta[_normalNeighbours.get_ID(i,j)];
        }
        boundaryTheta[i] = interpolationMatrix.solve()[0];
    }

    return boundaryTheta;
}


void Foam::myTemperatureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

     // Temperatur
     scalarField& boundaryTheta = *this;
     boundaryTheta = interpolateBoundaryTheta();

     // "updated" auf wahr setzen
    fixedValueFvPatchScalarField::updateCoeffs();
}

//---------------------------------------------------------------------------//

void Foam::myTemperatureFvPatchScalarField::write(Ostream& os) const
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
        myTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
