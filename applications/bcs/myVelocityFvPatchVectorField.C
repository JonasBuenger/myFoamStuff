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
#include "myVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myVelocityFvPatchVectorField::
myVelocityFvPatchVectorField
(
    const fvPatch& u,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(u, iF)
{
//	Info << "myVelocity-Konstructor 1" << endl;
}


Foam::myVelocityFvPatchVectorField::
myVelocityFvPatchVectorField
(
    const myVelocityFvPatchVectorField& ptf,
    const fvPatch& u,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(u, iF)
{
//	Info << "myVelocity-Konstructor 2" << endl;
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
 /*  fixedValueFvPatchVectorField::operator=
    (
        vectorField
        (
            ptf.alpha*ptf.patch().nf(),
            mapper
        )
    );
 */
}



Foam::myVelocityFvPatchVectorField::
myVelocityFvPatchVectorField
(
    const fvPatch& u,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(u, iF)
{
//	Info << "myVelocity-Konstructor 3" << endl;
//    fvPatchVectorField::operator=(alpha*patch().nf());
    if (dict.found("value"))
    {
        Info << "found" << endl;
        fixedValueFvPatchField<vector>::operator==
        (
             Field<vector>("value", dict, u.size())
        );
    }
}


Foam::myVelocityFvPatchVectorField::
myVelocityFvPatchVectorField
(
    const myVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf)
{
//	Info << "myVelocity-Konstructor 4" << endl;
}

Foam::myVelocityFvPatchVectorField::
myVelocityFvPatchVectorField
(
    const myVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{
//	Info << "myVelocity-Konstructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::myVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const myVelocityFvPatchVectorField& tiptf =
        refCast<const myVelocityFvPatchVectorField>(ptf);

}

void Foam::myVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Variablen/ Felder
    const volVectorField& u = db().lookupObject<volVectorField>("u");
    const volVectorField& uOld = u.oldTime();
    tmp<vectorField> tu_b = *this;
    vectorField& u_b = tu_b();

    const tmp<vectorField> tn = this->patch().nf();
    const vectorField& n = tn();

    tmp<volTensorField> tGradU = fvc::grad(u);
    const volTensorField gradU = tGradU();
    tmp<vectorField> tD = patch().delta();
    const vectorField& D = tD();
    const labelList& bcs = patch().faceCells();

    for(int i=0; i<patch().size(); i++){

        const label bc = bcs[i];
        const vector normal = n[i];
        const tensor Jac_u = gradU[bc];
        const vector u_bc = u[bc];
        const vector Di = D[i];

        vector u_b_extrapolated = u_bc + (Jac_u & Di);

        u_b[i] = u_b_extrapolated - (u_b_extrapolated & normal) * normal;

    }

    fixedValueFvPatchVectorField::updateCoeffs();     // updated_ = true

}

//---------------------------------------------------------------------------//

void Foam::myVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
