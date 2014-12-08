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

#include "fvc.H"
#include "myHeatFluxFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const fvPatch& Theta,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(Theta, iF),
    Theta_wall(Theta.size())
{
//	Info << "myHeatFlux-Konstructor 1" << endl;
}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& ptf,
    const fvPatch& Theta,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(Theta, iF),
    Theta_wall(ptf.Theta_wall, mapper)
{
//	Info << "myHeatFlux-Konstructor 2" << endl;
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

Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const fvPatch& Theta,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(Theta, iF),
    Theta_wall("Theta_wall", dict, Theta.size())
{
//	Info << "myHeatFlux-Konstructor 3" << endl;
//    fvPatchVectorField::operator=(alpha*patch().nf());
}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    Theta_wall(pivpvf.Theta_wall)

{
//	Info << "myHeatFlux-Konstructor 4" << endl;
}

Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    Theta_wall(pivpvf.Theta_wall)
{
//	Info << "myHeatFlux-Konstructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myHeatFluxFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    Theta_wall.autoMap(m);
    //g0.autoMap();
}


void Foam::myHeatFluxFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const myHeatFluxFvPatchVectorField& tiptf =
        refCast<const myHeatFluxFvPatchVectorField>(ptf);

    Theta_wall.rmap(tiptf.Theta_wall, addr);
    //g0.rmap(tiptf.gamma,addr);
}

void Foam::myHeatFluxFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }



    // VARIABLEN/FELDER
    scalarField sn(patch().size()),
                st(patch().size()),
                sct(patch().size()),
                gt(patch().size());
    vectorField t(patch().size()),
                n(patch().size());
    vectorField& s = *this;

    // Berechnung von Richtungen
    n = this->patch().nf();
//    for(int i=0; i<g0.size(); i++){
//        if(mag(g0[i]) == 0)
//            t[i] = vector(0,0,0);
//        else
//            t[i] = g0[i]/mag(g0[i]);
//    }
    for(int i=0; i<n.size(); i++){
        const vector vec = n[i];
        t[i] = vector(vec.y(),-vec.x(),0);
    }

    // Betrag in normalen Richtung n
    const volScalarField& Theta = db().lookupObject<volScalarField>("Theta");
    //const label patchID = patch().boundaryMesh().findPatchID(patch().name());
    const label patchID = patch().index();

    const scalarField& boundaryTheta = Theta.boundaryField()[patchID];

    //const volScalarField& Theta = db().lookupObject<volScalarField>("Theta");
    tmp<volVectorField> tThetaGradient = fvc::grad(Theta);
    const volVectorField ThetaGradient = tThetaGradient();
    //const vectorField ThetaGradient = fvc::grad(internalTheta,"leastSquares");          // klären, wie gradient-methode festgelegt werden kann
    tmp<vectorField> tpatchDeltas = patch().delta();
    const vectorField& patchDeltas = tpatchDeltas();

/*
    for(int i=0; i<patch().size(); i++){
        boundaryTheta[i] = Theta[patch().faceCells()[i]]
                            + ( ThetaGradient[patch().faceCells()[i]] & patchDeltas[i] );
    }
*/
    sn = alpha * (boundaryTheta - Theta_wall);

    // Betrag in tangentiale Richtung t (Approximation durch einseitige FD)
    vector e1(1,0,0);
    vector e2(0,1,0);
    vector e3(0,0,1);
    tmp<vectorField> tsc = this->patchInternalField();
    vectorField sc = tsc();
    scalarField d = 1.0/this->patch().deltaCoeffs();

    gt  = ( (g0 ^ n) & e3 );
    sct = ( (sc ^ n) & e3 );
    st = ( sct - d*gt ) / ( 1 + d*gamma );

    /*
    const volVectorField& s_int = db().lookupObject<volVectorField>("s");
    tmp<volTensorField> tgrads = fvc::grad(s_int);
    const volTensorField grads = tgrads();
    vectorField partialn_s(patch().size());
    const vectorField& normals = this->patch().nf();
    for(int i=0; i<patch().size(); i++){
        partialn_s[i] = grads[patch().faceCells()[i]] & normals[i];
    }

    st = ((partialn_s & t) - gt) / gamma;
*/
    // gesamter Wärmestrom s
    s = sn*n;// + st*t;

    fixedValueFvPatchVectorField::updateCoeffs();     // updated_ = true

}

//---------------------------------------------------------------------------//

void Foam::myHeatFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    alpha.writeEntry("alpha", os);
    Theta_wall.writeEntry("Theta_wall", os);
 	 this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myHeatFluxFvPatchVectorField
    );
}

// ************************************************************************* //
