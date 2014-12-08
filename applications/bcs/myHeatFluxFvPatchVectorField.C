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
    const fvPatch& s,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(s, iF),
    ThetaWall(s.size())
{
//	Info << "myHeatFlux-Konstructor 1" << endl;
}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& ptf,
    const fvPatch& s,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(s, iF),
    ThetaWall(ptf.ThetaWall, mapper)
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
    const fvPatch& s,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(s, iF),
    ThetaWall("ThetaWall", dict, s.size())
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
    ThetaWall(pivpvf.ThetaWall)

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
    ThetaWall(pivpvf.ThetaWall)
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
    ThetaWall.autoMap(m);
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

    ThetaWall.rmap(tiptf.ThetaWall, addr);
    //g0.rmap(tiptf.gamma,addr);
}

void Foam::myHeatFluxFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    // DIRECTIONS
    tmp<vectorField> tnormals = this->patch().nf();
    const vectorField& normals = tnormals();
    vectorField tangents1(normals.size());
    vectorField tangents2(normals.size());

    forAll(normals, i){

        const vector n = normals[i];
        vector& t1 = tangents1[i];
        vector& t2 = tangents2[i];
        const scalar tol = 1e-20;

        if(mag(n.x()) > tol){  // if not n_x = 0

            t1 = vector(-n.y(), n.x(), 0);
            t1 = 1/mag(t1) * t1;
            t2 = n ^ t1;   // cross-product

        }
        else if(n.y() > tol ){ // if n_x = 0 AND not n_y = 0

            t1 = vector(0, -n.z(), n.y());
            t1 = 1/mag(t1) * t1;
            t2 = n ^ t1;   // cross-product

        }
        else{   // n = (0 0 1) --> n_z != 0

            t1 = vector(1,0,0);
            t2 = vector(0,1,0);

        }

    }


    // CONSTANTS
    const dictionary& boundaryConstants = db().lookupObject<IOdictionary>("boundaryConstants");
    const dimensionedScalar alpha2(boundaryConstants.lookup("alpha2"));
    const dimensionedScalar beta2(boundaryConstants.lookup("beta2"));
    const dimensionedScalar gamma2(boundaryConstants.lookup("gamma2"));
    const dimensionedScalar delta2(boundaryConstants.lookup("delta2"));


    // BOUNDARY FIELDS
    const volScalarField& Theta = db().lookupObject<volScalarField>("Theta");
    const scalarField& ThetaBoundary = Theta.boundaryField()[patch().index()];

    const volSymmTensorField& sigma = db().lookupObject<volSymmTensorField>("sigma");
    const symmTensorField& sigmaBoundary = sigma.boundaryField()[patch().index()];

    const volVectorField& u = db().lookupObject<volVectorField>("u");
    const vectorField& uBoundary = u.boundaryField()[patch().index()];

    vectorField& sBoundary = *this;
    tmp<vectorField> tsPatchInternalField = this->patchInternalField();
    vectorField sPatchInternalField = tsPatchInternalField();


    // BOUNDARY DELTAS
    scalarField d = 1.0/this->patch().deltaCoeffs();



    // NORMAL HEATFLUX: sn
    scalarField sn =    alpha2.value() * (ThetaBoundary - ThetaWall)
                       + beta2.value() * ( normals & (sigmaBoundary & normals) );

    // TANGENTIAL HEATFLUX
    scalarField st1 = (     ( sPatchInternalField & tangents1 ) - ( d * delta2.value() * (uBoundary & tangents1) )     )
                         /
                            ( 1 + delta2.value() * d );

    scalarField st2 = (     ( sPatchInternalField & tangents2 ) - ( d * delta2.value() * (uBoundary & tangents2) )     )
                         /
                            ( 1 + gamma2.value() * d );


    // RESULTING HEATFLUX ...
    sBoundary = sn * normals + st1 * tangents1 + st2 * tangents2;












/*

    // VARIABLEN/FELDER
    scalarField st(patch().size()),
                sct(patch().size()),
                gt(patch().size());
    vectorField t(patch().size());
    vectorField& s = *this;

    // Berechnung von Richtungen
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
    //const volScalarField& Theta = db().lookupObject<volScalarField>("Theta");
    //const label patchID = patch().boundaryMesh().findPatchID(patch().name());


//    const scalarField& boundaryTheta = Theta.boundaryField()[patchID];

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
/*
    sn = alpha * (boundaryTheta - ThetaWall);

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
/*
    s = sn*n;// + st*t;
*/
    fixedValueFvPatchVectorField::updateCoeffs();     // updated_ = true

}

//---------------------------------------------------------------------------//

void Foam::myHeatFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    ThetaWall.writeEntry("ThetaWall", os);
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
