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
#include "myHeatFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const fvPatch& Theta,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(Theta, iF),
    alpha(Theta.size()),
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
    alpha(ptf.alpha, mapper),
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
    alpha("alpha", dict, Theta.size()),
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
    alpha(pivpvf.alpha),
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
    alpha(pivpvf.alpha),
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
    alpha.autoMap(m);
    Theta_wall.autoMap(m);
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

    alpha.rmap(tiptf.alpha, addr);
    Theta_wall.rmap(tiptf.Theta_wall, addr);
}

// *** Eigene Member-Function *** //

Foam::label Foam::myHeatFluxFvPatchVectorField::getNeighbour(
	const vector& direction,
	const label& cellLabel,
	const vectorField& cellCenters,
	const labelList& cellNeighbours
)
{

	if ((mag(direction) > 1.001) && (mag(direction) < 0.999))
		Info << "Warnung! Richtung ist nicht von Länge 1!" << endl;
	
	label neighbourLabel = -1;
	scalar maxInnerProduct = -1;

	for (int i=0; i<cellNeighbours.size(); i++){
		vector vec = cellCenters[cellNeighbours[i]] - cellCenters[cellLabel];
		vector normed_vec = vec / mag(vec);
		scalar innerProduct = normed_vec & direction;
		if ( innerProduct > maxInnerProduct ){
			maxInnerProduct = innerProduct;
			neighbourLabel = cellNeighbours[i];
		}
	}
	//Info << "maxInnerProduct: " << maxInnerProduct << endl;
	
	return neighbourLabel;
}

void Foam::myHeatFluxFvPatchVectorField::testPatchNF
(
	const vectorField& normals
)
{
	const vectorField& myNormals = patch().nf();
	for(int i=0; i<10; i++){
		if(normals[i] != myNormals[i]){
			Info << "Nicht Gleich!!!!\t" ;
			Info << normals[i] << " vs. " << myNormals[i] << endl;
		}
	}
}

Foam::fvPatchField<Foam::scalar> Foam::myHeatFluxFvPatchVectorField::interpolateTheta_i
(
	const vectorField& normals
)
{

	const scalar numPatchCells = patch().lookupPatchField<volScalarField, scalar>("Theta").size();		// Anzahl Randzellen

	// Feld zum Speichern der Werte von Theta_i
	const volScalarField& iTheta = db().lookupObject<volScalarField>("Theta");
 	const label patchID = patch().boundaryMesh().findPatchID(patch().name()); 
 	fvPatchField<scalar> Theta_i = iTheta.boundaryField()[patchID]; 

	const fvMesh& mesh = iTheta.mesh();  							// Mesh
	const labelListList& cellNeighbours = mesh.cellCells();	// Nachbarn jeder Zelle des Mesh
	const vectorField& cellCenters = mesh.C();					// Zellmittelpunkte
	//const vectorField& patchNormals = patch().nf();			// NormalenVektoren des Randes
	const vectorField& patchFaceCenters = patch().Cf();		// Feld mit Zellmittelpunkten der Randzellen
	const labelList& patchCellLabels = patch().faceCells();	// Label der Randzellen

	label patchCell, neighbour;
	scalar delta1 = 0, delta2 = 0;
	scalar Theta1 = 0, Theta2 = 0;
	
//	Info << "patchNormals[1] " << patchNormals[1] << endl;
	for(int i=0; i<numPatchCells; i++){

//		vector direction = (-1) * patchNormals[i];
		vector direction = (-1) * normals[i];
		patchCell = patchCellLabels[i];

		/*neighbour = getNeighbour( -patchNormals[i],
										  patchCell,
										  cellCenters,
										  cellNeighbours[patchCell] );
*/
		neighbour = getNeighbour( direction,
										  patchCell,
										  cellCenters,
										  cellNeighbours[patchCell] );
		delta1 = mag(patchFaceCenters[i] - cellCenters[patchCell]);
		delta2 = mag(cellCenters[patchCell] - cellCenters[neighbour]);
		Theta1 = iTheta[patchCell];
		Theta2 = iTheta[neighbour];
		/*if (patchCell == 0){ 
			Info << "patchCell: " << patchCell << "\t neighbour: " << neighbour << endl;
			Info << "Theta1: " << Theta1 << "\t Theta2: " << Theta2 << endl;
			Info << "delta1: " << delta1 << "\t delta2: " << delta2 << endl;
		*/
		// Lineare Interpolation
		Theta_i[i] = Theta1 + delta1 * ( Theta1 - Theta2 ) / delta2;
	}

	return Theta_i;
	
}

void Foam::myHeatFluxFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
	 //testPatchNF( patch().nf() );
	 
	 // Wärmestrom auf patch
	 vectorField& s = *this;
	
	 // Theta im inneren Zellwerte
	 const volScalarField& iTheta = 
			db().lookupObject<volScalarField>("Theta");

	 // Theta_i entspricht den Randwerten :-(
 	 label patchID = patch().boundaryMesh().findPatchID(patch().name()); 
 	 fvPatchField<scalar> Theta_i = iTheta.boundaryField()[patchID];

	 // Theta_i interpolieren
	 Theta_i = interpolateTheta_i( patch().nf() );

	 // Wärmestrom setzen
	 s = alpha * (Theta_i - Theta_wall)* patch().nf();

	 /*
	 for (int i=123; i<124; i++){
		 Info << "Theta_i[" << i << "]: " << Theta_i[i] << endl;
		 Info << "Theta_wall[" << i << "]: " << Theta_wall[i] << endl;
		 Info << "alpha[" << i << "]: " << alpha[i] << endl;
		 Info << "s[" << i << "]: " << s[i] << endl;
	 }
	 */

	 // "updated" auf wahr setzen
    fixedValueFvPatchVectorField::updateCoeffs();
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
