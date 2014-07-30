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
#include "simpleMatrix.H"

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
    Theta_wall(Theta.size()),
    gamma(Theta.size()),
    gt(Theta.size())
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
    Theta_wall(ptf.Theta_wall, mapper),
    gamma(ptf.gamma, mapper),
    gt(ptf.gt, mapper)
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
    Theta_wall("Theta_wall", dict, Theta.size()),
    gamma("gamma", dict, Theta.size()),
    gt("gt", dict, Theta.size())
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
    Theta_wall(pivpvf.Theta_wall),
    gamma(pivpvf.gamma),
    gt(pivpvf.gt)
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
    Theta_wall(pivpvf.Theta_wall),
    gamma(pivpvf.gamma),
    gt(pivpvf.gt)
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
    gamma.autoMap(m);
    //gt.autoMap();
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
    gamma.rmap(tiptf.gamma, addr);
    //gt.rmap(tiptf.gamma,addr);
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

	for(int i=0; i<numPatchCells; i++){

		vector direction = (-1) * normals[i];
		patchCell = patchCellLabels[i];

		neighbour = getNeighbour( direction,
										  patchCell,
										  cellCenters,
										  cellNeighbours[patchCell] );
		//Info << "patchFaceCenter: " << patchFaceCenters[i] << endl;
		//Info << "cellCenter:" << cellCenters[patchCell] << endl;
		//Info << "cellCenter:" << cellCenters[neighbour] << endl;
		delta1 = mag(patchFaceCenters[i] - cellCenters[patchCell]);
		delta2 = mag(cellCenters[patchCell] - cellCenters[neighbour]);
		Theta1 = iTheta[patchCell];
		Theta2 = iTheta[neighbour];
		//if (patchCell == numPatchCells){ 
			//Info << "patchCell: " << patchCell << "\t neighbour: " << neighbour << endl;
			//Info << "Theta1: " << Theta1 << "\t Theta2: " << Theta2 << endl;
			//Info << "delta1: " << delta1 << "\t delta2: " << delta2 << endl;
		//}

		// Lineare Interpolation
		Theta_i[i] = Theta1 + delta1 * ( Theta1 - Theta2 ) / delta2;
	}

	return Theta_i;
	
}

Foam::fvPatchField<Foam::vector> Foam::myHeatFluxFvPatchVectorField::interpolate_normalDerivative_st
(
    const vectorField& normals,
    const scalar& polynomialDegree
)
{

    // Variablen

    // Anzahl Randzellen
    const scalar numPatchCells = patch().lookupPatchField<volScalarField, scalar>("Theta").size();

    // Feld zum Speichern der Werte von normalDerivative_st
    const volVectorField& s = db().lookupObject<volVectorField>("s");
    const label patchID = patch().boundaryMesh().findPatchID(patch().name());
    fvPatchField<vector> partialn_st = s.boundaryField()[patchID];

    // mesh
    const fvMesh& mesh = s.mesh();                          // Mesh
    const labelListList& cellNeighbours = mesh.cellCells();	// Nachbarn jeder Zelle des Mesh
    const vectorField& cellCenters = mesh.C();              // Zellmittelpunkte
    const vectorField& patchFaceCenters = patch().Cf();		// Feld mit Zellmittelpunkten der Randzellen
    const labelList& patchCellLabels = patch().faceCells();	// Label der Randzellen

    // matrix D
    simpleMatrix<vector> D(polynomialDegree);
    simpleMatrix<scalar> Dx(polynomialDegree);
    simpleMatrix<scalar> Dy(polynomialDegree);
    simpleMatrix<scalar> Dz(polynomialDegree);
    vector st_loc;

    // di for temporal storage of distances d between cells and tangential heatfluxes st
    double di = 0;

    // loop over all patch cells
    for(int i=0; i<numPatchCells; i++){

        vector direction = (-1) * normals[i];
        label neighbour = patchCellLabels[i];

        for(int j=0; j<polynomialDegree; j++){
            neighbour = getNeighbour( direction, neighbour, cellCenters, cellNeighbours[neighbour] );
            di = mag(patchFaceCenters[i] - cellCenters[neighbour]);
            for(int c=0; c<polynomialDegree; c++){
                D[j][c] = Dx[j][c] = Dy[j][c] = Dz[j][c] = pow(di,c);
            }
            Info << "s[neighbour]: " << s[neighbour] << endl;
            Info << "direction: " << direction << endl;
            st_loc = s[neighbour] - ( s[neighbour] & direction ) * direction;
            Info << "st_loc" << st_loc << endl;
            Info << "(st_loc & direction): " << (st_loc & direction) << endl;
            D.source()[j] = st_loc;
            Dx.source()[j] = st_loc.x();
            Dy.source()[j] = st_loc.y();
            Dz.source()[j] = st_loc.z();
        }

//        for(int s=0; s<polynomialDegree; s++){
//            vector vec;
//            vec.x() = Dx.solve()[s];
//            vec.y() = Dy.solve()[s];
//            vec.z() = Dz.solve()[s];
//            if(D.solve()[s] != vec)
//                Info << "Nicht Gleich!!!" << endl;
//            else
//                Info << "Gleich!!!" << endl;
//        }

        //Info << D << endl << endl;
        //Info << D.solve() << endl << endl << endl;

        //Info << Dx.solve() << endl << endl << endl;

        partialn_st[i].x() = Dx.solve()[1];
        partialn_st[i].y() = Dy.solve()[1];
        partialn_st[i].z() = Dz.solve()[1];

        // Info << partialn_st[i] << endl;

    }

    return partialn_st;

}


void Foam::myHeatFluxFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
	 
	 // Wärmestrom auf patch
     vectorField sn, st;// = *this;
     vectorField& s = *this;
	
	 // Theta im inneren Zellwerte
     //const volScalarField& iTheta = db().lookupObject<volScalarField>("Theta");

     // Theta_i entspricht den Randwerten
     //label patchID = patch().boundaryMesh().findPatchID(patch().name());
     //fvPatchField<scalar> Theta_i = iTheta.boundaryField()[patchID];

     // in normalen Richtung
     fvPatchField<scalar> Theta_i = interpolateTheta_i( patch().nf() );
     sn = alpha * (Theta_i - Theta_wall)* patch().nf();

     // in tangentialer Richtung
     scalar polynomDegree = 3;
     fvPatchField<vector> partialn_st = interpolate_normalDerivative_st (patch().nf(), polynomDegree);
     st = (- partialn_st ) / gamma;
     vector vec;
     vec.x() = 1; vec.y() = 0; vec.z() = 0;
     st = vec - (( vec & patch().nf() ) * patch().nf());
     // gesamter Wärmestrom s
     s = sn + st;

	 //Info << "Ende updateCoeffs: HeatFlux" << endl;
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
