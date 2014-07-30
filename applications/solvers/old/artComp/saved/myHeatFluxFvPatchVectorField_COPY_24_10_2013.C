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
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    alpha(p.size())
{}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    alpha(ptf.alpha, mapper)
{
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
    fixedValueFvPatchVectorField::operator=
    (
        vectorField
        (
            ptf.alpha*ptf.patch().nf(),
            mapper
        )
    );
}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    alpha("alpha", dict, p.size())
{
    fvPatchVectorField::operator=(alpha*patch().nf());
}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    alpha(pivpvf.alpha)
{}


Foam::myHeatFluxFvPatchVectorField::
myHeatFluxFvPatchVectorField
(
    const myHeatFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    alpha(pivpvf.alpha)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myHeatFluxFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    alpha.autoMap(m);
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
}

//---------------------------------------------------------------------------//

//- Update the coefficients associated with the patch field
	 /* 
	 Berechnung des Wärmestroms als
	 
		  		s = alpha * (Theta_w - Theta_i) * normal
	 
	 hierbei ist: - alpha ein konstanter Faktor, 
					  - Theta_w, die Temperatur an der Wand,
					  - Theta_i, die Temperatur im inneren auf die Wand interpoliert
					  - normal, der Einheitsvektor senkrecht auf der Oberfläche nach außen zeigend
	 */
void Foam::myHeatFluxFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

	 // Wärmestrom auf patch
	 vectorField& s = *this;
	 // Theta im inneren Zellwerte
	 const volScalarField& iTheta = 
			db().lookupObject<volScalarField>("Theta");
	 // Randwerte von Theta am patch
	 const fvPatchField<scalar>& Theta_w =     
 			patch().lookupPatchField<volScalarField, scalar>("Theta");
 	 // Theta_i entspricht den Randwerten :-(
 	 label patchID = patch().boundaryMesh().findPatchID(patch().name()); 
 	 fvPatchField<scalar> Theta_i = iTheta.boundaryField()[patchID]; 
	 
	 // Mesh
	 const fvMesh& m = iTheta.mesh();	
	 // Zellnachbarn
	 const labelListList neighbours = m.cellCells();
	 // Zellmittelpunkte
	 const vectorField C = m.C();
	 // Normale Vektoren
	 const vectorField np = patch().nf();
	 // Zellmittelpunkte der PatchZellen
	 const vectorField Cp = patch().Cn();
	 // Label der PatchZellen 
	 const labelList lp = patch().faceCells();

 	 #define get_normal_neighbour(cell_num, direction, nb)			\
		for(int in=0; in<neighbours[cell_num].size(); in++){			\
			vector vec = C[neighbours[cell_num][in]] - C[cell_num];	\
			vector vec_n = vec / mag(vec);									\
			scalar inProd = vec_n & direction;								\
			if ( inProd>maxInProd ){											\
				maxInProd = inProd;												\
				nb = neighbours[cell_num][in];								\
			}																			\
		}

	 const int numCells = 4;
	 scalar Thetas[numCells] = {0};
	 const scalar factors[numCells] = {0.5, 0.5, 0, 0};
	 //scalar deltas[numCells-1];
	 label Neighbours[numCells];
		
	 for (int pc=0; pc<Theta_w.size(); pc++){
		// Initialisierung
		vector dir = -np[pc];
		Neighbours[0] = lp[pc];
		Thetas[0] = iTheta[ Neighbours[0] ];
		// Zell-Label und deltas bestimmen
		for (int c=1; c<numCells; c++){
			scalar maxInProd = -1;
			get_normal_neighbour(Neighbours[c-1], dir, Neighbours[c]);
			Thetas[c] = iTheta[Neighbours[c]];
			//Info << "maxInProd: " << maxInProd << endl;
		}
		
		Theta_i[pc] = 0;
		// Interpolate Theta_i
		for(int j=0; j<numCells; j++){
			Theta_i[pc] += factors[j]*Thetas[j];
		}
		
		/*Info << endl << endl;
		Info << "---------------------------------------------" << endl;
		Info << patch().name() << endl;
		Info << "---------------------------------------------" << endl;
		Info << "np[pc]: " << np[pc] << endl;
		
		// Wert in 1. Zelle
		scalar Theta01 = iTheta[ patch().faceCells()[pc] ];

		// Wert in 2. Zelle
		label normalNeighbour01 = -1;
		scalar maxInProd = -1;
		label normalNeighbour00 = m.findCell(Cp[pc]);				// (Sollte auch geschickter gehen)
		get_normal_neighbour(normalNeighbour00, -np[pc], normalNeighbour01)
		scalar Theta02 = iTheta[normalNeighbour01]; 
		//Info << "normalNeighbour01: " << normalNeighbour01 << endl;
		Info << "maxInProd: " << maxInProd << endl;

		// Wert in 3. Zelle
		label normalNeighbour02 =-1;
		maxInProd = -1;
		get_normal_neighbour(normalNeighbour01, -np[pc], normalNeighbour02)
		scalar Theta03 = iTheta[normalNeighbour02];
		//Info << "normalNeighbour02: " << normalNeighbour02 << endl;
		Info << "maxInProd: " << maxInProd << endl;

		// Wert in 4. Zelle
		label normalNeighbour03 =-1;
		maxInProd = -1;
		get_normal_neighbour(normalNeighbour02, -np[pc], normalNeighbour03)
		scalar Theta04 = iTheta[normalNeighbour03];
		//Info << "normalNeighbour03: " << normalNeighbour03 << endl;
		Info << "maxInProd: " << maxInProd << endl;
	//	Info << Theta01 << " " << Theta02 << " " << Theta03 << " " << Theta04 << endl;
		// Setze Temperatur gleich dem Wert der angrenzenden Zelle

		Theta_i[pc] = 0.5*(Theta01+Theta02);
		//Info << "Theta_i[pc] :" << Theta_i[pc] << endl;
*/
	 }
	
	// Wärmestrom setzen
	s = alpha * (Theta_i - Theta_w) * patch().nf();

//	Info << "patch: " << patch().name() << "; s[1]: " << s[1] << endl;

/*    Interpoliert zwischen einer Punkte-Menge und einem Volumenfeld. Nicht
		vom Inneren auf den Rand :-(

      primitivePatchInterpolation inter(patch().patch());
      tmp<Field<scalar> > Theta_i=inter.pointToFaceInterpolate(iTheta);
      scalarField Theta_i=inter.pointToFaceInterpolate(iTheta);

// 	Berechnet Werte der internen Faces, aber nicht die des Randes :-(
// 	const surfaceScalarField& Theta_i_total = fvc::interpolate(iTheta);
*/
    fixedValueFvPatchVectorField::updateCoeffs();

}

//---------------------------------------------------------------------------//

void Foam::myHeatFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    alpha.writeEntry("alpha", os);
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
