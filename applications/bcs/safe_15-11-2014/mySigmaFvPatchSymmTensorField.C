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
#include "mySigmaFvPatchSymmTensorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const fvPatch& u,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(u, iF)
{
//	Info << "myVelocity-Konstructor 1" << endl;
}


Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const mySigmaFvPatchSymmTensorField& ptf,
    const fvPatch& u,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(u, iF)
{
//	Info << "myVelocity-Konstructor 2" << endl;
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
 /*  fixedValueFvPatchSymmTensorField::operator=
    (
        symmTensorField
        (
            ptf.alpha*ptf.patch().nf(),
            mapper
        )
    );
 */
}



Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const fvPatch& u,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(u, iF)
{
//	Info << "myVelocity-Konstructor 3" << endl;
//    FvPatchSymmTensorField::operator=(alpha*patch().nf());
}


Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const mySigmaFvPatchSymmTensorField& pivpvf
)
:
    fixedValueFvPatchSymmTensorField(pivpvf)
{
//	Info << "myVelocity-Konstructor 4" << endl;
}

Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const mySigmaFvPatchSymmTensorField& pivpvf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(pivpvf, iF)
{
//	Info << "myVelocity-Konstructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mySigmaFvPatchSymmTensorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchSymmTensorField::autoMap(m);
}


void Foam::mySigmaFvPatchSymmTensorField::rmap
(
    const fvPatchSymmTensorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchSymmTensorField::rmap(ptf, addr);

    const mySigmaFvPatchSymmTensorField& tiptf =
        refCast<const mySigmaFvPatchSymmTensorField>(ptf);

}

void Foam::mySigmaFvPatchSymmTensorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Variablen/ Felder
    const volVectorField& u = db().lookupObject<volVectorField>("u");
    const vectorField& u_b = u.boundaryField()[patch().index()];

    const volSymmTensorField& sigma = db().lookupObject<volSymmTensorField>("sigma");
    //FvPatchSymmTensorField& u_b = *this;
    tmp<symmTensorField> tsigma_b = *this;
    symmTensorField& sigma_b = tsigma_b();
    tmp<vectorField> tn = patch().nf();
    const vectorField normals = tn();
    const scalarField& d = patch().deltaCoeffs();

    const vector e1(1,0,0);
    const vector e2(0,1,0);
    const vector e3(0,0,1);

    const double alpha1 = -0.01;
    const double alpha2 = -0.01;
    const double gamma1 = 0.01;

    // get directions
    vectorField tangents1(normals.size());
    vectorField tangents2(normals.size());
    forAll(normals, i){

        const vector n = normals[i];
        vector& t1 = tangents1[i];
        vector& t2 = tangents2[i];

        if(n.x() != 0){  // if not n_x = 0

            t1 = vector(-n.y(), n.x(), 0);
            t1 = 1/mag(t1) * t1;
            t2 = n ^ t1;   // cross-product

        }
        else if(n.y() != 0){ // if n_x = 0 AND not n_y = 0

            t1 = vector(0, -n.z(), n.y());
            t1 = 1/mag(t1) * t1;
            t2 = n ^ t1;   // cross-product

        }
        else{   // n_z != 0

            t1 = vector(1,0,0);
            t2 = vector(0,1,0);

        }

    }

    forAll(sigma_b, i){

        const vector n = normals[i];
        const vector t1 = tangents1[i];
        const vector t2 = tangents2[i];
        const vector ub = u_b[i];

        tensor T10;
        T10.xx() = n.x(); T10.xy() = t1.x(); T10.xz() = t2.x();
        T10.yx() = n.y(); T10.yy() = t1.y(); T10.yz() = t2.y();
        T10.zx() = n.z(); T10.zy() = t1.z(); T10.zz() = t2.z();

        const tensor T01 = inv(T10);

        const symmTensor sigma_bc_symm = sigma[patch().faceCells()[i]];
        tensor sigma_bc;
        sigma_bc.xx() = sigma_bc_symm.xx();
        sigma_bc.xy() = sigma_bc_symm.xy();
        sigma_bc.xz() = sigma_bc_symm.xz();
        sigma_bc.yx() = sigma_bc_symm.xy();
        sigma_bc.yy() = sigma_bc_symm.yy();
        sigma_bc.yz() = sigma_bc_symm.yz();
        sigma_bc.zx() = sigma_bc_symm.xz();
        sigma_bc.zy() = sigma_bc_symm.yz();
        sigma_bc.zz() = sigma_bc_symm.zz();

        tensor sigma_bc_t = T01 & (sigma_bc & inv(T01));

        sigma_bc_t.xx() = sigma_bc_t.xx() / (1 + gamma1*d[i]);
        sigma_bc_t.xy() = alpha1 * (ub & t1);
        sigma_bc_t.xz() = alpha2 * (ub & t2);
        sigma_bc_t.yx() = alpha1 * (ub & t1);
        sigma_bc_t.zx() = alpha2 * (ub & t2);

        sigma_bc = T10 & (sigma_bc_t & inv(T10));

        //Info << (sigma_bc & n) << "\t" << ub << endl;

        /*
        sigma_b[i].xx() = sigma_bc.xx();
        sigma_b[i].xy() = sigma_bc.xy();
        sigma_b[i].xz() = sigma_bc.xz();
        sigma_b[i].yy() = sigma_bc.yy();
        sigma_b[i].yz() = sigma_bc.yz();
        sigma_b[i].zz() = sigma_bc.zz();
        */

        sigma_b[i] = symm(sigma_bc);

        //Info << "sigma_bc: " << sigma_bc << endl;

        /*

        tensor sigma_b_tmp;
        sigma_b_tmp = sigma_bc
                - ( n & (sigma_bc & t1) ) * sigma_bc * nt1
                - ( n & (sigma_bc & t2) ) * sigma_bc * nt2;

        sigma_b_tmp += alpha1 * (t1 & ub) * t1n;
        sigma_b_tmp += alpha1 * (t2 & ub) * t2n;

        Info << (sigma_b_tmp & t1) << "\t" << ub << endl;
    */
    }


    /*
    for(int i=0; i<patch().size(); i++){

        const label faceCell = patch().faceCells()[i];

        if(n[i] == e1 || n[i] == -e1){
            sigma_b[i].component(symmTensor::XX) = sigma[faceCell].component(symmTensor::XX) /
                                                    (1 + gamma1*d[i]);
            sigma_b[i].component(symmTensor::XY) = alpha1 * u_b[i].component(vector::Y);
            sigma_b[i].component(symmTensor::XZ) = alpha2 * u_b[i].component(vector::Z);
            sigma_b[i].component(symmTensor::YY) = sigma[faceCell].component(symmTensor::YY);
            sigma_b[i].component(symmTensor::YZ) = sigma[faceCell].component(symmTensor::YZ);
            sigma_b[i].component(symmTensor::ZZ) = sigma[faceCell].component(symmTensor::ZZ);
        }
        else if(n[i] == e2 || n[i] == -e2){
            sigma_b[i].component(symmTensor::XX) = sigma[faceCell].component(symmTensor::XX);
            sigma_b[i].component(symmTensor::XY) = alpha1 * u_b[i].component(vector::X);
            sigma_b[i].component(symmTensor::XZ) = sigma[faceCell].component(symmTensor::XZ);
            sigma_b[i].component(symmTensor::YY) = sigma[faceCell].component(symmTensor::YY) /
                                                    (1 + gamma1*d[i]);
            sigma_b[i].component(symmTensor::YZ) = alpha2 * u_b[i].component(vector::Z);
            sigma_b[i].component(symmTensor::ZZ) = sigma[faceCell].component(symmTensor::ZZ);
        }
        else if(n[i] == e3 || n[i] == -e3){
            sigma_b[i].component(symmTensor::XX) = sigma[faceCell].component(symmTensor::XX);
            sigma_b[i].component(symmTensor::XY) = sigma[faceCell].component(symmTensor::XY);
            sigma_b[i].component(symmTensor::XZ) = alpha1 * u_b[i].component(vector::X);
            sigma_b[i].component(symmTensor::YY) = sigma[faceCell].component(symmTensor::YY);
            sigma_b[i].component(symmTensor::YZ) = alpha2 * u_b[i].component(vector::Y);
            sigma_b[i].component(symmTensor::ZZ) = sigma[faceCell].component(symmTensor::ZZ) /
                                                    (1 + gamma1*d[i]);
        }
        else{
            Info << "n[i]: " << n[i] << " ist kein Standard-Einheitsvektor!" << endl;
        }
    }
*/
    //Info << u_b << endl;
    //Info << sigma_b << endl;

    fixedValueFvPatchSymmTensorField::updateCoeffs();     // updated_ = true

}

//---------------------------------------------------------------------------//

void Foam::mySigmaFvPatchSymmTensorField::write(Ostream& os) const
{
    fvPatchSymmTensorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchSymmTensorField,
        mySigmaFvPatchSymmTensorField
    );
}

// ************************************************************************* //
