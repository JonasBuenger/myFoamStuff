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
    const fvPatch& sigma,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(sigma, iF)//,
    //ThetaWall(sigma.size())
{
//	Info << "myVelocity-Konstructor 1" << endl;
}


Foam::mySigmaFvPatchSymmTensorField::
mySigmaFvPatchSymmTensorField
(
    const mySigmaFvPatchSymmTensorField& ptf,
    const fvPatch& sigma,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(sigma, iF)//,
    //ThetaWall(ptf.ThetaWall, mapper)
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
    const fvPatch& sigma,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(sigma, iF)//,
    //ThetaWall("ThetaWall", dict, sigma.size())
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
    fixedValueFvPatchSymmTensorField(pivpvf)//,
    //ThetaWall(pivpvf.ThetaWall)
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
    fixedValueFvPatchSymmTensorField(pivpvf, iF)//,
    //ThetaWall(pivpvf.ThetaWall)
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
    //ThetaWall.autoMap(m);
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

    //ThetaWall.rmap(tiptf.ThetaWall, addr);
}

void Foam::mySigmaFvPatchSymmTensorField::updateCoeffs()
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
        const scalar tol = 0;

        if(mag(n.x()) > tol){  // if not n_x = 0

            t1 = vector(-n.y(), n.x(), 0);
            t1 = 1/mag(t1) * t1;
            t2 = n ^ t1;   // cross-product

        }
        else if(mag(n.y()) > tol ){ // if n_x = 0 AND not n_y = 0

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
    const dimensionedScalar alpha1(boundaryConstants.lookup("alpha1"));
    const dimensionedScalar beta1(boundaryConstants.lookup("beta1"));
    const dimensionedScalar gamma1(boundaryConstants.lookup("gamma1"));
    const dimensionedScalar delta1(boundaryConstants.lookup("delta1"));

    // BOUNDARY FIELDS
    //const volScalarField& Theta = db().lookupObject<volScalarField>("Theta");
    //const scalarField& ThetaBoundary = Theta.boundaryField()[patch().index()];

    const volVectorField& u = db().lookupObject<volVectorField>("u");
    const vectorField& uBoundary = u.boundaryField()[patch().index()];

    //const volVectorField& s = db().lookupObject<volVectorField>("s");
    //const vectorField& sBoundary = s.boundaryField()[patch().index()];

    const volSymmTensorField& sigmaInternalField = db().lookupObject<volSymmTensorField>("sigma");

    symmTensorField& sigmaBoundary = *this;
    tmp<symmTensorField> tsigmaPatchInternalField = this->patchInternalField();
    symmTensorField sigmaPatchInternalField = tsigmaPatchInternalField();


    // BOUNDARY DELTAS
    scalarField d = 1.0/this->patch().deltaCoeffs();
    tmp<vectorField> tpatchDeltas = patch().delta();
    const vectorField& patchDeltas = tpatchDeltas();


    // GRADIENT-FIELD OF SIGMA FOR EXTRAPOLATION
    const volVectorField  grad_sigmaInternalField_XX = fvc::grad(sigmaInternalField.component(symmTensor::XX));
    const volVectorField  grad_sigmaInternalField_XY = fvc::grad(sigmaInternalField.component(symmTensor::XY));
    const volVectorField  grad_sigmaInternalField_XZ = fvc::grad(sigmaInternalField.component(symmTensor::XZ));
    const volVectorField& grad_sigmaInternalField_YX =      grad_sigmaInternalField_XY;
    const volVectorField  grad_sigmaInternalField_YY = fvc::grad(sigmaInternalField.component(symmTensor::YY));
    const volVectorField  grad_sigmaInternalField_YZ = fvc::grad(sigmaInternalField.component(symmTensor::YZ));
    const volVectorField& grad_sigmaInternalField_ZX =      grad_sigmaInternalField_XZ;
    const volVectorField& grad_sigmaInternalField_ZY =      grad_sigmaInternalField_YZ;
    const volVectorField  grad_sigmaInternalField_ZZ = fvc::grad(sigmaInternalField.component(symmTensor::ZZ));


    // SET SIGMA AT BOUNDARY
    forAll(sigmaBoundary, i){

        // generate transformation Matrixes T between localBoundary coordinates (l) and global coordinates (g)
        const vector n = normals[i];
        const vector t1 = tangents1[i];
        const vector t2 = tangents2[i];
        const vector ub = uBoundary[i];
        //const vector sb = sBoundary[i];
        //const scalar Thetab = ThetaBoundary[i];

        tensor T_l2g;
        T_l2g.xx() = n.x(); T_l2g.xy() = t1.x(); T_l2g.xz() = t2.x();
        T_l2g.yx() = n.y(); T_l2g.yy() = t1.y(); T_l2g.yz() = t2.y();
        T_l2g.zx() = n.z(); T_l2g.zy() = t1.z(); T_l2g.zz() = t2.z();
        const tensor T_g2l = inv(T_l2g);


        // extrapolate internal sigma-Field onto boundary
        tensor sigmaExtrapolatedGlobal;

        const label faceCell = patch().faceCells()[i];

        sigmaExtrapolatedGlobal.component(tensor::XX) =  sigmaInternalField[faceCell].component(symmTensor::XX) + ( grad_sigmaInternalField_XX[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::XY) =  sigmaInternalField[faceCell].component(symmTensor::XY) + ( grad_sigmaInternalField_XY[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::XZ) =  sigmaInternalField[faceCell].component(symmTensor::XZ) + ( grad_sigmaInternalField_XZ[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::YX) =  sigmaInternalField[faceCell].component(symmTensor::XY) + ( grad_sigmaInternalField_YX[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::YY) =  sigmaInternalField[faceCell].component(symmTensor::YY) + ( grad_sigmaInternalField_YY[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::YZ) =  sigmaInternalField[faceCell].component(symmTensor::YZ) + ( grad_sigmaInternalField_YZ[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::ZX) =  sigmaInternalField[faceCell].component(symmTensor::XZ) + ( grad_sigmaInternalField_ZX[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::ZY) =  sigmaInternalField[faceCell].component(symmTensor::YZ) + ( grad_sigmaInternalField_ZY[faceCell] & patchDeltas[i]);
        sigmaExtrapolatedGlobal.component(tensor::ZZ) =  sigmaInternalField[faceCell].component(symmTensor::ZZ) + ( grad_sigmaInternalField_ZZ[faceCell] & patchDeltas[i]);

//        Info << "sigmaInternalField: " << endl << sigmaInternalField[faceCell] << endl;
//        Info << "sigmaExtrapolatedGlobal: " << endl << sigmaExtrapolatedGlobal << endl;

        // sigmaPatchInternal --> sigma at patchCell
        tensor sigmaPatchInternalGlobal;

//        sigmaPatchInternalGlobal.xx() = sigmaPatchInternalField[i].xx();
//        sigmaPatchInternalGlobal.xy() = sigmaPatchInternalField[i].xy();
//        sigmaPatchInternalGlobal.xz() = sigmaPatchInternalField[i].xz();
//        sigmaPatchInternalGlobal.yx() = sigmaPatchInternalField[i].xy();
//        sigmaPatchInternalGlobal.yy() = sigmaPatchInternalField[i].yy();
//        sigmaPatchInternalGlobal.yz() = sigmaPatchInternalField[i].yz();
//        sigmaPatchInternalGlobal.zx() = sigmaPatchInternalField[i].xz();
//        sigmaPatchInternalGlobal.zy() = sigmaPatchInternalField[i].yz();
//        sigmaPatchInternalGlobal.zz() = sigmaPatchInternalField[i].zz();


        // transformed to local coordinates sigmas ...
//        tensor sigmaExtrapolatedLocal = T_g2l & (sigmaExtrapolatedGlobal & T_l2g);
//        tensor sigmaPatchInternalLocal = T_g2l & (sigmaPatchInternalGlobal & T_l2g);


        // VALIDATION MOMENTUM TRANSPORT

//        if(this->patch().name() == "in" || this->patch().name() == "insideCircle"){

//            sigmaExtrapolatedLocal.xy() = 0.1 * (ub & t1);
//            sigmaExtrapolatedLocal.yx() = sigmaExtrapolatedLocal.xy();
//            sigmaExtrapolatedLocal.xx() = sigmaPatchInternalLocal.xx() / (1 + d[i]);
//            sigmaExtrapolatedLocal.yy() = -1 * sigmaExtrapolatedLocal.xx(); // sigma_tt = -sigma_nn

//            //sigmaExtrapolatedLocal.xx() = sigmaPatchInternalLocal.xx() / (1 - d[i]); // partial_n sigma_nn = -1 * sigma_nn
//            //sigmaExtrapolatedLocal.xx() = -sigmaPatchInternalLocal.xx();   // partial_n sigma_nn = 0;
//            //sigmaExtrapolatedLocal.xx() = 0; // sigma_nn = 0;
//            //sigmaExtrapolatedLocal.xy() = 0;
//            //sigmaExtrapolatedLocal.yx() = sigmaExtrapolatedLocal.xy();
//            //sigmaExtrapolatedLocal.yy() = -1 * sigmaExtrapolatedLocal.xx(); // sigma_tt = -sigma_nn
//            //sigmaExtrapolatedLocal.zz() = 0;


//        }
//        else if (this->patch().name() == "out" || this->patch().name() == "outsideCircle"){

//            sigmaExtrapolatedLocal.xx() = sigmaPatchInternalLocal.xx() / (1 + d[i]);
//            //sigmaExtrapolatedLocal.xx() = -sigmaPatchInternalLocal.xx();   // partial_n sigma_nn = 0;
//            sigmaExtrapolatedLocal.yy() = -1 * sigmaExtrapolatedLocal.xx(); // sigma_tt = -sigma_nn

//            //sigmaExtrapolatedLocal.xx() = sigmaPatchInternalLocal.xx() / (1 - d[i]); // partial_n sigma_nn = -1 * sigma_nn
//            //sigmaExtrapolatedLocal.xy() = 0;
//            //sigmaExtrapolatedLocal.yx() = sigmaExtrapolatedLocal.xy();
//            //sigmaExtrapolatedLocal.yy() = -1 * sigmaExtrapolatedLocal.xx(); // sigma_tt = -sigma_nn
//            //sigmaExtrapolatedLocal.zz() = 0;

//        }
//        else
//            Info << "Unknown boundary!!!" << endl;

//        // Set components of transformed sigma that are influenced by boundary condition (simple in local coordinates) ...
//        // ... sigma_nn


//        //        sigmaExtrapolatedLocal.xx() = (     sigmaExtrapolatedLocal.xx() + d[i] * delta1.value() * ( Thetab - ThetaWall[i] ) )
// //                            /
// //                          (     1 - gamma1.value()*d[i] );
//        // ... sigma_nt1
//        //sigmaExtrapolatedLocal.xy() = alpha1.value() * (ub & t1); // + beta1.value() * (sb & t1);
//        //sigmaExtrapolatedLocal.yx() = sigmaExtrapolatedLocal.xy();
//        // ... sigma_nt2
//        //sigmaExtrapolatedLocal.xz() = alpha1.value() * (ub & t2); // + beta1.value() * (sb & t2);
//        //sigmaExtrapolatedLocal.zx() = sigmaExtrapolatedLocal.xz();


        // transform back to global coordinates
//        sigmaExtrapolatedGlobal = T_l2g & (sigmaExtrapolatedLocal & T_g2l);

//        sigmaExtrapolatedGlobal.xz() = 0;
//        sigmaExtrapolatedGlobal.zx() = 0;
//        sigmaExtrapolatedGlobal.yz() = 0;
//        sigmaExtrapolatedGlobal.zy() = 0;
//        sigmaExtrapolatedGlobal.zz() = 0;

        // set boundary-value with the symmetric part
        sigmaBoundary[i] = symm(sigmaExtrapolatedGlobal);
        //Info << "sigmaBoundary[i]: " << endl << sigmaBoundary[i] << endl;

    }

    // Set updated = true
    fixedValueFvPatchSymmTensorField::updateCoeffs();

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
