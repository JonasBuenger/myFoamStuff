/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myImplicitVelocityFvPatchVectorField.H"
#include "implicitExtrapolationFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

myImplicitVelocityFvPatchVectorField::myImplicitVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
):
    implicitExtrapolationFvPatchVectorField(p, iF),
    internalCoeffsDiagTensor(p.size()),
    internalCoeffsLowerTensor(p.size()),
    internalCoeffsUpperTensor(p.size())
{
    //Info << "constructor 1" << endl;
}


myImplicitVelocityFvPatchVectorField::myImplicitVelocityFvPatchVectorField
(
    const myImplicitVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
):
    implicitExtrapolationFvPatchVectorField(p, iF),
    internalCoeffsDiagTensor(p.size()),
    internalCoeffsLowerTensor(p.size()),
    internalCoeffsUpperTensor(p.size())
{
    //Info << "constructor 2" << endl;
}

myImplicitVelocityFvPatchVectorField::myImplicitVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
):
    implicitExtrapolationFvPatchVectorField(p, iF, dict),
    internalCoeffsDiagTensor(p.size()),
    internalCoeffsLowerTensor(p.size()),
    internalCoeffsUpperTensor(p.size())
{
    //Info << "constructor 3" << endl;
}


myImplicitVelocityFvPatchVectorField::myImplicitVelocityFvPatchVectorField
(
    const myImplicitVelocityFvPatchVectorField& zgpf
):
    implicitExtrapolationFvPatchVectorField(zgpf),
    internalCoeffsDiagTensor(zgpf.size()),
    internalCoeffsLowerTensor(zgpf.size()),
    internalCoeffsUpperTensor(zgpf.size())
{
    //Info << "constructor 4" << endl;
}


myImplicitVelocityFvPatchVectorField::myImplicitVelocityFvPatchVectorField
(
    const myImplicitVelocityFvPatchVectorField& zgpf,
    const DimensionedField<vector, volMesh>& iF
):
    implicitExtrapolationFvPatchVectorField(zgpf, iF),
    internalCoeffsDiagTensor(zgpf.size()),
    internalCoeffsLowerTensor(zgpf.size()),
    internalCoeffsUpperTensor(zgpf.size())
{
    //Info << "constructor 5" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void myImplicitVelocityFvPatchVectorField::updateCoeffs(){

    implicitExtrapolationFvPatchVectorField::updateCoeffs();

    tmp<vectorField> tn = this->patch().nf();
    const vectorField& n = tn();

    tensorField& icd = internalCoeffsDiagTensor;
    tensorField& icu = internalCoeffsUpperTensor;
    tensorField& icl = internalCoeffsLowerTensor;

    forAll(icd,i){

        const tensor nn = n[i] * n[i];
        const vector ic = internalCoeffsDiag[i];

        icd[i].xx() = ic.x() - ic.x() * nn.xx();
        icd[i].xy() =        - ic.y() * nn.xy();
        icd[i].xz() =        - ic.z() * nn.xz();
        icd[i].yx() =        - ic.x() * nn.yx();
        icd[i].yy() = ic.y() - ic.y() * nn.yy();
        icd[i].yz() =        - ic.z() * nn.yz();
        icd[i].zx() =        - ic.x() * nn.zx();
        icd[i].zy() =        - ic.y() * nn.zy();
        icd[i].zz() = ic.z() - ic.z() * nn.zz();

    }

    forAll(icu,i){

        const tensor nn = n[i] * n[i];
        const vector ic = internalCoeffsUpper[i];

        icu[i].xx() = ic.x() - ic.x() * nn.xx();
        icu[i].xy() =        - ic.y() * nn.xy();
        icu[i].xz() =        - ic.z() * nn.xz();
        icu[i].yx() =        - ic.x() * nn.yx();
        icu[i].yy() = ic.y() - ic.y() * nn.yy();
        icu[i].yz() =        - ic.z() * nn.yz();
        icu[i].zx() =        - ic.x() * nn.zx();
        icu[i].zy() =        - ic.y() * nn.zy();
        icu[i].zz() = ic.z() - ic.z() * nn.zz();

    }

    forAll(icl,i){

        const tensor nn = n[i] * n[i];
        const vector ic = internalCoeffsLower[i];

        icl[i].xx() = ic.x() - ic.x() * nn.xx();
        icl[i].xy() =        - ic.y() * nn.xy();
        icl[i].xz() =        - ic.z() * nn.xz();
        icl[i].yx() =        - ic.x() * nn.yx();
        icl[i].yy() = ic.y() - ic.y() * nn.yy();
        icl[i].yz() =        - ic.z() * nn.yz();
        icl[i].zx() =        - ic.x() * nn.zx();
        icl[i].zy() =        - ic.y() * nn.zy();
        icl[i].zz() = ic.z() - ic.z() * nn.zz();

    }

    extrapolateUsingOldValues();

    fvPatchField<vector>::updateCoeffs();

}


void myImplicitVelocityFvPatchVectorField::extrapolateUsingOldValues()
{
    // Werte aus altem Zeitschritt setzen -> werden für die Auswertung von Sigma am Rand verwendet.

    tmp<vectorField> tboundaryField = *this;
    vectorField& boundaryField = tboundaryField();
    const labelList& secondFace = secondFaces();
    const labelList& fc = patch().faceCells();
    tmp<vectorField> tu = this->internalField();
    const vectorField u = tu();
    const labelList& o = this->patch().boundaryMesh().mesh().owner();
    const labelList& n = this->patch().boundaryMesh().mesh().neighbour();

    const tensorField& icds = internalCoeffsDiagTensor;
    const tensorField& icus = internalCoeffsUpperTensor;
    const tensorField& icls = internalCoeffsLowerTensor;

    forAll(boundaryField, i){

        vector b(0,0,0);
        const tensor icd = icds[i];
        const tensor icl = icls[i];
        const tensor icu = icus[i];

        b += icd & u[fc[i]];
        b += icu & u[o[secondFace[i]]];
        b += icl & u[n[secondFace[i]]];

        boundaryField[i] = b;
        //Info << "b: " << b << endl;

    }

}

tmp<Field<tensor> > myImplicitVelocityFvPatchVectorField::valueInternalCoeffsTensor
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsDiagTensor;
}

// JONAS
tmp<Field<tensor> > myImplicitVelocityFvPatchVectorField::valueInternalCoeffsUpperTensor
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsUpperTensor;
}

// JONAS
tmp<Field<tensor> > myImplicitVelocityFvPatchVectorField::valueInternalCoeffsLowerTensor
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsLowerTensor;
}

tmp<Field<vector> > myImplicitVelocityFvPatchVectorField::gradientInternalCoeffs() const
{
    return -pTraits<vector>::one*this->patch().deltaCoeffs();
}

tmp<Field<vector> > myImplicitVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myImplicitVelocityFvPatchVectorField
    )
}

// ************************************************************************* //


// ************************************************************************* //
