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

#include "implicitExtrapolationFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    secondNormalCellIsOwner(p.size()),
    secondDeltas(p.size()),
    secondFacesIDs(p.size()),
    internalCoeffsUpper(p.size()),
    internalCoeffsLower(p.size()),
    internalCoeffsDiag(p.size())
{
    //Info << "constructor 1" << endl;
    setExtraData(p);
}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const implicitExtrapolationFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    secondNormalCellIsOwner(p.size()),
    secondDeltas(p.size()),
    secondFacesIDs(p.size()),
    internalCoeffsUpper(p.size()),
    internalCoeffsLower(p.size()),
    internalCoeffsDiag(p.size())
{
    //Info << "constructor 2" << endl;
    setExtraData(p);
}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    secondNormalCellIsOwner(p.size()),
    secondDeltas(p.size()),
    secondFacesIDs(p.size()),
    internalCoeffsUpper(p.size()),
    internalCoeffsLower(p.size()),
    internalCoeffsDiag(p.size())
{
    //Info << "constructor 3" << endl;
    setExtraData(p);
    fvPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const implicitExtrapolationFvPatchField& zgpf
)
:
    fvPatchField<Type>(zgpf),
    secondNormalCellIsOwner(zgpf.size()),
    secondDeltas(zgpf.size()),
    secondFacesIDs(zgpf.size()),
    internalCoeffsUpper(zgpf.size()),
    internalCoeffsLower(zgpf.size()),
    internalCoeffsDiag(zgpf.size())
{
    //Info << "constructor 4" << endl;
    setExtraData(zgpf.patch());
}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const implicitExtrapolationFvPatchField& zgpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(zgpf, iF),
    secondNormalCellIsOwner(zgpf.size()),
    secondDeltas(zgpf.size()),
    secondFacesIDs(zgpf.size()),
    internalCoeffsUpper(zgpf.size()),
    internalCoeffsLower(zgpf.size()),
    internalCoeffsDiag(zgpf.size())
{
    //Info << "constructor 5" << endl;
    setExtraData(zgpf.patch());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void implicitExtrapolationFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::operator==(this->patchInternalField());
    fvPatchField<Type>::evaluate();
}

template<class Type>
const List<label>& implicitExtrapolationFvPatchField<Type>::secondFaces() const{

    return secondFacesIDs;

}

template<class Type>
void implicitExtrapolationFvPatchField<Type>::updateCoeffs(){


    Field<Type>& boundaryField = *this;                         // will be set from old time-step
    const Field<Type>& internalField = this->internalField();
    const fvMesh& m = this->patch().boundaryMesh().mesh();
    const scalarField& d1 = this->patch().deltaCoeffs();
    const scalarField& d2 = secondDeltas;
    const scalarField d1_dividedBy_d2 = (1/d1)/d2;
    const Field<Type> oneField = Field<Type>(this->size(), pTraits<Type>::one);
    const labelList& faceCells = this->patch().faceCells();

    internalCoeffsDiag = (1 + d1_dividedBy_d2 ) * oneField;

    forAll(boundaryField, i){

        if(secondNormalCellIsOwner[i]){
            internalCoeffsLower[i] = - d1_dividedBy_d2[i] * oneField[i];
            internalCoeffsUpper[i] = 0 * oneField[i];
            const label boundaryCell = faceCells[i];
            const label secondCell = m.owner()[secondFaces()[i]];
            boundaryField[i] = cmptMultiply(internalCoeffsDiag[i],internalField[boundaryCell])
                             + cmptMultiply(internalCoeffsLower[i],internalField[secondCell]);
        }
        else{
            internalCoeffsLower[i] = 0 * oneField[i];
            internalCoeffsUpper[i] = - d1_dividedBy_d2[i] * oneField[i];
            const label boundaryCell = faceCells[i];
            const label secondCell = m.neighbour()[secondFaces()[i]];
            boundaryField[i] = cmptMultiply(internalCoeffsDiag[i],internalField[boundaryCell])
                             + cmptMultiply(internalCoeffsUpper[i],internalField[secondCell]);
        }

    }

    fvPatchField<Type>::updateCoeffs();

}


template<class Type>
tmp<Field<Type> > implicitExtrapolationFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsDiag;
}

// JONAS
template<class Type>
tmp<Field<Type> > implicitExtrapolationFvPatchField<Type>::valueInternalCoeffsUpper
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsUpper;
}

// JONAS
template<class Type>
tmp<Field<Type> > implicitExtrapolationFvPatchField<Type>::valueInternalCoeffsLower
(
    const tmp<scalarField>&
) const
{
    return internalCoeffsLower;
}

template<class Type>
void implicitExtrapolationFvPatchField<Type>::setExtraData(const fvPatch& p){

    const fvMesh& m = p.boundaryMesh().mesh();
    const faceList& meshFaces = m.faces();
    const cellList cells = m.cells();
    const labelList& faceCells = p.faceCells();
    const volVectorField& centers = m.C();

    forAll(faceCells,i){

        const label curCellLabel = faceCells[i];
        const label curPatchFaceLabel = p.patch().start() + i;

        const label secondFace = cells[curCellLabel].opposingFaceLabel(curPatchFaceLabel,meshFaces);
        secondFacesIDs[i] = secondFace;
        const label ownerSecondFace = m.owner()[secondFace];
        const label neighbourSecondFace = m.neighbour()[secondFace];

        if (curCellLabel == ownerSecondFace){
            secondNormalCellIsOwner[i] = 0;
        }
        else{
            secondNormalCellIsOwner[i] = 1;
        }

        secondDeltas[i] = mag(centers[ownerSecondFace] - centers[neighbourSecondFace]);

    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
