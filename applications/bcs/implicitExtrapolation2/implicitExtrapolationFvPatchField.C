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
    fvPatchField<Type>(p, iF)
{

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
    fvPatchField<Type>(ptf, p, iF, mapper)
{

}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict)
{
    fvPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const implicitExtrapolationFvPatchField& zgpf
)
:
    fvPatchField<Type>(zgpf)
{

}


template<class Type>
implicitExtrapolationFvPatchField<Type>::implicitExtrapolationFvPatchField
(
    const implicitExtrapolationFvPatchField& zgpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(zgpf, iF)
{

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
tmp<Field<Type> > implicitExtrapolationFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    Info << "noch implementieren --> gibt erstmal 0 zurück" << endl;
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > implicitExtrapolationFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    Info << "noch implementieren --> gibt erstmal 0 zurück" << endl;
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> >
implicitExtrapolationFvPatchField<Type>::gradientInternalCoeffs() const
{
    Info << "noch implementieren --> gibt erstmal 0 zurück" << endl;
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> >
implicitExtrapolationFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    Info << "noch implementieren --> gibt erstmal 0 zurück" << endl;
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

/*
// JONAS
template<class Type>
void implicitExtrapolationFvPatchField<Type>::setExtraData(const fvPatch& p){


    const fvMesh& m = p.boundaryMesh().mesh();
    const faceList& meshFaces = m.faces();
    const cellList cells = m.cells();
    const labelList& faceCells = p.faceCells();
    const volVectorField& centers = m.C();
    const labelList& neighbours = m.neighbour();
    const labelList& owners = m.owner();

    forAll(faceCells,i){

        const label curCellLabel = faceCells[i];
        const label curPatchFaceLabel = p.patch().start() + i;

        const label secondFace = cells[curCellLabel].opposingFaceLabel(curPatchFaceLabel,meshFaces);
        secondFaces[i] = secondFace;

        if (curCellLabel == m.owner()[secondFace]){

            secondNormalFaceIsOwner[i] = false;
            secondDeltas[i] = mag(centers[curCellLabel] - centers[neighbours[secondFace]]);

        }
        else{

            secondNormalFaceIsOwner[i] = true;
            secondDeltas[i] = mag(centers[curCellLabel] - centers[owners[secondFace]]);

        }

    }

}

*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
