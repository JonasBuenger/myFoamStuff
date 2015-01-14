#include"fvCFD.H"
#include"snGradScheme.H"

namespace myFvm{

template<class vectorType>
void gaussLaplacian(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    Foam::fvScalarMatrix m = fvm::laplacian(vf);
    m.addBoundaryDiag(m.diag(),0);
    m.addBoundarySource(m.source(),false);

    if (!updateOnlyRHS){

        forAll(m.diag(), i){
            diag[i](dirI, dirJ) += f*m.diag()[i];
        }

        forAll(m.upper(), i){
            upper[i](dirI, dirJ) += f*m.upper()[i];
        }

        forAll(m.lower(), i){
            lower[i](dirI, dirJ) += f*m.lower()[i];
        }

    }

    forAll(m.source(), i){
        B[i](dirI) += f*m.source()[i];
    }

}
/*
template<class vectorType>
void gaussLaplacian(volVectorField& vf, scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){


    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    Foam::fvVectorMatrix m = fvm::laplacian(vf);
    m.addBoundaryDiag(m.diag(),0);
    m.addBoundarySource(m.source(),false);

    //f = f*2;

    if (!updateOnlyRHS){

        forAll(m.diag(), i){
            diag[i](dirI+0, dirJ+0) += f*m.diag()[i];
            diag[i](dirI+1, dirJ+1) += f*m.diag()[i];
            diag[i](dirI+2, dirJ+2) += f*m.diag()[i];
        }

        forAll(m.upper(), i){
            upper[i](dirI+0, dirJ+0) += f*m.upper()[i];
            upper[i](dirI+1, dirJ+1) += f*m.upper()[i];
            upper[i](dirI+2, dirJ+2) += f*m.upper()[i];
        }

        forAll(m.lower(), i){
            lower[i](dirI+0, dirJ+0) += f*m.lower()[i];
            lower[i](dirI+1, dirJ+1) += f*m.lower()[i];
            lower[i](dirI+2, dirJ+2) += f*m.lower()[i];
        }

    }

    forAll(m.source(), i){
        B[i](dirI+0) += f*m.source()[i].x();
        B[i](dirI+1) += f*m.source()[i].y();
        B[i](dirI+2) += f*m.source()[i].z();
    }


    tmp<volVectorField> tlaplacianS = fvc::laplacian(vf);
    volVectorField& laplacianS = tlaplacianS();
    //blockMatrixTools::blockAdd(dirI+0, 0.5 * laplacianS.internalField().component(0), B);
    //blockMatrixTools::blockAdd(dirI+1, 0.5 * laplacianS.internalField().component(1), B);
    //blockMatrixTools::blockAdd(dirI+2, 0.5 * laplacianS.internalField().component(2), B);

    vectorField& internalField = laplacianS.internalField();

    forAll(laplacianS.internalField(), i){
        //B[i](dirI+0) += -0.5*f*internalField[i].component(0);
        //B[i](dirI+1) += -0.5*f*internalField[i].component(1);
        //B[i](dirI+2) += -0.5*f*internalField[i].component(2);
    }


}
*/


template<class vectorType>
void gaussLaplacian(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const fvMesh& mesh = vf.mesh();

    surfaceScalarField gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf
    (
        gamma*mesh.magSf()
    );

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    // Laplacian Uncorrected
    tmp<Foam::fv::snGradScheme<vector> > tsnGradScheme =
            Foam::fv::snGradScheme<vector>::New(
                vf.mesh(),
                vf.mesh().schemesDict().snGradScheme("snGrad(" + vf.name() + ")")
            );

    tmp<surfaceScalarField> tdeltaCoeffs = tsnGradScheme().deltaCoeffs(vf);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    const labelList& o = vf.mesh().owner();
    const labelList& n = vf.mesh().neighbour();

    const scalarField gammaMagSfInternalField = gammaMagSf.internalField();
    const scalarField deltaCoeffsInternalField = deltaCoeffs.internalField();

    if (!updateOnlyRHS){

        // Internal Field
        forAll(upper, faceI){

            const label owner = o[faceI];
            const label neighbour = n[faceI];

            diag[owner](dirI+0, dirJ+0) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+1, dirJ+1) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+2, dirJ+2) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];

            diag[neighbour](dirI+0, dirJ+0) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+1, dirJ+1) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+2, dirJ+2) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];

            upper[faceI](dirI+0, dirJ+0) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+1, dirJ+1) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+2, dirJ+2) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];

            lower[faceI](dirI+0, dirJ+0) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+1, dirJ+1) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+2, dirJ+2) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];

        }
    }

    // Boundary Field
    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tgradIC = patchVolField.gradientInternalCoeffs();
        const vectorField& gradIC = tgradIC();
        tmp<vectorField> tgradBC = patchVolField.gradientBoundaryCoeffs();
        const vectorField& gradBC = tgradBC();

        const fvsPatchScalarField& patchGamma = gammaMagSf.boundaryField()[patchI];

        forAll(patchVolField, faceI){

            const label c = patch.faceCells()[faceI];

            if(!updateOnlyRHS){

                diag[c](dirI+0, dirJ+0) += f*gradIC[faceI].component(vector::X) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+1, dirJ+1) += f*gradIC[faceI].component(vector::Y) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+2, dirJ+2) += f*gradIC[faceI].component(vector::Z) * patchGamma[faceI] * volInv[c];

            }

            B[c](dirI+0) -= f*gradBC[faceI].component(vector::X) * patchGamma[faceI] * volInv[c];
            B[c](dirI+1) -= f*gradBC[faceI].component(vector::Y) * patchGamma[faceI] * volInv[c];
            B[c](dirI+2) -= f*gradBC[faceI].component(vector::Z) * patchGamma[faceI] * volInv[c];

        }

    }

    // Correction
    if (tsnGradScheme().corrected())
        {
            if (mesh.schemesDict().fluxRequired(vf.name()))
            {
                const surfaceVectorField* faceFluxCorrectionPtr = new
                surfaceVectorField
                (
                    gammaMagSf*tsnGradScheme().correction(vf)
                );

                const scalarField V = mesh.V();
                const vectorField divFaceFluxCorr = fvc::div
                        (
                            *faceFluxCorrectionPtr
                        )().internalField();


                forAll(B,c){

                    B[c](dirI + 0) -= V[c] * divFaceFluxCorr[c].component(vector::X) * volInv[c];
                    B[c](dirI + 1) -= V[c] * divFaceFluxCorr[c].component(vector::Y) * volInv[c];
                    B[c](dirI + 2) -= V[c] * divFaceFluxCorr[c].component(vector::Z) * volInv[c];

                }

            }
            else
            {

                const scalarField V = mesh.V();
                const vectorField divFaceFluxCorr = fvc::div
                        (
                            gammaMagSf*tsnGradScheme().correction(vf)
                        )().internalField();

                forAll(B,c){

                    B[c](dirI + 0) -= V[c] * divFaceFluxCorr[c].component(vector::X);
                    B[c](dirI + 1) -= V[c] * divFaceFluxCorr[c].component(vector::Y);
                    B[c](dirI + 2) -= V[c] * divFaceFluxCorr[c].component(vector::Z);

                }
            }
        }


}


template<class vectorType>
void gaussLaplacian(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const fvMesh& mesh = vf.mesh();

    surfaceScalarField gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf
    (
        gamma*mesh.magSf()
    );

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    // Laplacian Uncorrected
    tmp<Foam::fv::snGradScheme<symmTensor> > tsnGradScheme =
            Foam::fv::snGradScheme<symmTensor>::New(
                vf.mesh(),
                vf.mesh().schemesDict().snGradScheme("snGrad(" + vf.name() + ")")
            );

    tmp<surfaceScalarField> tdeltaCoeffs = tsnGradScheme().deltaCoeffs(vf);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    const labelList& o = vf.mesh().owner();
    const labelList& n = vf.mesh().neighbour();

    const scalarField gammaMagSfInternalField = gammaMagSf.internalField();
    const scalarField deltaCoeffsInternalField = deltaCoeffs.internalField();

    if(!updateOnlyRHS){

        // Internal Field
        forAll(upper, faceI){

            const label owner = o[faceI];
            const label neighbour = n[faceI];

            diag[owner](dirI+0, dirJ+0) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+1, dirJ+1) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+2, dirJ+2) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+3, dirJ+3) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+4, dirJ+4) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            diag[owner](dirI+5, dirJ+5) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];

            diag[neighbour](dirI+0, dirJ+0) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+1, dirJ+1) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+2, dirJ+2) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+3, dirJ+3) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+4, dirJ+4) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            diag[neighbour](dirI+5, dirJ+5) -= f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];

            upper[faceI](dirI+0, dirJ+0) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+1, dirJ+1) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+2, dirJ+2) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+3, dirJ+3) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+4, dirJ+4) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];
            upper[faceI](dirI+5, dirJ+5) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[owner];

            lower[faceI](dirI+0, dirJ+0) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+1, dirJ+1) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+2, dirJ+2) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+3, dirJ+3) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+4, dirJ+4) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];
            lower[faceI](dirI+5, dirJ+5) += f*deltaCoeffsInternalField[faceI] * gammaMagSfInternalField[faceI] * volInv[neighbour];

        }

    }

    // Boundary Field
    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<symmTensor>& patchVolField = vf.boundaryField()[patchI];
        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<symmTensorField> tgradIC = patchVolField.gradientInternalCoeffs();
        const symmTensorField& gradIC = tgradIC();
        tmp<symmTensorField> tgradBC = patchVolField.gradientBoundaryCoeffs();
        const symmTensorField& gradBC = tgradBC();

        const fvsPatchScalarField& patchGamma = gammaMagSf.boundaryField()[patchI];

        forAll(patchVolField, faceI){

            const label c = patch.faceCells()[faceI];

            if (!updateOnlyRHS){

                diag[c](dirI+0, dirJ+0) += f*gradIC[faceI].component(symmTensor::XX) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+1, dirJ+1) += f*gradIC[faceI].component(symmTensor::XY) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+2, dirJ+2) += f*gradIC[faceI].component(symmTensor::XZ) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+3, dirJ+3) += f*gradIC[faceI].component(symmTensor::YY) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+4, dirJ+4) += f*gradIC[faceI].component(symmTensor::YZ) * patchGamma[faceI] * volInv[c];
                diag[c](dirI+5, dirJ+5) += f*gradIC[faceI].component(symmTensor::ZZ) * patchGamma[faceI] * volInv[c];

            }

            B[c](dirI+0) -= f*gradBC[faceI].component(symmTensor::XX) * patchGamma[faceI] * volInv[c];
            B[c](dirI+1) -= f*gradBC[faceI].component(symmTensor::XY) * patchGamma[faceI] * volInv[c];
            B[c](dirI+2) -= f*gradBC[faceI].component(symmTensor::XZ) * patchGamma[faceI] * volInv[c];
            B[c](dirI+3) -= f*gradBC[faceI].component(symmTensor::YY) * patchGamma[faceI] * volInv[c];
            B[c](dirI+4) -= f*gradBC[faceI].component(symmTensor::YZ) * patchGamma[faceI] * volInv[c];
            B[c](dirI+5) -= f*gradBC[faceI].component(symmTensor::ZZ) * patchGamma[faceI] * volInv[c];

        }

    }

    // Correction
    if (tsnGradScheme().corrected())
        {
            if (mesh.schemesDict().fluxRequired(vf.name()))
            {
                const surfaceSymmTensorField* faceFluxCorrectionPtr = new
                surfaceSymmTensorField
                (
                    gammaMagSf*tsnGradScheme().correction(vf)
                );

                const scalarField V = mesh.V();
                const symmTensorField divFaceFluxCorr = fvc::div
                        (
                            *faceFluxCorrectionPtr
                        )().internalField();


                forAll(B,c){

                    B[c](dirI + 0) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XX) * volInv[c];
                    B[c](dirI + 1) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XY) * volInv[c];
                    B[c](dirI + 2) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XZ) * volInv[c];
                    B[c](dirI + 3) -= V[c] * divFaceFluxCorr[c].component(symmTensor::YY) * volInv[c];
                    B[c](dirI + 4) -= V[c] * divFaceFluxCorr[c].component(symmTensor::YZ) * volInv[c];
                    B[c](dirI + 5) -= V[c] * divFaceFluxCorr[c].component(symmTensor::ZZ) * volInv[c];

                }

            }
            else
            {

                const scalarField V = mesh.V();
                const symmTensorField divFaceFluxCorr = fvc::div
                        (
                            gammaMagSf*tsnGradScheme().correction(vf)
                        )().internalField();

                forAll(B,c){

                    B[c](dirI + 0) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XX);
                    B[c](dirI + 1) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XY);
                    B[c](dirI + 2) -= V[c] * divFaceFluxCorr[c].component(symmTensor::XZ);
                    B[c](dirI + 3) -= V[c] * divFaceFluxCorr[c].component(symmTensor::YY);
                    B[c](dirI + 4) -= V[c] * divFaceFluxCorr[c].component(symmTensor::YZ);
                    B[c](dirI + 5) -= V[c] * divFaceFluxCorr[c].component(symmTensor::ZZ);

                }
            }
        }


}

}
