#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{

template<class vectorType>
void gaussGradDivSymm(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("gradDivSymm(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv2 = 1 / (tvol()*tvol());

    const int XX = 0;
    const int XY = 1;
    const int XZ = 2;
    const int YX = XY;
    const int YY = 3;
    const int YZ = 4;
    const int ZX = XZ;
    const int ZY = YZ;
    const int ZZ = 5;


    if (!updateOnlyRHS){

        // contribution of internal Field ...
        for(int i=0; i<vf.mesh().owner().size(); i++)
        {
            int o = vf.mesh().owner()[i];
            int n = vf.mesh().neighbour()[i];
            scalar w = weights.internalField()[i];
            vector s = vf.mesh().Sf()[i];

            // gradDivSymm(sigma)_ij = 1/V^2 * sum_f 1/2 * ( sigma_ik * S_k * Sj + sigma_jk * S_k * S_i )

            diag[o](dirI + XX, dirJ + XX) += f * w * s.x() * s.x() * volInv2[o];
            diag[o](dirI + XX, dirJ + XY) += f * w * s.y() * s.x() * volInv2[o];
            diag[o](dirI + XX, dirJ + XZ) += f * w * s.z() * s.x() * volInv2[o];

            diag[o](dirI + XY, dirJ + XX) += f * w * 0.5 * s.x() * s.y() * volInv2[o];
            diag[o](dirI + XY, dirJ + XY) += f * w * 0.5 * s.y() * s.y() * volInv2[o];
            diag[o](dirI + XY, dirJ + XZ) += f * w * 0.5 * s.z() * s.y() * volInv2[o];
            diag[o](dirI + XY, dirJ + YX) += f * w * 0.5 * s.x() * s.x() * volInv2[o];
            diag[o](dirI + XY, dirJ + YY) += f * w * 0.5 * s.y() * s.x() * volInv2[o];
            diag[o](dirI + XY, dirJ + YZ) += f * w * 0.5 * s.z() * s.x() * volInv2[o];

            diag[o](dirI + XZ, dirJ + XX) += f * w * 0.5 * s.x() * s.z() * volInv2[o];
            diag[o](dirI + XZ, dirJ + XY) += f * w * 0.5 * s.y() * s.z() * volInv2[o];
            diag[o](dirI + XZ, dirJ + XZ) += f * w * 0.5 * s.z() * s.z() * volInv2[o];
            diag[o](dirI + XZ, dirJ + ZX) += f * w * 0.5 * s.x() * s.x() * volInv2[o];
            diag[o](dirI + XZ, dirJ + ZY) += f * w * 0.5 * s.y() * s.x() * volInv2[o];
            diag[o](dirI + XZ, dirJ + ZZ) += f * w * 0.5 * s.z() * s.x() * volInv2[o];

            diag[o](dirI + YY, dirJ + YX) += f * w * s.x() * s.y() * volInv2[o];
            diag[o](dirI + YY, dirJ + YY) += f * w * s.y() * s.y() * volInv2[o];
            diag[o](dirI + YY, dirJ + YZ) += f * w * s.z() * s.y() * volInv2[o];

            diag[o](dirI + YZ, dirJ + YX) += f * w * 0.5 * s.x() * s.z() * volInv2[o];
            diag[o](dirI + YZ, dirJ + YY) += f * w * 0.5 * s.y() * s.z() * volInv2[o];
            diag[o](dirI + YZ, dirJ + YZ) += f * w * 0.5 * s.z() * s.z() * volInv2[o];
            diag[o](dirI + YZ, dirJ + ZX) += f * w * 0.5 * s.x() * s.y() * volInv2[o];
            diag[o](dirI + YZ, dirJ + ZY) += f * w * 0.5 * s.y() * s.y() * volInv2[o];
            diag[o](dirI + YZ, dirJ + ZZ) += f * w * 0.5 * s.z() * s.y() * volInv2[o];

            diag[o](dirI + ZZ, dirJ + ZX) += f * w * s.x() * s.z() * volInv2[o];
            diag[o](dirI + ZZ, dirJ + ZY) += f * w * s.y() * s.z() * volInv2[o];
            diag[o](dirI + ZZ, dirJ + ZZ) += f * w * s.z() * s.z() * volInv2[o];


            // -1 * -1 = 1  -->  +=
            diag[n](dirI + XX, dirJ + XX) += f * (1-w) * s.x() * s.x() * volInv2[o];
            diag[n](dirI + XX, dirJ + XY) += f * (1-w) * s.y() * s.x() * volInv2[o];
            diag[n](dirI + XX, dirJ + XZ) += f * (1-w) * s.z() * s.x() * volInv2[o];

            diag[n](dirI + XY, dirJ + XX) += f * (1-w) * 0.5 * s.x() * s.y() * volInv2[o];
            diag[n](dirI + XY, dirJ + XY) += f * (1-w) * 0.5 * s.y() * s.y() * volInv2[o];
            diag[n](dirI + XY, dirJ + XZ) += f * (1-w) * 0.5 * s.z() * s.y() * volInv2[o];
            diag[n](dirI + XY, dirJ + YX) += f * (1-w) * 0.5 * s.x() * s.x() * volInv2[o];
            diag[n](dirI + XY, dirJ + YY) += f * (1-w) * 0.5 * s.y() * s.x() * volInv2[o];
            diag[n](dirI + XY, dirJ + YZ) += f * (1-w) * 0.5 * s.z() * s.x() * volInv2[o];

            diag[n](dirI + XZ, dirJ + XX) += f * (1-w) * 0.5 * s.x() * s.z() * volInv2[o];
            diag[n](dirI + XZ, dirJ + XY) += f * (1-w) * 0.5 * s.y() * s.z() * volInv2[o];
            diag[n](dirI + XZ, dirJ + XZ) += f * (1-w) * 0.5 * s.z() * s.z() * volInv2[o];
            diag[n](dirI + XZ, dirJ + ZX) += f * (1-w) * 0.5 * s.x() * s.x() * volInv2[o];
            diag[n](dirI + XZ, dirJ + ZY) += f * (1-w) * 0.5 * s.y() * s.x() * volInv2[o];
            diag[n](dirI + XZ, dirJ + ZZ) += f * (1-w) * 0.5 * s.z() * s.x() * volInv2[o];

            diag[n](dirI + YY, dirJ + YX) += f * (1-w) * s.x() * s.y() * volInv2[o];
            diag[n](dirI + YY, dirJ + YY) += f * (1-w) * s.y() * s.y() * volInv2[o];
            diag[n](dirI + YY, dirJ + YZ) += f * (1-w) * s.z() * s.y() * volInv2[o];

            diag[n](dirI + YZ, dirJ + YX) += f * (1-w) * 0.5 * s.x() * s.z() * volInv2[o];
            diag[n](dirI + YZ, dirJ + YY) += f * (1-w) * 0.5 * s.y() * s.z() * volInv2[o];
            diag[n](dirI + YZ, dirJ + YZ) += f * (1-w) * 0.5 * s.z() * s.z() * volInv2[o];
            diag[n](dirI + YZ, dirJ + ZX) += f * (1-w) * 0.5 * s.x() * s.y() * volInv2[o];
            diag[n](dirI + YZ, dirJ + ZY) += f * (1-w) * 0.5 * s.y() * s.y() * volInv2[o];
            diag[n](dirI + YZ, dirJ + ZZ) += f * (1-w) * 0.5 * s.z() * s.y() * volInv2[o];

            diag[n](dirI + ZZ, dirJ + ZX) += f * (1-w) * s.x() * s.z() * volInv2[o];
            diag[n](dirI + ZZ, dirJ + ZY) += f * (1-w) * s.y() * s.z() * volInv2[o];
            diag[n](dirI + ZZ, dirJ + ZZ) += f * (1-w) * s.z() * s.z() * volInv2[o];


            upper[i](dirI + XX, dirJ + XX) += f * (1-w) * s.x() * s.x() * volInv2[o];
            upper[i](dirI + XX, dirJ + XY) += f * (1-w) * s.y() * s.x() * volInv2[o];
            upper[i](dirI + XX, dirJ + XZ) += f * (1-w) * s.z() * s.x() * volInv2[o];

            upper[i](dirI + XY, dirJ + XX) += f * (1-w) * 0.5 * s.x() * s.y() * volInv2[o];
            upper[i](dirI + XY, dirJ + XY) += f * (1-w) * 0.5 * s.y() * s.y() * volInv2[o];
            upper[i](dirI + XY, dirJ + XZ) += f * (1-w) * 0.5 * s.z() * s.y() * volInv2[o];
            upper[i](dirI + XY, dirJ + YX) += f * (1-w) * 0.5 * s.x() * s.x() * volInv2[o];
            upper[i](dirI + XY, dirJ + YY) += f * (1-w) * 0.5 * s.y() * s.x() * volInv2[o];
            upper[i](dirI + XY, dirJ + YZ) += f * (1-w) * 0.5 * s.z() * s.x() * volInv2[o];

            upper[i](dirI + XZ, dirJ + XX) += f * (1-w) * 0.5 * s.x() * s.z() * volInv2[o];
            upper[i](dirI + XZ, dirJ + XY) += f * (1-w) * 0.5 * s.y() * s.z() * volInv2[o];
            upper[i](dirI + XZ, dirJ + XZ) += f * (1-w) * 0.5 * s.z() * s.z() * volInv2[o];
            upper[i](dirI + XZ, dirJ + ZX) += f * (1-w) * 0.5 * s.x() * s.x() * volInv2[o];
            upper[i](dirI + XZ, dirJ + ZY) += f * (1-w) * 0.5 * s.y() * s.x() * volInv2[o];
            upper[i](dirI + XZ, dirJ + ZZ) += f * (1-w) * 0.5 * s.z() * s.x() * volInv2[o];

            upper[i](dirI + YY, dirJ + YX) += f * (1-w) * s.x() * s.y() * volInv2[o];
            upper[i](dirI + YY, dirJ + YY) += f * (1-w) * s.y() * s.y() * volInv2[o];
            upper[i](dirI + YY, dirJ + YZ) += f * (1-w) * s.z() * s.y() * volInv2[o];

            upper[i](dirI + YZ, dirJ + YX) += f * (1-w) * 0.5 * s.x() * s.z() * volInv2[o];
            upper[i](dirI + YZ, dirJ + YY) += f * (1-w) * 0.5 * s.y() * s.z() * volInv2[o];
            upper[i](dirI + YZ, dirJ + YZ) += f * (1-w) * 0.5 * s.z() * s.z() * volInv2[o];
            upper[i](dirI + YZ, dirJ + ZX) += f * (1-w) * 0.5 * s.x() * s.y() * volInv2[o];
            upper[i](dirI + YZ, dirJ + ZY) += f * (1-w) * 0.5 * s.y() * s.y() * volInv2[o];
            upper[i](dirI + YZ, dirJ + ZZ) += f * (1-w) * 0.5 * s.z() * s.y() * volInv2[o];

            upper[i](dirI + ZZ, dirJ + ZX) += f * (1-w) * s.x() * s.z() * volInv2[o];
            upper[i](dirI + ZZ, dirJ + ZY) += f * (1-w) * s.y() * s.z() * volInv2[o];
            upper[i](dirI + ZZ, dirJ + ZZ) += f * (1-w) * s.z() * s.z() * volInv2[o];


            lower[i](dirI + XX, dirJ + XX) += f * w * s.x() * s.x() * volInv2[n];
            lower[i](dirI + XX, dirJ + XY) += f * w * s.y() * s.x() * volInv2[n];
            lower[i](dirI + XX, dirJ + XZ) += f * w * s.z() * s.x() * volInv2[n];

            lower[i](dirI + XY, dirJ + XX) += f * w * 0.5 * s.x() * s.y() * volInv2[n];
            lower[i](dirI + XY, dirJ + XY) += f * w * 0.5 * s.y() * s.y() * volInv2[n];
            lower[i](dirI + XY, dirJ + XZ) += f * w * 0.5 * s.z() * s.y() * volInv2[n];
            lower[i](dirI + XY, dirJ + YX) += f * w * 0.5 * s.x() * s.x() * volInv2[n];
            lower[i](dirI + XY, dirJ + YY) += f * w * 0.5 * s.y() * s.x() * volInv2[n];
            lower[i](dirI + XY, dirJ + YZ) += f * w * 0.5 * s.z() * s.x() * volInv2[n];

            lower[i](dirI + XZ, dirJ + XX) += f * w * 0.5 * s.x() * s.z() * volInv2[n];
            lower[i](dirI + XZ, dirJ + XY) += f * w * 0.5 * s.y() * s.z() * volInv2[n];
            lower[i](dirI + XZ, dirJ + XZ) += f * w * 0.5 * s.z() * s.z() * volInv2[n];
            lower[i](dirI + XZ, dirJ + ZX) += f * w * 0.5 * s.x() * s.x() * volInv2[n];
            lower[i](dirI + XZ, dirJ + ZY) += f * w * 0.5 * s.y() * s.x() * volInv2[n];
            lower[i](dirI + XZ, dirJ + ZZ) += f * w * 0.5 * s.z() * s.x() * volInv2[n];

            lower[i](dirI + YY, dirJ + YX) += f * w * s.x() * s.y() * volInv2[n];
            lower[i](dirI + YY, dirJ + YY) += f * w * s.y() * s.y() * volInv2[n];
            lower[i](dirI + YY, dirJ + YZ) += f * w * s.z() * s.y() * volInv2[n];

            lower[i](dirI + YZ, dirJ + YX) += f * w * 0.5 * s.x() * s.z() * volInv2[n];
            lower[i](dirI + YZ, dirJ + YY) += f * w * 0.5 * s.y() * s.z() * volInv2[n];
            lower[i](dirI + YZ, dirJ + YZ) += f * w * 0.5 * s.z() * s.z() * volInv2[n];
            lower[i](dirI + YZ, dirJ + ZX) += f * w * 0.5 * s.x() * s.y() * volInv2[n];
            lower[i](dirI + YZ, dirJ + ZY) += f * w * 0.5 * s.y() * s.y() * volInv2[n];
            lower[i](dirI + YZ, dirJ + ZZ) += f * w * 0.5 * s.z() * s.y() * volInv2[n];

            lower[i](dirI + ZZ, dirJ + ZX) += f * w * s.x() * s.z() * volInv2[n];
            lower[i](dirI + ZZ, dirJ + ZY) += f * w * s.y() * s.z() * volInv2[n];
            lower[i](dirI + ZZ, dirJ + ZZ) += f * w * s.z() * s.z() * volInv2[n];

        }

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();

    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<symmTensor>& patchVolField = vf.boundaryField()[patchI];
        const fvPatch& patch = patchVolField.patch();   // reference to patch
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        tmp<symmTensorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<symmTensorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const symmTensorField& ic = tic();     // internal coefficient
        const symmTensorField& bc = tbc();     // boundary coefficient

        forAll(patchVolField, faceI){   // loop all faces

            label c = patch.faceCells()[faceI];     // boundary cell
            const symmTensor ci = ic[faceI];
            const symmTensor cb = bc[faceI];
            const vector s = sn[faceI];

            if(!updateOnlyRHS){

                diag[c](dirI + XX, dirJ + XX) += f * ci.component(XX) * s.x() * s.x() * volInv2[c];
                diag[c](dirI + XX, dirJ + XY) += f * ci.component(XY) * s.y() * s.x() * volInv2[c];
                diag[c](dirI + XX, dirJ + XZ) += f * ci.component(XZ) * s.z() * s.x() * volInv2[c];

                diag[c](dirI + XY, dirJ + XX) += f * ci.component(XX) * 0.5 * s.x() * s.y() * volInv2[c];
                diag[c](dirI + XY, dirJ + XY) += f * ci.component(XY) * 0.5 * s.y() * s.y() * volInv2[c];
                diag[c](dirI + XY, dirJ + XZ) += f * ci.component(XZ) * 0.5 * s.z() * s.y() * volInv2[c];
                diag[c](dirI + XY, dirJ + YX) += f * ci.component(YX) * 0.5 * s.x() * s.x() * volInv2[c];
                diag[c](dirI + XY, dirJ + YY) += f * ci.component(YY) * 0.5 * s.y() * s.x() * volInv2[c];
                diag[c](dirI + XY, dirJ + YZ) += f * ci.component(YZ) * 0.5 * s.z() * s.x() * volInv2[c];

                diag[c](dirI + XZ, dirJ + XX) += f * ci.component(XX) * 0.5 * s.x() * s.z() * volInv2[c];
                diag[c](dirI + XZ, dirJ + XY) += f * ci.component(XY) * 0.5 * s.y() * s.z() * volInv2[c];
                diag[c](dirI + XZ, dirJ + XZ) += f * ci.component(XZ) * 0.5 * s.z() * s.z() * volInv2[c];
                diag[c](dirI + XZ, dirJ + ZX) += f * ci.component(ZX) * 0.5 * s.x() * s.x() * volInv2[c];
                diag[c](dirI + XZ, dirJ + ZY) += f * ci.component(ZY) * 0.5 * s.y() * s.x() * volInv2[c];
                diag[c](dirI + XZ, dirJ + ZZ) += f * ci.component(ZZ) * 0.5 * s.z() * s.x() * volInv2[c];

                diag[c](dirI + YY, dirJ + YX) += f * ci.component(YX) * s.x() * s.y() * volInv2[c];
                diag[c](dirI + YY, dirJ + YY) += f * ci.component(YY) * s.y() * s.y() * volInv2[c];
                diag[c](dirI + YY, dirJ + YZ) += f * ci.component(YZ) * s.z() * s.y() * volInv2[c];

                diag[c](dirI + YZ, dirJ + YX) += f * ci.component(YX) * 0.5 * s.x() * s.z() * volInv2[c];
                diag[c](dirI + YZ, dirJ + YY) += f * ci.component(YY) * 0.5 * s.y() * s.z() * volInv2[c];
                diag[c](dirI + YZ, dirJ + YZ) += f * ci.component(YZ) * 0.5 * s.z() * s.z() * volInv2[c];
                diag[c](dirI + YZ, dirJ + ZX) += f * ci.component(ZX) * 0.5 * s.x() * s.y() * volInv2[c];
                diag[c](dirI + YZ, dirJ + ZY) += f * ci.component(ZY) * 0.5 * s.y() * s.y() * volInv2[c];
                diag[c](dirI + YZ, dirJ + ZZ) += f * ci.component(ZZ) * 0.5 * s.z() * s.y() * volInv2[c];

                diag[c](dirI + ZZ, dirJ + ZX) += f * ci.component(ZX) * s.x() * s.z() * volInv2[c];
                diag[c](dirI + ZZ, dirJ + ZY) += f * ci.component(ZY) * s.y() * s.z() * volInv2[c];
                diag[c](dirI + ZZ, dirJ + ZZ) += f * ci.component(ZZ) * s.z() * s.z() * volInv2[c];

            }


            B[c](dirI + XX) -= f * cb.component(XX) * s.x() * s.x() * volInv2[c];
            B[c](dirI + XX) -= f * cb.component(XY) * s.y() * s.x() * volInv2[c];
            B[c](dirI + XX) -= f * cb.component(XZ) * s.z() * s.x() * volInv2[c];

            B[c](dirI + XY) -= f * 0.5 * cb.component(XX) * s.x() * s.y() * volInv2[c];
            B[c](dirI + XY) -= f * 0.5 * cb.component(XY) * s.y() * s.y() * volInv2[c];
            B[c](dirI + XY) -= f * 0.5 * cb.component(XZ) * s.z() * s.y() * volInv2[c];
            B[c](dirI + XY) -= f * 0.5 * cb.component(YX) * s.x() * s.x() * volInv2[c];
            B[c](dirI + XY) -= f * 0.5 * cb.component(YY) * s.y() * s.x() * volInv2[c];
            B[c](dirI + XY) -= f * 0.5 * cb.component(YZ) * s.z() * s.x() * volInv2[c];

            B[c](dirI + XZ) -= f * 0.5 * cb.component(XX) * s.x() * s.z() * volInv2[c];
            B[c](dirI + XZ) -= f * 0.5 * cb.component(XY) * s.y() * s.z() * volInv2[c];
            B[c](dirI + XZ) -= f * 0.5 * cb.component(XZ) * s.z() * s.z() * volInv2[c];
            B[c](dirI + XZ) -= f * 0.5 * cb.component(ZX) * s.x() * s.x() * volInv2[c];
            B[c](dirI + XZ) -= f * 0.5 * cb.component(ZY) * s.y() * s.x() * volInv2[c];
            B[c](dirI + XZ) -= f * 0.5 * cb.component(ZZ) * s.z() * s.x() * volInv2[c];

            B[c](dirI + YY) -= f * cb.component(YX) * s.x() * s.y() * volInv2[c];
            B[c](dirI + YY) -= f * cb.component(YY) * s.y() * s.y() * volInv2[c];
            B[c](dirI + YY) -= f * cb.component(YZ) * s.z() * s.y() * volInv2[c];

            B[c](dirI + YZ) -= f * 0.5 * cb.component(YX) * s.x() * s.z() * volInv2[c];
            B[c](dirI + YZ) -= f * 0.5 * cb.component(YY) * s.y() * s.z() * volInv2[c];
            B[c](dirI + YZ) -= f * 0.5 * cb.component(YZ) * s.z() * s.z() * volInv2[c];
            B[c](dirI + YZ) -= f * 0.5 * cb.component(ZX) * s.x() * s.y() * volInv2[c];
            B[c](dirI + YZ) -= f * 0.5 * cb.component(ZY) * s.y() * s.y() * volInv2[c];
            B[c](dirI + YZ) -= f * 0.5 * cb.component(ZZ) * s.z() * s.y() * volInv2[c];

            B[c](dirI + ZZ) -= f * cb.component(ZX) * s.x() * s.z() * volInv2[c];
            B[c](dirI + ZZ) -= f * cb.component(ZY) * s.y() * s.z() * volInv2[c];
            B[c](dirI + ZZ) -= f * cb.component(ZZ) * s.z() * s.z() * volInv2[c];

        }

    }

}

}

