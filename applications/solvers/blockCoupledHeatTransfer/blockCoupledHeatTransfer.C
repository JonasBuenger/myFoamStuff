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

Application
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fieldTypes.H"
#include "Time.H"
#include "fvMesh.H"

#include "blockLduSolvers.H"
#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockVectorNMatrices.H"
#include "blockMatrixTools.H"

#include "myImplicitSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating heat transport\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        int nNonOrthCorr = 0;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

            /**************************************************************
             * CONSTRUCT BLOCK-MATRIX M AND BLOCK-SOURCE-VECTOR S
             *************************************************************/

            // ... block matrix
            BlockLduMatrix<vector4> M(mesh);
            Field<tensor4>& d = M.diag().asSquare();
            Field<tensor4>& u = M.upper().asSquare();
            Field<tensor4>& l = M.lower().asSquare();

            // ... block source
            Field<vector4> S(mesh.nCells(), vector4::zero);


            /**************************************************************
             * SET MATRIX COEFFICIENTS
             *************************************************************/

            // 1)TRANSPORT-EQUATION

            // 1.1) heat flux s

            s.boundaryField().updateCoeffs();
            fvVectorMatrix sEqn
            (
                fvm::Sp(corrDim, s)
              - fvm::laplacian(nu,s)
            );
            tmp<scalarField> tdiag = sEqn.D();   // Function D() includes the boundary diagonal contribution
            scalarField& sDiag = tdiag();
            scalarField& sUpper = sEqn.upper();
            scalarField& sLower = sEqn.lower();
            vectorField& sSource = sEqn.source();
            sEqn.addBoundarySource(sSource, false);   // add contribution of boundary to source


            // 1.2) temperatur T
            //
            // ... there is no implicit discretisation for the gradient operator ...
            // ... we will have to discretize implicitly ourselfes
            //
            // \int_V grad(T) dV  = \int_S T dS        (Gauss ?)
            //                   ~= \sum_f T_f S_f
            //
            // Wir benötigen die Temperatur auf den Faces und Rändern. Openfoam
            // interpoliert auf eine Face f
            //
            //      T_f = w*T_O + (1-w)*T_N
            //
            // Die Gewichte w erhält man über ein InterpolationsObjekt ...

            tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
                    surfaceInterpolationScheme<scalar>::New(
                            Theta.mesh(),
                            mesh.schemesDict().interpolationScheme("grad(Theta)")
                    )
            );
            tmp<surfaceScalarField> tweights = tinterpScheme_().weights(Theta);
            const surfaceScalarField& weights = tweights();


            // Create 4 fields (diag, lower, upper, source) for the coefficients of an implicit gradient operator
            tmp<vectorField> tTDiag = tmp<vectorField>(new vectorField(sDiag.size(), pTraits<vector>::zero));
            vectorField& TDiag = tTDiag();
            tmp<vectorField> tTUpper = tmp<vectorField>(new vectorField(sUpper.size(), pTraits<vector>::zero));
            vectorField& TUpper = tTUpper();
            tmp<vectorField> tTLower = tmp<vectorField>(new vectorField(sLower.size(), pTraits<vector>::zero));
            vectorField& TLower = tTLower();
            tmp<vectorField> tTSource = tmp<vectorField>(new vectorField(sSource.size(), pTraits<vector>::zero));
            vectorField& TSource = tTSource();

            // contribution of internal Field ...
            for(int i=0; i<mesh.owner().size(); i++){
                int o = mesh.owner()[i];
                int n = mesh.neighbour()[i];
                scalar w = weights.internalField()[i];
                vector s = mesh.Sf()[i];

                TDiag[o] +=  s*w;
                TDiag[n] -=  s*(1-w);
                TLower[i] = -s*w;
                TUpper[i] =  s*(1-w);

            }

            // contribution of boundary Field ...
            Theta.boundaryField().updateCoeffs();
            forAll(Theta.boundaryField(), patchI){

                // finite volume patch field of Theta of present patch
                fvPatchField<scalar>& pT = Theta.boundaryField()[patchI];
                // weights, that the boundary values are multiplied --> here: [1 1 1 ... 1]
                const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

                // Diagonal- and Source-Contribution of Theta at boundary
                //
                // T_patchFace = ic*T_N + bc
                //
                // ic und bc werden von der Randbedingung vorgegeben und werden über
                // ... valueInternalCoeffs und valueBoundaryCoeffs abgefragt
                //
                tmp<Field<scalar> > tic = pT.valueInternalCoeffs(pw);
                tmp<Field<scalar> > tbc = pT.valueBoundaryCoeffs(pw);
                const Field<scalar>& ic = tic();            // internal coefficient
                const Field<scalar>& bc = tbc();            // boundary coefficient

                // reference to patch
                const fvPatch& patch = pT.patch();

                // reference to patch normals
                tmp<Field<vector> > tsn = patch.Sf();
                Field<vector> sn = tsn();                   // patch normal vectors

                forAll(pT, faceI){
                    label c = patch.faceCells()[faceI];     // boundary cell

                    TDiag[c]   += ic[faceI]*sn[faceI];
                    TSource[c] -= bc[faceI]*sn[faceI];
                }

            }


            // set blockMatrix coefficients for equations 0,1 und 2 ...

            // ...diagonal elements
            forAll(d,i){
                d[i](0,0) = sDiag[i];
                d[i](1,1) = sDiag[i];
                d[i](2,2) = sDiag[i];

                d[i](0,3) = TDiag[i].x();
                d[i](1,3) = TDiag[i].y();
                d[i](2,3) = TDiag[i].z();
            }

            // ...lower elements
            forAll(l,i){
                l[i](0,0) = sLower[i];
                l[i](1,1) = sLower[i];
                l[i](2,2) = sLower[i];

                l[i](0,3) = TLower[i].x();
                l[i](1,3) = TLower[i].y();
                l[i](2,3) = TLower[i].z();
            }

            // ...upper elements
            forAll(u,i){
                u[i](0,0) = sUpper[i];
                u[i](1,1) = sUpper[i];
                u[i](2,2) = sUpper[i];

                u[i](0,3) = TUpper[i].x();
                u[i](1,3) = TUpper[i].y();
                u[i](2,3) = TUpper[i].z();
            }

            // ...source elements
            forAll(S,i){
                S[i](0) = sSource[i].x() + TSource[i].x();
                S[i](1) = sSource[i].y() + TSource[i].y();
                S[i](2) = sSource[i].z() + TSource[i].z();
            }



            // 2) continuity-equation
            //
            // \int_V nabla * s dV  = \int_S s * dS
            //                     ~= \sum_F s_F * S_F
            //
            // ... Temperature-Heat-Coupling via Rhie-Chow interpolation <<???>>
            //
            // \sum_F -D_F * grad(T_F) * S_F  +  \sum_F s_F * S_F     = \sum_F -D_F * grad(T_F) * S_F
            //
            // ------- term 1 (impl.) ------  +  -- term 2 (impl.) -- = ------ term 3 (expl.) -------
            //
            // Di_F = [ #1 0 0 ; 0 #2 0 ; 0 0 #3 ]   mit #dir = Volume_i / diagCoeffOfsEqn for dir


            // ... term 1 and 3

            tmp<volScalarField> tA = sEqn.A();      // .A() returns central coefficients (Wann wir mit dem Volumen multipliziert?)
            volScalarField& A = tA();

            tmp<volVectorField> texpGradTheta = fvc::grad(Theta);    // explicit gradient
            volVectorField& expGradTheta = texpGradTheta();
            tmp<volVectorField> texpGradTheta2 = expGradTheta/A;              // gradient divided by central coefficients
            volVectorField expGradTheta2 = texpGradTheta2();
            tmp<fvScalarMatrix> MEqnTpart1
            (
                - fvm::laplacian(1/A,Theta)
                  ==
                - fvc::div(expGradTheta2)
            );
            // ... unterteilen in term 1 und 3 --> anschaulicher

            scalarField& TMdiag = MEqnTpart1().diag();
            scalarField& TMupper = MEqnTpart1().upper();
            scalarField& TMlower = MEqnTpart1().lower();

            // Add boundary contribution to TMdiag
            MEqnTpart1().addBoundaryDiag(TMdiag, 0);

            // Add boundary contribution to source
            scalarField& TMsource = MEqnTpart1().source();
            MEqnTpart1().addBoundarySource(TMsource, false);


            // term 2
            //
            // as again no implicit divergence discretization is provided by Foam --> one has to discretize

            // interpolationScheme
            tmp<surfaceInterpolationScheme<scalar> > stinterpScheme_(
                surfaceInterpolationScheme<scalar>::New(s.mesh(),mesh.schemesDict().interpolationScheme("div(s)(implicit)"))
            );

            // set up fields to store coefficients
            tmp<vectorField> tsDiag2 = tmp<vectorField>(new vectorField(sDiag.size(), pTraits<vector>::zero));
            vectorField& sDiag2 = tsDiag2();
            tmp<vectorField> tsUpper2 = tmp<vectorField>(new vectorField(sUpper.size(), pTraits<vector>::zero));
            vectorField& sUpper2 = tsUpper2();
            tmp<vectorField> tsLower2 = tmp<vectorField>(new vectorField(sLower.size(), pTraits<vector>::zero));
            vectorField& sLower2 = tsLower2();
            tmp<vectorField> tsSource2 = tmp<vectorField>(new vectorField(sSource.component(0)().size(), pTraits<vector>::zero));
            vectorField& sSource2 = tsSource2();

            // interpolation weights
            tmp<surfaceScalarField> tsWeights = stinterpScheme_().weights(mag(s));
            const surfaceScalarField& sWeights = tsWeights();

            // coeffs of internal field
            for(int i=0; i<mesh.owner().size(); i++){
                int o = mesh.owner()[i];
                int n = mesh.neighbour()[i];
                scalar w = sWeights.internalField()[i];
                vector s = mesh.Sf()[i];

                sDiag2[o]  +=  s*w;
                sDiag2[n]  -=  s*(1-w);
                sLower2[i]  = -s*w;
                sUpper2[i]  =  s*(1-w);
            }

            // coeffs and source of boundary
            s.boundaryField().updateCoeffs();
            forAll(s.boundaryField(), patchI){

                // present fvPatchField
                fvPatchField<vector>& ps = s.boundaryField()[patchI];

                // boundary weights
                const fvsPatchScalarField pw = sWeights.boundaryField()[patchI];

                // contribution from boundary to diagonal and source
                tmp<Field<vector> > tic = ps.valueInternalCoeffs(pw);
                Field<vector> ic = tic();
                tmp<Field<vector> > tbc = ps.valueBoundaryCoeffs(pw);
                Field<vector> bc = tbc();

                // fvPatch
                const fvPatch& patch = ps.patch();

                // surface normals for this patch
                tmp<Field<vector> > tsn = patch.Sf();
                const Field<vector> sn = tsn();

                // manually do work of addBoundaryDiag(), addBoundarySource()
                forAll(ps, faceI){
                    label c = patch.faceCells()[faceI];

                    sDiag2[c]   += cmptMultiply(ic[faceI], sn[faceI]);
                    sSource2[c] -= cmptMultiply(bc[faceI], sn[faceI]);
                }
            }


            // transfer coeffs to block matrix

            // ...diagonal elements
            forAll(d,i){
                d[i](3,0) = sDiag2[i].x();
                d[i](3,1) = sDiag2[i].y();
                d[i](3,2) = sDiag2[i].z();
                d[i](3,3) = TMdiag[i];
            }

            // ...lower elements
            forAll(l,i){
                l[i](3,0) = sLower2[i].x();
                l[i](3,1) = sLower2[i].y();
                l[i](3,2) = sLower2[i].z();
                l[i](3,3) = TMlower[i];
            }

            // ...upper elements
            forAll(u,i){
                u[i](3,0) = sUpper2[i].x();
                u[i](3,1) = sUpper2[i].y();
                u[i](3,2) = sUpper2[i].z();
                u[i](3,3) = TMupper[i];
            }

            // ...source elements
            forAll(S,i){
                S[i](3) =   sSource2[i].x()
                          + sSource2[i].y()
                          + sSource2[i].z()
                          + TMsource[i];
            }



            /**************************************************************
             * SOLVE SYSTEM
             *************************************************************/

            BlockSolverPerformance<vector4> solverPerf =
                BlockLduSolver<vector4>::New
                (
                    block_sT.name(),
                    M,
                    mesh.solutionDict().solver(block_sT.name())
                )->solve(block_sT, S);

            solverPerf.print();


            /**************************************************************
             * RETRIEVE SOLUTION
             *************************************************************/

            tmp<scalarField> tsx = s.internalField().component(0);
            tmp<scalarField> tsy = s.internalField().component(1);
            tmp<scalarField> tsz = s.internalField().component(2);
            scalarField& sx = tsx();
            scalarField& sy = tsy();
            scalarField& sz = tsz();
            blockMatrixTools::blockRetrieve(0, sx, block_sT);
            blockMatrixTools::blockRetrieve(1, sy, block_sT);
            blockMatrixTools::blockRetrieve(2, sz, block_sT);
            s.internalField().replace(0,sx);
            s.internalField().replace(1,sy);
            s.internalField().replace(2,sz);

            blockMatrixTools::blockRetrieve(3, Theta.internalField(), block_sT);

            s.correctBoundaryConditions();
            Theta.correctBoundaryConditions();

        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
