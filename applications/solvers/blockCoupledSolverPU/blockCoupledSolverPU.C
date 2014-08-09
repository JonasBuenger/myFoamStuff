#include"fvCFD.H"
//#include"singlePhaseTransportModel.H"
//#include"RASModel.H"

//Add the support for arbitry vector size types
#include"VectorNFieldTypes.H"
#include"volVectorNFields.H"

//Add the support for block matrices
#include"blockLduSolvers.H"
#include"blockVectorNMatrices.H"

//Add the utilities written for block Matrixt ransformations
#include"blockMatrixTools.H"

//Surface interpolation schemes to get the weights
#include"surfaceInterpolationScheme.H"
//*************************************//

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info << "\n Starting timeloop \n" << endl;

    while(runTime.loop())
    {
        Info << "Time=" << runTime.timeName() << nl << endl;

        //Theta.storePrevIter();
        //s.storePrevIter();

        const surfaceVectorField& Sf = Theta.mesh().Sf();
        const unallocLabelList& owner = mesh.owner();
        const unallocLabelList& neighbour = mesh.neighbour();

        {

            //===========================================================================//
            //Velocity equation (LHS) matrix
            //===========================================================================//
            tmp<fvVectorMatrix> UEqnLHS
            (
                fvm::ddt(s)
              - fvm::laplacian(nu,s)
              + fvm::Sp(corrDim,s)
              //+ turbulence->divDevReff(s)
            );

            //Matrixblock
            BlockLduMatrix<vector4> blockM(mesh);

            //Diagonalissetseparately
            Field<tensor4>& d = blockM.diag().asSquare();

            //Off-diagonalalsoassquare
            Field<tensor4>& u = blockM.upper().asSquare();
            Field<tensor4>& l = blockM.lower().asSquare();

            //Sourcetermfortheblockmatrix
            Field<vector4> blockB (mesh.nCells(),vector4::zero);

            //Add the boundary contributions for the velocity equation
            tmp<scalarField> tdiag = UEqnLHS().D();
            scalarField& diag = tdiag();
            scalarField& upper = UEqnLHS().upper();
            scalarField& lower = UEqnLHS().lower();

            //Add diagonal boundary contribution
            //This is automatically done when you do UEqnLHS().D();
            //UEqnLHS().addBoundaryDiag(diag,0);
            //Add source boundary contribution
            vectorField& source = UEqnLHS().source();
            UEqnLHS().addBoundarySource(source, false);

            //===========================================================================//
            //Pressuregradientmatrix
            //===========================================================================//
            //Interpolationschemeforthepressureweights
            tmp<surfaceInterpolationScheme<scalar> >
            tinterpScheme_
            (
                surfaceInterpolationScheme<scalar>::New
                (
                    mesh,
                    mesh.schemesDict().interpolationScheme("grad(Theta)")
                )
            );

            //Pressure gradient contributions - corresponds to an imThetaLicit
            //gradient operator
            tmp<vectorField> tThetaUv = tmp<vectorField>
            (
                new vectorField(upper.size(), pTraits<vector>::zero)
            );
            vectorField& ThetaUv = tThetaUv();
            tmp<vectorField> tThetaLv = tmp<vectorField>
            (
                new vectorField(lower.size(), pTraits<vector>::zero)
            );
            vectorField& ThetaLv = tThetaLv();
            tmp<vectorField> tThetaSv = tmp<vectorField>
            (
                new vectorField(source.size(), pTraits<vector>::zero)
            );
            vectorField& ThetaSv = tThetaSv();

            tmp<vectorField>tThetaDv=tmp<vectorField>
            (
                new vectorField(diag.size(), pTraits<vector>::zero)
            );
            vectorField& ThetaDv = tThetaDv();

            //2)Use interpolationweights to assemble the contributions
            tmp<surfaceScalarField> tweights = tinterpScheme_().weights(Theta);
            const surfaceScalarField& weights = tweights();
            for(int i=0; i<owner.size(); i++)
            {
                int o = owner[i];
                int n = neighbour[i];
                scalar w = weights.internalField()[i];
                vector s = Sf[i];

                ThetaDv[o] += s*w;
                ThetaDv[n] -= s*(1-w);
                ThetaLv[i] =- s*w;
                ThetaUv[i] = s*(1-w);
            }

            //Get boundary condition contributions for pressure grad(Theta)
            Theta.boundaryField().updateCoeffs();
            forAll(Theta.boundaryField(),patchI)
            {
                //PresentfvPatchField
                fvPatchField<scalar>& fv = Theta.boundaryField()[patchI];

                //Retrievetheweightsfortheboundary
                const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

                //Contributionsfromtheboundarycoefficients
                tmp<Field<scalar> > tic = fv.valueInternalCoeffs(pw);
                Field<scalar>& ic = tic();
                tmp<Field<scalar> > tbc = fv.valueBoundaryCoeffs(pw);
                Field<scalar>& bc = tbc();

                //GetthefvPatchonly
                const fvPatch& patch = fv.patch();

                //Surfacenormalsforthispatch
                tmp<Field<vector> > tsn = patch.Sf();
                Field<vector> sn = tsn();

                //Manuallyaddthecontributionsfromtheboundary
                //ThiswhathappenswithaddBoundaryDiag,addBoundarySource
                forAll(fv,facei)
                {
                    label c = patch.faceCells()[facei];
                    ThetaDv[c] += ic[facei]*sn[facei];
                    ThetaSv[c] -= bc[facei]*sn[facei];
                }
        }
        //===========================================================================//
        //Assemblemomentumequation
        //===========================================================================//
        //Assemblethemomentumequationcontributions


        forAll(d,i)
        {
            d[i](0,0) = diag[i];
            d[i](1,1) = diag[i];
            d[i](2,2) = diag[i];
            d[i](0,3) = ThetaDv[i].x();
            d[i](1,3) = ThetaDv[i].y();
            d[i](2,3) = ThetaDv[i].z();
        }
        forAll(l,i)
        {
            l[i](0,0) = lower[i];
            l[i](1,1) = lower[i];
            l[i](2,2) = lower[i];
            l[i](0,3) = ThetaLv[i].x();
            l[i](1,3) = ThetaLv[i].y();
            l[i](2,3) = ThetaLv[i].z();
        }
        forAll(u,i)
        {
            u[i](0,0) = upper[i];
            u[i](1,1) = upper[i];
            u[i](2,2) = upper[i];
            u[i](0,3) = ThetaUv[i].x();
            u[i](1,3) = ThetaUv[i].y();
            u[i](2,3) = ThetaUv[i].z();
        }
        forAll(blockB,i)
        {
            blockB[i](0) = source[i].x() + ThetaSv[i].x();
            blockB[i](1) = source[i].y() + ThetaSv[i].y();
            blockB[i](2) = source[i].z() + ThetaSv[i].z();
        }

        //===========================================================================//
        //Create implicit velocity (LHS) for continuity equation
        //===========================================================================//
        //Again an implicit version not existing, now the div operator
        tmp<surfaceInterpolationScheme<scalar> >
        UtinterpScheme_
        (
            surfaceInterpolationScheme<scalar>::New
            (
                mesh,
                mesh.schemesDict().interpolationScheme("div(s)(imThetaLicit)")
            )
        );

        //1)Set up diagonal, source, upper and lower
        tmp<vectorField> tMUpper = tmp<vectorField>
            ( new vectorField(upper.size(), pTraits<vector>::zero));
        vectorField& MUpper = tMUpper();
        tmp<vectorField> tMLower = tmp<vectorField>
            ( new vectorField(lower.size(), pTraits<vector>::zero));
        vectorField& MLower = tMLower();
        tmp<vectorField> tMDiag = tmp<vectorField>
            ( new vectorField(diag.size(), pTraits<vector>::zero));
        vectorField& MDiag = tMDiag();
        tmp<vectorField> tMSource = tmp<vectorField>
        (
            new vectorField
            (
                source.component(0)().size(), pTraits<vector>::zero
            )
        );
        vectorField& MSource = tMSource();

        //2) Use interpolationweights to assemble the contributions
        const surfaceInterpolationScheme<double>& UinterpScheme_ = UtinterpScheme_();
        tmp<surfaceScalarField> tMweights = UinterpScheme_.weights(mag(s));
        const surfaceScalarField& Mweights = tMweights();
        for(int i=0; i<owner.size(); i++)
        {
            int o = owner[i];
            int n = neighbour[i];
            scalar w = Mweights.internalField()[i];
            vector sn = Sf[i];
            MDiag[o] += sn*w;
            MDiag[n] -= sn*(1-w);
            MLower[i] = -sn*w;
            MUpper[i] = sn*(1-w);
        }

        //Get boundary condition contributions for the pressure grad(Theta)
        s.boundaryField().updateCoeffs();
        forAll(s.boundaryField(),patchI)
        {
            //Present fvPatchField
            fvPatchField<vector>& fv = s.boundaryField()[patchI];

            //Retrieve the weights for the boundary
            const fvsPatchScalarField& Mw =
                Mweights.boundaryField()[patchI];

            //Contributions from the boundary coefficients
            tmp<Field<vector> > tic = fv.valueInternalCoeffs(Mw);
            Field<vector>& ic = tic();
            tmp<Field<vector> > tbc = fv.valueBoundaryCoeffs(Mw);
            Field<vector>& bc = tbc();
            //Get the fvPatch only
            const fvPatch& patch = fv.patch();
            //Surfacenormals for this patch
            tmp<Field<vector> > tsn = patch.Sf();
            Field<vector> sn = tsn();

            //Manually add the contributions from the boundary
            //This what happens with addBoundaryDiag, addBoundarySource
            forAll(fv,facei)
            {
                label c = patch.faceCells()[facei];
                MDiag[c] += cmptMultiply(ic[facei],sn[facei]);
                MSource[c] -= cmptMultiply(bc[facei],sn[facei]);
            }
        }

        //===========================================================================//
        //Create explicit and implict pressure parts for continuity equation
        //===========================================================================//
        //Pressure parts of the mass equation
        tmp<volScalarField> tA = UEqnLHS().A();
        volScalarField& A = tA();
        tmp<volVectorField> texp = fvc::grad(Theta);
        volVectorField& exp = texp();
        tmp<volVectorField> texp2 = exp/A;
        volVectorField exp2 = texp2();
        tmp<fvScalarMatrix> MEqnLHSp
        (
                (1/beta) *
                fvm::ddt(Theta)
              //  ==
              //- fvc::div(exp2)
        );

        //Add the boundary contributions
        scalarField& pMdiag = MEqnLHSp().diag();
        scalarField& pMupper = MEqnLHSp().upper();
        scalarField& pMlower = MEqnLHSp().lower();

        //Add diagonal boundary contribution
        MEqnLHSp().addBoundaryDiag(pMdiag,0);

        //Add source boundary contribution
        scalarField& pMsource = MEqnLHSp().source();
        MEqnLHSp().addBoundarySource(pMsource,false);

        //===========================================================================//
        //Assemble mass equation
        //===========================================================================//
        forAll(d,i)
        {
            d[i](3,0) = MDiag[i].x();
            d[i](3,1) = MDiag[i].y();
            d[i](3,2) = MDiag[i].z();
            d[i](3,3) = pMdiag[i];
        }
        forAll(l,i)
        {
            l[i](3,0) = MLower[i].x();
            l[i](3,1) = MLower[i].y();
            l[i](3,2) = MLower[i].z();
            l[i](3,3) = pMlower[i];
        }
        forAll(u,i)
        {
            u[i](3,0) = MUpper[i].x();
            u[i](3,1) = MUpper[i].y();
            u[i](3,2) = MUpper[i].z();
            u[i](3,3) = pMupper[i];
        }
        forAll(blockB,i)
        {
            blockB[i](3) = MSource[i].x()
            + MSource[i].y()
            + MSource[i].z()
            + pMsource[i];
        }

        //===========================================================================//
        //Solve the block matrix
        //===========================================================================//
        BlockSolverPerformance<vector4> solverPerf =
            BlockLduSolver<vector4>::New
            (
                block_sT.name(),
                blockM,
                mesh.solutionDict().solver(block_sT.name())
            )->solve(block_sT,blockB);
        solverPerf.print();

        //===========================================================================//
        //Retrieve the solution and update for next iteration
        //===========================================================================//

        tmp<scalarField> tUx = s.internalField().component(0);
        scalarField& Ux = tUx();
        blockMatrixTools::blockRetrieve(0, Ux, block_sT);
        s.internalField().replace(0, Ux);

        tmp<scalarField> tUy = s.internalField().component(1);
        scalarField& Uy = tUy();
        blockMatrixTools::blockRetrieve(1, Uy, block_sT);
        s.internalField().replace(1, Uy);

        tmp<scalarField> tUz = s.internalField().component(2);
        scalarField& Uz = tUz();
        blockMatrixTools::blockRetrieve(2, Uz, block_sT);
        s.internalField().replace(2, Uz);

        blockMatrixTools::blockRetrieve(3, Theta.internalField(), block_sT);

        // UEqnLHS.clear();

        // Theta.relax();

        s.correctBoundaryConditions();
        Theta.correctBoundaryConditions();

        }

        runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
            << "ClockTime= " << runTime.elapsedClockTime() << "s"
            << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}
