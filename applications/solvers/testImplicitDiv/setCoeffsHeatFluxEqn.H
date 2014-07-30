#ifndef SETCOEFFSHEATFLUXEQN_H
#define SETCOEFFSHEATFLUXEQN_H

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
forAll(l,i){
    u[i](0,0) = sUpper[i];
    u[i](1,1) = sUpper[i];
    u[i](2,2) = sUpper[i];

    u[i](0,3) = TUpper[i].x();
    u[i](1,3) = TUpper[i].y();
    u[i](2,3) = TUpper[i].z();
}

// ...source elements
forAll(l,i){
    S[i](0) = sSource[i].x() + TSource[i].x();
    S[i](1) = sSource[i].y() + TSource[i].y();
    S[i](2) = sSource[i].z() + TSource[i].z();
}


#endif // SETCOEFFSHEATFLUXEQN_H
