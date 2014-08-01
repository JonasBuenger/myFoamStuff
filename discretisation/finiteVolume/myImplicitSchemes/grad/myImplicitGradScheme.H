#ifndef MYIMPLICITGRADSCHEME_H
#define MYIMPLICITGRADSCHEME_H

#include "fvCFD.H"

namespace Foam{

class myImplicitGradScheme{

private:
    volScalarField& vf_;
    const fvMesh& mesh_;
    vectorField diag_;
    vectorField upper_;
    vectorField lower_;
    vectorField source_;

public:
    myImplicitGradScheme(volScalarField&);
    void updateCoeffs();
    const vectorField& coeffsDiag();
    const vectorField& coeffsUpper();
    const vectorField& coeffsLower();
    const vectorField& coeffsSource();

};


}
#endif // MYIMPLICITGRADSCHEME_H
