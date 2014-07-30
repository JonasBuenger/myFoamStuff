#include "patchNormalNeighbours.H"

namespace Foam{

patchNormalNeighbours::patchNormalNeighbours(const fvPatch& p):
    _patch(p),
    _mesh(p.boundaryMesh().mesh())
{
    _levels = 4;
    _size = _patch.size();
    _ID = new label*[_size];
    _C = new vector*[_size];
    _d = new double*[_size];
    _interpolationCoeffs = new double**[_size];
    for(int j=0; j<_size; j++){
        _ID[j] = new label[_levels];
        _C[j] = new vector[_levels];
        _d[j] = new double[_levels];
        _interpolationCoeffs[j] = new double*[_levels];
        for(int s=0; s<_levels; s++){
            _interpolationCoeffs[j][s] = new double[_levels];
        }
    }
    updateData();
}

void patchNormalNeighbours::set_ID(){
    vectorField normalVectors = _patch.nf();
    for(int i=0; i<_size; i++){
        _ID[i][0] = _patch.faceCells()[i];
        vector vec_n = normalVectors[i]/mag(normalVectors[i]);
        for(int j=1; j<_levels; j++){
            _ID[i][j] = getNeighbour(vec_n, _ID[i][j-1]);
        }
    }
}

void patchNormalNeighbours::set_C(){
    const vectorField centers = _mesh.C();
    for(int i=0; i<_size; i++){
        _C[i][0] = centers[i];//_ID[i][0]];
        for(int j=1; j<_levels; j++){
            _C[i][j] = centers[_ID[i][j]];
        }
    }
}

void patchNormalNeighbours::set_d(){
    vectorField patchCenters = _patch.Cf();
    for(int i=0; i<_size; i++){
        _d[i][0] = mag(patchCenters[i] - _C[i][0]);
        for(int j=1; j<_levels; j++){
            _d[i][j] = mag(patchCenters[i] - _C[i][j]);
        }
    }
}

void patchNormalNeighbours::set_interpolationCoeffs(){
    for(int i=0; i<_size; i++){
        for(int j=0; j<_levels; j++){
            for(int s=0; s<_levels; s++){
                _interpolationCoeffs[i][j][s] = pow(_d[i][j], s);
            }
        }
    }
}

void patchNormalNeighbours::updateData(){
    set_ID();
    set_C();
    set_d();
    set_interpolationCoeffs();
}

label patchNormalNeighbours::getNeighbour(const vector direction, const label cell_ID){

    if ((mag(direction) > 1.001) && (mag(direction) < 0.999))
        Info << "Warnung! Richtung ist nicht von Länge 1!" << endl;

    label normalNeighbour = -1;
    const labelList& neighbours = _mesh.cellCells()[cell_ID];
    scalar maxInnerProduct = -1;

    for (int i=0; i<neighbours.size(); i++){
        vector connectingVector = _mesh.C()[neighbours[i]] - _mesh.C()[cell_ID];
        scalar innerProduct = (connectingVector & direction);
        if ( innerProduct > maxInnerProduct ){
            maxInnerProduct = innerProduct;
            normalNeighbour = neighbours[i];
        }
    }

    return normalNeighbour;
}

double patchNormalNeighbours::get_d(const int patchFace, const int level){
    return _d[patchFace][level];
}

label patchNormalNeighbours::get_ID(const int patchFace, const int level){
    return _ID[patchFace][level];
}

simpleMatrix<scalar> patchNormalNeighbours::get_interpolationMatrix(const int patchFace, const int order){

    if(order>_levels)
        Info << "order to high!" << endl;

    simpleMatrix<scalar> mat(order);
    for(int r=0; r<order; r++){
        for(int c=0; c<order; c++){
            mat[r][c] = _interpolationCoeffs[patchFace][r][c];
        }
    }
    return mat;
}

}
