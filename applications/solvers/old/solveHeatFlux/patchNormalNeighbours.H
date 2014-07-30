#ifndef PATCHNORMALNEIGHBOURS_H
#define PATCHNORMALNEIGHBOURS_H

#include "fvPatch.H"
#include "volFields.H"
#include <utility>
#include "fvMesh.H"
#include "simpleMatrix.H"


namespace Foam{

class patchNormalNeighbours
{
private:
    int _size;
    int _levels;
    const fvMesh& _mesh;
    const fvPatch& _patch;
    label** _ID;
    vector** _C;
    double** _d;
    double*** _interpolationCoeffs;

public:
    patchNormalNeighbours(const fvPatch& p);
    void updateData();
    void set_ID();
    void set_C();
    void set_d();
    void set_interpolationCoeffs();
    label getNeighbour(const vector direction, const label cellID);
    double get_d(const int patchFace, const int level);
    label get_ID(const int patchFace, const int level);
    simpleMatrix<scalar> get_interpolationMatrix(const int patchFace, const int order);
};

}
#endif // PATCHNORMALNEIGHBOURS_H
