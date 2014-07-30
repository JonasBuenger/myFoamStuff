
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        /*
         *
         *
         ALGORITHMUS
         *
         *
         */

        runTime.write();
    }

    return 0;
}
