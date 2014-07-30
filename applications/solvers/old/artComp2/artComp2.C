/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    artComp2 

Description
    artificial "Compressibility" solver for Heatflux s and Temperatur Theta

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    // erzeugt Zeit-Objekt?
    #include "createTime.H"
    // Netz generieren (?)
    #include "createMesh.H"
    // Felder für Theta und s erzeugen (?)
    #include "createFields.H"
    // (?)
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
   
    const double epsilon = 1e(-3);
    double Theta_residual,
	   s_residual,
           max_residual;
 
    while (max_residual > epsilon) 
        runTime++;
        //counter++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

	// Gleichung (2)
        fvScalarMatrix ThetaEqn
        (
             (1/beta)*fvm::ddt(Theta)
           + fvc::div(s)
        );
	// kein Plan was diese Zeile macht
        ThetaEqn.setReference(ThetaRefCell, ThetaRefValue);
        volScalarField Theta_old = Theta;
	// Update Theta
        ThetaEqn.solve();
        Theta_residual = max(mag(Theta_old - Theta)).value();

	// Gleichung (1)
        fvVectorMatrix sEqn
        (
              fvm::ddt(s)
            - fvm::laplacian(s)
        );

        volVectorField s_old = s;
	// Gleichung (1) lösen mit neuem Theta --> Update s
        solve(sEqn == -fvc::grad(s));
	// maximale Veränderun
        s_residual = max(mag(s_old - s)).value();
        Info << "Residual in s: " << s_residual << " Residual in Theta: "
        << Theta_residual << endl;

        max_residual = max(s_residual, Theta_residual);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
