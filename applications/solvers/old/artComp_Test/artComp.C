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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H" // includes a lot of useful functions like fvm, fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const double epsilon = 0.000001;
    double 	s_residual = 2*epsilon,
		 		Theta_residual = 2*epsilon,
		 		residual_max = 2*epsilon,
 		 		div_s_MAX = 2*epsilon;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop() && epsilon < residual_max)
    {
       Info << "Time = " << runTime.timeName() << nl << endl;

       fvVectorMatrix sEqn
       (
           fvm::ddt(s)
         - fvm::laplacian(nu, s)
       );
		// save previous s-field
		volVectorField s_old = s;
		// Update s
   	solve(sEqn == -fvc::grad(Theta));
		// calculate residual
	   s_residual = max(mag(s_old-s)).value();

		volScalarField div_s = fvc::div(s);
		div_s_MAX = max(mag(div_s)).value();
	   //	Info << "div_s.internalField()[1]: " << div_s.internalField()[1] << endl;
		fvScalarMatrix ThetaEqn
		(
		      (1/beta)*fvm::ddt(Theta)
		);
   	ThetaEqn.setReference(ThetaRefCell, ThetaRefValue);
		// save previous Theta-field
		volScalarField Theta_old = Theta;	
		// Update Theta
		solve(ThetaEqn == -fvc::div(s));
		// calculate residual
		Theta_residual = max(mag(Theta_old-Theta)).value();
	
		Info << "Residual in s: " << s_residual << endl <<
		"Residual in Theta: " << Theta_residual << endl <<
		"div_s.internalField_MAX: " << div_s_MAX << endl;

		residual_max = max(Theta_residual, max(s_residual, div_s_MAX));
	
	        runTime.write();

	        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
  	          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
  	          << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
