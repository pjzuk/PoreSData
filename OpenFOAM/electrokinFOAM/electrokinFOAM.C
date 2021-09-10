/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    electrokinFoam

Description
    Solver for Poisson-Nernst-Planck equations.
    Written by Pawel Jan Zuk 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControls.H"
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensionedScalar nEE = (e/epsilon)*nRef;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Time = " << runTime.timeName() << endl;
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;
        Info<< "deltaT = " << runTime.deltaT().value() << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"

        scalar deltaTFact = 1.0;    

        if (adjustTimeStep)
        {
            if (CoNum > maxCo) deltaTFact = maxCo/(CoNum + SMALL)/2.; 
            if (CoNum + SMALL < maxCo/2.) deltaTFact = 1.01;
            runTime.setDeltaT
            (
                    min(runTime.deltaT().value()*deltaTFact,maxDeltaT)
            );
        }

        // --- PIMPLElike loop for ion dynamics
        for (int corr=0; corr<nOuterCorrIons; corr++)
        {
	
    	    snPhiGrad = fvc::snGrad(ePhi) * mesh.magSf();

	    muPlusSnPhiGrad = snPhiGrad*muPlus*ZPlus;   
	    fvScalarMatrix nPlusEqn
            (
                fvm::ddt(nPlus)
    	      - fvm::laplacian(DPlus,nPlus)
              - fvc::div(muPlusSnPhiGrad,nPlus)
            );    
            if (corr == nOuterCorrIons -1)
            {
                nPlusEqn.relax(1);
            }
            else
            {
                nPlusEqn.relax();
            }
	    nPlusEqn.solve();
	    nPlus.correctBoundaryConditions();

       	    Info << "nPlus equation solved \n";

            muMinusSnPhiGrad = snPhiGrad*muMinus*ZMinus;
            fvScalarMatrix nMinusEqn
            (
                fvm::ddt(nMinus)
              - fvm::laplacian(DMinus,nMinus)
              + fvc::div(muMinusSnPhiGrad,nMinus)
            );
            if (corr == nOuterCorrIons -1)
            {
                nMinusEqn.relax(1);
            }
            else
            {
                nMinusEqn.relax();
            }
	    nMinusEqn.solve();
	    nMinus.correctBoundaryConditions();

    	    Info << "nMinus equation solved \n";

            // --- PISOlike loop for electric potential
            for (int inCorr=0; inCorr<nCorrIons; inCorr++)
            {
                fvScalarMatrix PhiEqn
                (
                    fvm::laplacian(ePhi)
                  ==
                    nEE*( ZMinus*nMinus - ZPlus*nPlus )
                );
                if (corr != nOuterCorrIons -1)
                {
                    PhiEqn.relax();
                }

                PhiEqn.solve();

                Info << "ePhi equation solved \n";
            }
           // --- end of PISOlike loop for ion dynamics

        }
        // --- end of PIMPLElike loop for ion dynamics


        // --- other usefull fields calculations
        PhiGrad = fvc::grad(ePhi);
        nMinusGrad = fvc::grad(nMinus);
    	nMinusSnGrad = fvc::interpolate(nMinusGrad) & mesh.Sf();
        nPlusGrad = fvc::grad(nPlus);
    	nPlusSnGrad = fvc::interpolate(nPlusGrad) & mesh.Sf();
	    nPlusCurrent = - nPlusGrad*DPlus - nPlus*ZPlus*muPlus*PhiGrad;
	    nMinusCurrent = - nMinusGrad*DMinus + nMinus*ZMinus*muMinus*PhiGrad;
	    nPlusSnVelocity = - nPlusSnGrad*DPlus/max(mag(fvc::interpolate(nPlus)),nMinimals) - muPlusSnPhiGrad;
	    nMinusSnVelocity = - nMinusSnGrad*DMinus/max(mag(fvc::interpolate(nMinus)),nMinimals) +  muMinusSnPhiGrad;
	    netCharge = nRef*e*(nPlus*ZPlus-nMinus*ZMinus);
        // --- end of other usefull fields calculations

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
