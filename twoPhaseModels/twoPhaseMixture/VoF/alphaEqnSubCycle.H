if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), 0)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
            #include "alphaEqn1.H"
	    #include "alphaEqn2.H"
	    #include "alphaEqn3.H"
	    #include "alphaEqn4.H"
	    #include "alphaEqn5.H"  //jhsu

	    #include "alphaEqn6.H"  //jhsu

	    #include "alphaEqn7.H"  //jhsu

	    #include "alphaEqn8.H"  //jhsu
	alpha1 = max(bubble1, bubble2);
	alpha1 = max(alpha1, bubble3);
	alpha1 = max(alpha1, bubble4);
	alpha1 = max(alpha1, bubble5);
	alpha1 = max(alpha1, bubble6);
	alpha1 = max(alpha1, bubble7);
	alpha1 = max(alpha1, bubble8);

	alpha2 = 1.0 - alpha1;
	#include "rhofs.H"
	volScalarField alphaRho = (alpha1 - (bubble1 + bubble2 + bubble3 + bubble4 + bubble5 + bubble6 + bubble7 + bubble8)) * (rho1f - rho2f);

	rhoPhi += fvc::interpolate(alphaRho) * phi;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
            #include "alphaEqn1.H"
            #include "alphaEqn2.H"
            #include "alphaEqn3.H"
            #include "alphaEqn4.H"
            #include "alphaEqn5.H"  //jhsu

            #include "alphaEqn6.H"  //jhsu

            #include "alphaEqn7.H"  //jhsu

            #include "alphaEqn8.H"  //jhsu
        alpha1 = max(bubble1, bubble2);
        alpha1 = max(alpha1, bubble3);
        alpha1 = max(alpha1, bubble4);
        alpha1 = max(alpha1, bubble5);
        alpha1 = max(alpha1, bubble6);
        alpha1 = max(alpha1, bubble7);
        alpha1 = max(alpha1, bubble8);

        alpha2 = 1.0 - alpha1;
        #include "rhofs.H"
        volScalarField alphaRho = (alpha1 - (bubble1 + bubble2 + bubble3 + bubble4 + bubble5 + bubble6 + bubble7 + bubble8)) * (rho1f - rho2f);

        rhoPhi += fvc::interpolate(alphaRho) * phi;

}


mixture.correct();

rho == alpha1*rho1 + alpha2*rho2;


 Info<< "Phase-1 volume fraction = "
     << alpha1.weightedAverage(mesh.Vsc()).value()
     << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
     << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
     <<  endl;

