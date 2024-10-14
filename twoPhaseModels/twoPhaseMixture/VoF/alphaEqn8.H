{
    #include "alphaScheme.H"

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    if (MULESCorr)
    {
        #include "alphaSuSp.H"

        fvScalarMatrix bubble8Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(bubble8)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(bubble8)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, bubble8)
       // - fvm::Sp(fvc::ddt(dimensionedScalar(dimless, 1), mesh)
       //           + fvc::div(phiCN), bubble8)
         ==
            Su + fvm::Sp(Sp + divU, bubble8)
        );

        bubble8Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << bubble8.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << bubble8.name() << ") = " << min(bubble8).value()
            << "  Max(" << bubble8.name() << ") = " << max(bubble8).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(bubble8Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                bubble8,
                alphaPhi10,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;
	alpha1 =  bubble8;
        alpha2 = 1.0 - bubble8;

        mixture.correct();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        #include "alphaSuSp.H"

        // Split operator
        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*bubble8 + (1.0 - cnCoeff)*bubble8.oldTime())(),
                compressionScheme.rewind()
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField bubble80("bubble80", bubble8);

            MULES::correct
            (
                geometricOneField(),
                bubble8,
                talphaPhi1Un(),
                talphaPhi1Corr.ref(),
                Sp,
                (-Sp*bubble8)(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                bubble8 = 0.5*bubble8 + 0.5*bubble80;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                bubble8,
                phiCN,
                alphaPhi10,
                Sp,
                (Su + divU*min(bubble8(), scalar(1)))(),
                oneField(),
                zeroField()
            );
        }

	alpha1 =  bubble8;
	
        alpha2 = 1.0 - bubble8;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi += alphaPhi10*(rho1f - rho2f);
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi += alphaPhi10*(rho1f - rho2f);
    }

    Info<< "Phase-1 volume fraction = "
        << bubble8.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << bubble8.name() << ") = " << min(bubble8).value()
        << "  Max(" << bubble8.name() << ") = " << max(bubble8).value()
        << endl;
	CSF += mixture.surfaceTensionForce();
}
