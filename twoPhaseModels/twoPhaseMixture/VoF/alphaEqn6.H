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

        fvScalarMatrix bubble6Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(bubble6)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(bubble6)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, bubble6)
       // - fvm::Sp(fvc::ddt(dimensionedScalar(dimless, 1), mesh)
       //           + fvc::div(phiCN), bubble6)
         ==
            Su + fvm::Sp(Sp + divU, bubble6)
        );

        bubble6Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << bubble6.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << bubble6.name() << ") = " << min(bubble6).value()
            << "  Max(" << bubble6.name() << ") = " << max(bubble6).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(bubble6Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                bubble6,
                alphaPhi10,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;
	alpha1 =  bubble6;
        alpha2 = 1.0 - bubble6;

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
                (cnCoeff*bubble6 + (1.0 - cnCoeff)*bubble6.oldTime())(),
                compressionScheme.rewind()
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField bubble60("bubble60", bubble6);

            MULES::correct
            (
                geometricOneField(),
                bubble6,
                talphaPhi1Un(),
                talphaPhi1Corr.ref(),
                Sp,
                (-Sp*bubble6)(),
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
                bubble6 = 0.5*bubble6 + 0.5*bubble60;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                bubble6,
                phiCN,
                alphaPhi10,
                Sp,
                (Su + divU*min(bubble6(), scalar(1)))(),
                oneField(),
                zeroField()
            );
        }

	alpha1 =  bubble6;
	
        alpha2 = 1.0 - bubble6;

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
        << bubble6.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << bubble6.name() << ") = " << min(bubble6).value()
        << "  Max(" << bubble6.name() << ") = " << max(bubble6).value()
        << endl;
	CSF += mixture.surfaceTensionForce();
}
