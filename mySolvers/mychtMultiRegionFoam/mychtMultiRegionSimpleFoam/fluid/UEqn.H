    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    UEqn =
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevRhoReff(U)
     ==
        pow(scalar(1)/(Porosityv+SMALL), 2)*fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
        ==
            fvc::reconstruct
            (
                (
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
