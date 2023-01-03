    
    gradAlpha = fvc::grad(1-alpha1);              // air is 1 so, 1-alpha1
    nHat = gradAlpha/(mag(gradAlpha) + mixture.deltaN());
    kappa = mixture.K();
    
    forAll(mesh.C(), i)
    {
		if(alpha2[i] > alpha3[i])
			{ 
			 factor1[i]  = (scalar(2.0)*rho[i])/(rho1.value() + rho2.value()); 
			}
		     else
			{  
			factor1[i]  = (scalar(2.0)*rho[i])/(rho1.value() + rho3.value());
			}
    }
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      - fvm::Sp(A,U)
      - dsigmadT2 * (fvc::grad(T) - nHat * (nHat & fvc::grad(T))) * mag(fvc::grad(alpha1)) * alpha2 * factor1
      - dsigmadT3 * (fvc::grad(T) - nHat * (nHat & fvc::grad(T))) * mag(fvc::grad(alpha1)) * alpha3 * factor1
      - 0.54 * P_0 * exp(( L_v2 * M2 * (T - T_v2))/( R * T * T_v2)) * nHat * mag(fvc::grad(alpha1)) * alpha2* factor1
      - 0.54 * P_0 * exp(( L_v3 * M3 * (T - T_v3))/( R * T * T_v3)) * nHat * mag(fvc::grad(alpha1))* alpha3 * factor1  
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho* beta * (T - Tref)) * (g & mesh.Sf())
                +
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
