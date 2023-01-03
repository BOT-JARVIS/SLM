
cp2 = cp2s*(1-lf) + cp2l*lf;
cp3 = cp3s*(1-lf) + cp3l*lf;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

{
    volScalarField rhoNew(alpha1*rho1 + alpha2*rho2 + alpha3*rho3);

    //solve(fvm::ddt(rho) + fvc::div(rhoPhi));
    //Info<< "density error = "
    //    << max((mag(rho - rhoNew)/mag(rhoNew))().primitiveField()) << endl;

    rho == rhoNew;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


rhoCp = alpha1*rho1*cp1 +  alpha2*rho2*cp2 + alpha3*rho3*cp3;
k2 = k2s*(1-lf) + k2l*lf;
k3 = k3s*(1-lf) + k3l*lf;
k = alpha1*k1 +  alpha2*k2 + alpha3*k3;
kappaf =  fvc::interpolate(k);
beta = alpha1*beta1 +  alpha2*beta2 +  alpha3*beta3;

forAll(mesh.C(), i)
        {
		if(alpha1[i] > 0.9)
		{
		   B[i] = 1;
		}
	     else  
	        {
		  B[i] = lf[i];
	         }    
        }
    
A =   -C*sqr(scalar(1)-B)/ (pow(B,scalar(3))+b);






