
#include "propertiesT.H"

////////////////////////////////////////////////////////////          Radiation          /////////////////////////////////////////////////////////////////// 
      
      forAll(mesh.C(), i)
       {
		if(alpha1[i] > 0.01)
			{ 
			 s1[i]  = alpha1[i]; 
			}
                   else 
			{
			   s1[i]  = 0;
			}
	}
      
      alphaT = 1-s1;
      radiation->correct();
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    
      int iter         = 0;
      scalar residual  = 1;
      int minTCorr     = 1;
      int maxTCorr     = 20;           //20
      scalar alphaTol = 0.01; 
        
      volScalarField lfPRELOOP = lf;  
       
            
      volScalarField Tl = alpha1*Tl2 + alpha2*Tl2 + alpha3*Tl3;

      volScalarField Ts = alpha1*Ts2 + alpha2*Ts2 + alpha3*Ts3;
      
      volScalarField DT = Tl - Ts;
       
      volScalarField L = alpha1*L2 + alpha2*L2 + alpha3*L3;
        do
        {        
        iter++;
        
        volScalarField lfold = lf;  
        volScalarField dFdT = 1 / (Tl - Ts);
        volScalarField invF = lfold * (Tl - Ts) + Ts;
	volScalarField dH =  (1-alpha1)*rho * L * pos(Tl - T) * pos(T - Ts);
	volScalarField Sp = - dH/runTime.deltaT() * dFdT;
	volScalarField Sc = dH/runTime.deltaT() * (lfPRELOOP - lf + dFdT * invF); 

        fvScalarMatrix TEqn
        (
             fvm::ddt(rhoCp, T)
          +  fvm::div(rhoCpphi, T)
          -  fvm::Sp(Sp, T)
          == fvm::laplacian(kappaf, T) + Sc + radiation->ST(T) - (vaploss+radloss+convloss)*mag(fvc::grad(alpha1)) * factor2
        );         

        TEqn.solve();
        
        lf = (lfold + (dFdT*(T-invF)));
	lf = max(lf,scalar(0));
	lf = min(scalar(1),lf);
               
        residual = max(mag((lf-lfold))).value();
 
        Info << "Iteration :"<<iter<<", max residual = " << residual << endl;
        
        
        } while ((iter < minTCorr || residual > alphaTol ) && iter <= maxTCorr);

                             
        Info << "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;
        Info << "min/max(Tl) = " << min(Tl).value() << ", " << max(Tl).value() << endl;
        Info << "min/max(Ts) = " << min(Ts).value() << ", " << max(Ts).value() << endl;
        Info << "min/max(DT) = " << min(DT).value() << ", " << max(DT).value() << endl;
        Info << "samay " << runTime.value() << endl;
        
        
        
                       
                
                
                
