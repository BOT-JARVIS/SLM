
        forAll(mesh.C(), i)
        {
	    if(alpha2[i] > alpha3[i])
	    { 
	     factor2[i]  = (scalar(2.0)*rhoCp[i])/(rho1.value()*cp1.value() + rho2.value()*cp2[i]); 
	    }
	    else
	    {  
	     factor2[i]  = (scalar(2.0)*rhoCp[i])/(rho1.value()*cp1.value() + rho3.value()*cp3[i]);
	    }
        }
        
        vaploss2 = 0.82 * ((L_v2 * M2) / (sqrt(2 * 3.14 * M2 * R * T))) * 0.54 * P_0 * exp(( L_v2 * M2 * (T - T_v2))/( R * T * T_v2));
        vaploss3 = 0.82 * ((L_v3 * M3) / (sqrt(2 * 3.14 * M3 * R * T))) * 0.54 * P_0 * exp(( L_v3 * M3 * (T - T_v3))/( R * T * T_v3));
        vaploss = alpha2*vaploss2 + alpha3*vaploss3;
        radloss2 = E2*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0)));
        radloss3 = E3*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0)));
        radloss = alpha2*radloss2 + alpha3*radloss3;
        convloss = h * (T - Tref);
        
  /*      forAll(mesh.C(), i)
        {
	    if(alpha2[i] > alpha3[i])
	    { 
	     vaploss[i] = vaploss2[i];
	     radloss[i] = radloss2[i];
	    }
	    else
	    {  
	     vaploss[i] = vaploss3[i];
	     radloss[i] = radloss3[i];
	    }
        }  */        
          
          
        
