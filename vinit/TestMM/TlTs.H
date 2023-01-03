
       forAll(mesh.C(), i)
        {
        
              if ((alpha3[i] >= 0) && (alpha3[i] < 0.32)) 
	          {
		    Tl[i] = -350.78*alpha3[i]+936.6;                 // line1
	          }
	          
	       else if ((alpha3[i] >= 0.32) && (alpha3[i] < 0.53)) 
	          {
		    Tl[i] = 246.66*alpha3[i]+745.41;                 // line2
	          }
	          
	       else if ((alpha3[i] >= 0.53) && (alpha3[i] < 0.84)) 
	          {
		    Tl[i] = 1548.38*alpha3[i]+55.50;                 // line3 
	          }
	          
	       else if ((alpha3[i] >= 0.84) && (alpha3[i] < 0.92)) 
	          {
		    Tl[i] = -250*alpha3[i]+1566.15;                  // line4
	          }
	          
	       else if ((alpha3[i] >= 0.92) && (alpha3[i] <= 1)) 
	          {
		    Tl[i] = 310.87*alpha3[i]+1050.145;               // line5
	          }                
        }
        
        forAll(mesh.C(), i)
        {
	       
	       if ((alpha3[i] >= 0) && (alpha3[i] < 0.06)) 
	          {
		    Ts[i] = -1871.66*alpha3[i]+933.6;                // line1           
	          }
	          
	       else if ((alpha3[i] >= 0.06) && (alpha3[i] < 0.52)) 
	          {
		    Ts[i] = 821.3;                                   // line2
	          }
	          
	       else if ((alpha3[i] >= 0.52) && (alpha3[i] < 0.53)) 
	          {
		    Ts[i] = 5185*alpha3[i]-1874.9;                   // line3
	          }
	          
	       else if ((alpha3[i] >= 0.53) && (alpha3[i] < 0.7)) 
	          {
		    Ts[i] = 873.15;                                  // line4
	          }
	          
	       else if ((alpha3[i] >= 0.7) && (alpha3[i] < 0.72)) 
	          {
		    Ts[i] = 2000*alpha3[i]-526.85;                   // line5 
	          }
	          
	       else if ((alpha3[i] >= 0.72) && (alpha3[i] < 0.75)) 
	          {
		    Ts[i] = 913.15;                                  // line6 
	          }
	          
	       else if ((alpha3[i] >= 0.75) && (alpha3[i] < 0.8)) 
	          {
		    Ts[i] = 6800*alpha3[i]-4186.85;                  // line7
	          }
	          
	       else if ((alpha3[i] >= 0.8) && (alpha3[i] < 0.84)) 
	          {
		    Ts[i] = 2500*alpha3[i]-746.85;                   // line8
	          }
	          
	       else if ((alpha3[i] >= 0.84) && (alpha3[i] < 0.92)) 
	          {
		    Ts[i] = -250*alpha3[i]+1563.15;                                // line9
	          }
	          
	       else if ((alpha3[i] >= 0.92) && (alpha3[i] <= 1)) 
	          {
		    Ts[i] = 310.875*alpha3[i]+1047.145;                                // line10
	          }                               
        }
        
        
        
        volScalarField DT = Tl-Ts;
               
          
          
        
