# include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define a       1		 // parameter a - relevant for the strength of the toroidal field
#define Bsurf 	1e12         // surface magnetic field in G
#define Rstar   10            // star radius in  km
#define Rcore   9              // star's core radius in km

#define  y        (M_PI/3.)		// inclination angle
//#define r	       ( sqrt( varpi*varpi + z*z ) )   // spherical radius as function of cylindrical coordinates


double Theta(double varpi, double phi, double z)
{
//    double ex1 = pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.);
//    double ex2 =  z * cos (y) - varpi * sin (y) * sin (phi) ;
    return (sqrt(pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.))/(z * cos (y) - varpi * sin (y) * sin (phi)));
}

double Torus(double r, double theta) // Function that appears in the regions, where the toroidal field exists
{
  if (   r >= 7.5      && r <= 8.5      && atan(theta) > (M_PI / 2. - 0.05)      && atan(theta) < (M_PI / 2. + 0.05))
      {return  a * Bsurf ;}
      else
      {return 0.;}
}


// Components of the magnetic field in cylindrical coordinates
double Bvarpi(double varpi, double phi, double z, double r, double tor) // polar component of the magnetic field
{
//   r =  sqrt( varpi*varpi + z*z )  ;
    return  ( ( cos(phi)*(-(( Bsurf*(z*cos(y) - varpi*sin(y)*sin(phi)) * pow((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.))/ pow(r,2.),1.5)*sqrt( pow(r,2.))*
             (6.* pow(r,2.) - 5.* pow(Rstar,2.))* varpi* cos(phi))/( pow((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.)),1.5)* pow(Rstar,2.)))\
         - ( Bsurf*(z*cos(y) - varpi*sin(y)*sin(phi)) *(-3.* pow(r,2.) + 5.* pow(Rstar,2.))* varpi* cos(phi))/
         ( pow(r,2.)* pow(Rstar,2.)) -
        ( 15.*tor*(z*sin( y) +  varpi* cos( y)*sin(phi)))/sqrt((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.)))))/15. +
   sin(phi)*(( Bsurf*(sqrt((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.)))*sqrt((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.))/ pow(r,2.))*sqrt( pow(r,2.))*
            (-6.* pow(r,2.) + 5.* pow(Rstar,2.)) +
            pow((z*cos(y) - varpi*sin(y)*sin(phi)) ,2.)*(-3.* pow(r,2.) + 5.* pow(Rstar,2.)))*sin( y))/
       (15.* pow(r,2.)* pow(Rstar,2.)) +
       cos( y)*(( varpi* cos(phi)*tor /
          sqrt((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.))) - ( Bsurf*(z*cos(y) - varpi*sin(y)*sin(phi)) * pow((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.))/ pow(r,2.),1.5)*sqrt( pow(r,2.))*
            (6.* pow(r,2.) - 5.* pow(Rstar,2.))*(z*sin( y) +  varpi* cos( y)*sin(phi)))/
          (15.* pow((pow(varpi*cos(phi),2.) + pow (z * sin(y)+ varpi*cos(y)*sin(phi),2.)),1.5)* pow(Rstar,2.)) -
         ( Bsurf*(z*cos(y) - varpi*sin(y)*sin(phi)) *(-3.* pow(r,2.) + 5.* pow(Rstar,2.))*(z*sin( y) +  varpi* cos( y)*sin(phi)))/
          (15.* pow(r,2.)* pow(Rstar,2.)))) ));
}


double Bphi(double varpi, double phi, double z, double r, double tor)  // toroidal component of the magnetic field
{
 //  r =  sqrt( varpi*varpi + z*z )  ;
          return ( -(sin(phi)*(-((tor*(z*sin(y) + varpi*cos(y)*sin(phi)))/
           sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))) +
        (Bsurf*varpi*cos(phi)*(z*cos(y) - varpi*sin(y)*sin(phi))*
           (3./pow(Rstar,2.) - 5./
              (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))))/15. +
        (Bsurf*varpi*cos(phi)*(z*cos(y) - varpi*sin(y)*sin(phi))*
           pow((pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))/
             (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
               pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)),1.5)*
           sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
             pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))*
           (0.3333333333333333 - (2.*
                (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                  pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)))/(5.*pow(Rstar,2.))))/
         pow(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.),1.5))
      ) + cos(phi)*((Bsurf*sin(y)*(sqrt(pow(varpi,2.)*pow(cos(phi),2.) +
              pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))*
            sqrt((pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))/
              (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)))*
            sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
              pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))*
            (5.*pow(Rstar,2.) - 6.*(pow(varpi,2.)*pow(cos(phi),2.) +
                 pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                 pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))) +
           pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)*
            (5.*pow(Rstar,2.) - 3.*(pow(varpi,2.)*pow(cos(phi),2.) +
                 pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                 pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)))))/
       (15.*pow(Rstar,2.)*(pow(varpi,2.)*pow(cos(phi),2.) +
           pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) + pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))
         ) + cos(y)*((varpi*cos(phi)*tor)/
          sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.)) +
         (Bsurf*(z*sin(y) + varpi*cos(y)*sin(phi))*(z*cos(y) - varpi*sin(y)*sin(phi))*
            (3./pow(Rstar,2.) - 5./
               (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                 pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))))/15. +
         (Bsurf*(z*sin(y) + varpi*cos(y)*sin(phi))*(z*cos(y) - varpi*sin(y)*sin(phi))*
            pow((pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))/
              (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)),1.5)*
            sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
              pow(z*cos(y) - varpi*sin(y)*sin(phi),2.))*
            (0.3333333333333333 - (2.*
                 (pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.) +
                   pow(z*cos(y) - varpi*sin(y)*sin(phi),2.)))/(5.*pow(Rstar,2.))))/
          pow(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.),1.5)
         ))

          ) ;
}
double Bzeta(double varpi, double phi, double z, double r, double tor) // z-component of the magnetic field
{
 //  r =  sqrt( varpi*varpi + z*z )  ;
    return  (-(Bsurf*cos(y)*((5.*pow(Rstar,2.) - 3.*(pow(z,2.) + pow(varpi,2.)))*
          pow(z*cos(y) - varpi*sin(y)*sin(phi),2.) +
         sqrt(pow(z,2.) + pow(varpi,2.))*(5.*pow(Rstar,2.) - 6.*(pow(z,2.) + pow(varpi,2.)))*
          sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))*
          sqrt((pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))/
            (pow(z,2.) + pow(varpi,2.)))))/(15.*pow(Rstar,2.)*(pow(z,2.) + pow(varpi,2.)))\
    + sin(y)*(-(Bsurf*(5.*pow(Rstar,2.) - 3.*(pow(z,2.) + pow(varpi,2.)))*
          (z*sin(y) + varpi*cos(y)*sin(phi))*(z*cos(y) - varpi*sin(y)*sin(phi)))/
       (15.*pow(Rstar,2.)*(pow(z,2.) + pow(varpi,2.))) +
      (varpi*cos(phi)*tor)/
       sqrt(pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.)) -
      (Bsurf*sqrt(pow(z,2.) + pow(varpi,2.))*
         (-5.*pow(Rstar,2.) + 6.*(pow(z,2.) + pow(varpi,2.)))*(z*sin(y) + varpi*cos(y)*sin(phi))*
         (z*cos(y) - varpi*sin(y)*sin(phi))*
         pow((pow(varpi,2.)*pow(cos(phi),2.) + pow(z*sin(y) + varpi*cos(y)*sin(phi),2.))/
           (pow(z,2.) + pow(varpi,2.)),1.5))/
       (15.*pow(Rstar,2.)*pow(pow(varpi,2.)*pow(cos(phi),2.) +
           pow(z*sin(y) + varpi*cos(y)*sin(phi),2.),1.5)))
           );
 }

int main (void)
{
    double varpi  ; 			   // Calculated in km
    double phi;				// Cylindrical angle
    double z ;   	   			// Calculated in km
    double length;			// vortex-line length
    double h = 0.05;  		// Step of the calculation, in km
    double Bp ;				// rename the polar-component of the magnetic field
    double Bz ;					// rename the z-component of the magnetic field
    double Bt ;					// rename the toroidal magnetic field
    double Bmod ;		//    Calculate the magnetic field modulus
    double tor;
    double theta;  		 // the polar angle of spherical coordinates as function of cylindrical coordinates
    double xx, yy ;

    double om_crit; 		// Initial value of the \Delta\Omega_{crit}
    double av_lag = 0.; 		//  Initial value of the average lag \Delta \bar{\Omega}


    double r;


	FILE *fp = fopen("AITDp3a0.txt", "w"); // Average lag / Isotropic / Tilted Dipole / angle=pi/3 /parameter a=0
	if (fp == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}

   double cot_psi;    				// cotangent of angle y between z-axis and magnetic field = Bz/Bp
   double psi = 0.1;				// minimum angle allowed to prevent infinite pinning or no pinning
   double x =1./tan(psi);


        phi = 0. ;

  for  (varpi=0. ; varpi<=Rstar ; varpi = varpi + h)
  { 		om_crit = 0.;
    	//for  (phi=0.01 ; phi<= 2.*M_PI ; phi = phi + h)
    	   {
    	           om_crit = 0.;
    	   	for (z = -pow( Rcore*Rcore - varpi*varpi,0.5) ; z <= pow( Rcore*Rcore - varpi*varpi,0.5)  ; z = z + h)
   		   {   	r =  sqrt( varpi*varpi + z*z )  ;
  			theta = Theta(varpi, phi, z);
  			tor = 0.; //Torus(r, theta);
 			Bp = Bvarpi(varpi,phi, z, r, tor); 			// rename the polar-component of the magnetic field
 			Bz = Bzeta(varpi, phi, z, r, tor);			// rename the z-component of the magnetic field
			Bt = Bphi(varpi, phi, z, r, tor);			// rename the toroidal magnetic field

   			Bmod= sqrt(pow(Bp,2.) + pow(Bz,2.) + pow(Bt,2.) );

// 			cot_psi = (Bz*cos(y) - Bp*sin(y))/(Bp*cos(y) + Bz*sin(y));
//   		        	if (cot_psi > x) 				// positive cutoff angle
//   				{cot_psi = x; }
//   				if(cot_psi <-x)					// negative cutoff angle
//   				{cot_psi = -x;}
   			if (r <= 0.9)  		// make sure we're inside the core
   			{
   			  om_crit  =  om_crit + h*0.1*sqrt(Bmod*(1.e-12)); //*(cot_psi)   ;   // Calculate the critical lag \Delta\Omega_{crit}
// 			  printf("omega crit = %e	theta = %lf		Bmod = %e \n", om_crit, theta(varpi,phi,z) , Bmod);
   		         }

    		   }
	length = 2.*sqrt( pow(Rstar,2.) - pow(varpi,2.) );	 // Calculate the length of the vortex
   	av_lag = om_crit/length   ;		// Calculate the average lag \Delta\Omega_{crit}
        xx = varpi*cos(phi);
        yy = varpi*sin(phi);
	fprintf(fp, "%e     %e   %e    \n", varpi, z, av_lag);
	}
//      fprintf(fp,"\n");
 }
fclose(fp);

return 0;
}
