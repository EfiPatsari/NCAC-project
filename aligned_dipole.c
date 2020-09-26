//-----------------------------------------------------------------------------------------------------------------------------------
// 		 EFI PATSARI  --   vortex - flux tube interaction
// Expressions of the magnetic field components in spherical coordinates can be found here: https://academic.oup.com/mnras/article/385/1/531/1033831#18455886
// (Equations A13 - A14)
// The average velocity lag expressions are explained here: https://academic.oup.com/mnras/article/421/3/2682/1747056#31308559
//-----------------------------------------------------------------------------------------------------------------------------------

# include <stdio.h>
#include <stdlib.h>
#include <math.h>
//-----------------------------------------------------------------------------------------------------------------------------------
// 		   CONSTANTS DECLARATION
//-----------------------------------------------------------------------------------------------------------------------------------
#define a       0	      // parameter a - relevant for the strength of the toroidal field
#define Bsurf 	1e12          // surface magnetic field in G
#define Rstar   10            // star radius in  km
#define Rcore   9             // star's core radius in km
#define  theta     (atan2(varpi, z))		// the polar angle of spherical coordinates as function of cylindrical coordinates
#define  r         (sqrt((pow(varpi, 2.) + pow(z, 2.))/pow(Rstar,2.)))      // spherical radius as function of cylindrical coordinates
#define angle       0.05 		// cutoff angle
//-----------------------------------------------------------------------------------------------------------------------------------
//		  FUNCTIONS DECLARATION
//-----------------------------------------------------------------------------------------------------------------------------------
// Expressing the components of the dipole magnetic field in cylindrical coordinates

double Bvarpi(double varpi, double z) 	// polar component of the magnetic field
{
    return  -Bsurf*varpi*z*0.2 /pow(Rstar,2)    ;
}
double Bzeta(double varpi, double z) 	// z-component of the magnetic field
{
    return  Bsurf*(-5.*pow(Rstar,2.) + 6.*pow(varpi,2) +3.*pow(z,2))/(15*pow(Rstar,2))		;
}
double Bphi(double varpi, double z) 	 // toroidal component of the magnetic field
{
    if ( r<=0.85	&&	r>=0.75 	&&	theta>=(M_PI/2. - angle)	&&	theta<=(M_PI/2. + angle) )
     		{
     		   return ( Bsurf*a );
     		 }
       else   {   return (0.);   }
}
//-----------------------------------------------------------------------------------------------------------------------------------
//		   MAIN 
//-----------------------------------------------------------------------------------------------------------------------------------
int main (void)
{
    double varpi  ; 		// polar radius - calculated in km
    double z ;   	   	// z cylindrical distance - calculated in km
    double length;		// vortex-line length
    double h = 0.01;  		// step of the calculation, in km
    double Bp ;			// rename the polar-component of the magnetic field
    double Bz ;			// rename the z-component of the magnetic field
    double Bt ;			// rename the toroidal magnetic field
   double Bmod ;		// calculate the magnetic field modulus
   double om_crit; 		// initial value of the \Delta\Omega_{crit}
   double av_lag = 0.; 		// initial value of the average lag \Delta \bar{\Omega}
//-----------------------------------------------------------------------------------------------------------------------------------
	FILE *fp = fopen("AISD.txt", "w");  // abbreviation for: Average lag / Isotropic angles / Standard Dipole
	if (fp == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}
//-----------------------------------------------------------------------------------------------------------------------------------
   double cot_psi;    		// cotangent of angle y between z-axis and magnetic field = Bz/Bp
   double psi = 0.1;		// minimum angle allowed - to exclude infinite pinning or no pinning at all
   double x =1./tan(psi);
//-----------------------------------------------------------------------------------------------------------------------------------
  for  (varpi=0. ; varpi<=Rstar ; varpi = varpi + h)
  { 		om_crit = 0.;
  	for (z = 0. ; z <= pow( Rcore*Rcore - varpi*varpi,0.5)  ; z = z+h)
   		{
 			Bp = Bvarpi(varpi,z); 				// rename the polar-component of the magnetic field
 			Bz = Bzeta(varpi, z);				// rename the z-component of the magnetic field
 			Bt = 0.;  //Bphi(varpi, z);					// rename the toroidal magnetic field
   			Bmod= sqrt(pow(Bp,2.) + pow(Bz,2.) + pow(Bt,2.) );
   				cot_psi = Bz/Bp;
   				if (cot_psi > x) 				// positive cutoff angle
   				{cot_psi = x; }
//   				if(cot_psi <-x)					// negative cutoff angle
//   				{cot_psi = -x;}
   			if (r <= 0.9)  		// make sure we're inside the core
   			{
   			  om_crit  =  om_crit + h*0.1*sqrt(Bmod*(1.e-12));//*(cot_psi)   ;   // Calculate the critical lag \Delta\Omega_{crit}
   		         }

    		}
	length = sqrt( pow(Rstar,2.) - pow(varpi,2.) );	 // Calculate the length of the vortex
   	av_lag = om_crit/length   ;			 // Calculate the average lag \Delta\Omega_{crit}
	fprintf(fp, "%e  	 %e     \n",  varpi, av_lag);
   }
//-----------------------------------------------------------------------------------------------------------------------------------	
fclose(fp);

return 0;
}
