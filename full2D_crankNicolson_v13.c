/*
January 2015

Covariant metric structure:

g11 = f[i][j]
g12 = g21 = h[i][j]*r = h[i][j]*i*dr
g22 = g[i][j]*r*r = g[i][j]*i*i*dr*dr

Contravariant metric structure:
g11upper = -g[i][j]/den
g12upper = g21upper = h[i][j]/(den*r) = h[i][j]/(den*i*dr)
g22upper = -f[i][j]/(den*r*r) = -f[i][j]/(den*i*i*dr*dr)

   where den = h[i][j]*h[i][j] - f[i][j]*g[i][j] = -det(metric)

 */

/******************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************/

const int rmax = 100;        //max radial coordinate
const int thetamax = 100;        //max radial coordinate
//const int tmax = 11;
// const int tmax = 1001;
// const int tmax = 10001;
// const int tmax = 100001;
const int tmax = 750001;
//const double dt = 0.001;  
//const double dt = 0.0001;  
//const double dt = 0.00001;  
//const double dt = 0.000001;  
//const double dt = 0.0000001;  
const double dt = 0.000101;     // regular dt (10^-4)
// const double dt = 0.0000501;    //small dt (5*10^-5)
//const double dt = 0.00000101;  //extra small dt (10^-6)
//const double dt = 0.000000501;  //xxsmall dt (5*10^-7)    
// const double dr = 0.0503;       
// const double dtheta = 0.0503;
const double dr = 0.103;       
const double dtheta = 0.103;

double fn[rmax+1][thetamax+1];       //metric component being calculated (t+1) - covariant 
double f[rmax+1][thetamax+1];        //current metric component (t) - covariant
double fp[rmax+1][thetamax+1];       //past metric (t-1) - covariant
double hn[rmax+1][thetamax+1];       //metric component being calculated (t+1) - covariant 
double h[rmax+1][thetamax+1];        //current metric component (t) - covariant
double hp[rmax+1][thetamax+1];       //past metric (t-1) - covariant
double gn[rmax+1][thetamax+1];       //metric component being calculated - covariant
double g[rmax+1][thetamax+1];        //current metric component - covariant
double gp[rmax+1][thetamax+1];       //past metric (t-1) - covariant
double v1n[rmax+1][thetamax+1];        //new radial velocity - contravariant
double v1[rmax+1][thetamax+1];         //current radial velocity - contravariant
double v1p[rmax+1][thetamax+1];        //past radial velocity - contravariant
double v2n[rmax+1][thetamax+1];        //new theta velocity - contravariant
double v2[rmax+1][thetamax+1];         //current theta velocity - contravariant
double v2p[rmax+1][thetamax+1];        //past theta velocity - contravariant
double gamma111[rmax+1][thetamax+1];  //connection coefficients
double gamma211[rmax+1][thetamax+1];  //connection coefficients
double gamma112[rmax+1][thetamax+1];  //connection coefficients
double gamma212[rmax+1][thetamax+1];  //connection coefficients
double gamma122[rmax+1][thetamax+1];  //connection coefficients
double gamma222[rmax+1][thetamax+1];  //connection coefficients
double T11[rmax+1][thetamax+1];        //Growth tensor component 
double T12[rmax+1][thetamax+1];        //Growth tensor component 
double T22[rmax+1][thetamax+1];        //Growth tensor component 
double K[rmax+1][thetamax+1];          //Scalar curvature

FILE *dynamicsData_thetaPi0,*dynamicsData_thetaPi4,*dynamicsData_thetaPi2,*dynamicsData_thetaPi,*plantLength,*info,*twoDdata,*modgrowth_Pi0,*modgrowth_Pi4,*modgrowth_Pi2,*modgrowth_Pi;
char filename[30];
char filetag[10] = "%s%d%s%d%s";

/******************************************************************/
/******************************************************************/


// Silk velocity profile
double SilkVel(double dblx)
{
  double q,n;
  n = 5.0;
  //q = 1.0 + exp(0.2*dblx-2.0);
  //q = 1.0 + exp(0.2*dblx-5.0);
  q = 1.0 + exp(0.2*dblx-7.0);
  //q = 1.0 + exp(0.2*dblx-9.0);
  //q = 1.0 + exp(0.2*dblx-11.0);
  //q = 1.0 + exp(0.2*dblx-13.0);
  //q = 1.0 + exp(0.2*dblx-15.0); 
  //return 0.1*pow(q,-1.0/n);
  //return 0.5*pow(q,-1.0/n);
  return 0.9*pow(q,-1.0/n);
}

// The InitialConditions function exploits the fact that the metric and 
// velocity fields are global variables. The switch is set in main(), in
// the function call to PlantModel.
void InitialConditions(int metricFunction, int velocityProfile)
{
  int l,m,ntheta;
  double metricParts,logisticPart,gaussianPart;

  // initialize the metric function using switches.
  switch (metricFunction) {
  case 0:        // flat metric
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 1.0;
        h[l][m] = 0.0;  
        g[l][m] = 1.0;
      }
    }
    break;
  case 1:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;    
  case 2:         // parabolic 
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = dr*dr*l*l + 1.0; 
        h[l][m] = 0.0;  
        g[l][m] = 1.0;
      }
    }
    break;  
  case 3:     // Silk velocity profile
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = SilkVel(l); 
        h[l][m] = 0.0;  
        g[l][m] = 1.0;
      }
    }
    break;  
  case 4:         // centered gaussian 
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-l*l*dr*dr/2.5) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;  
  case 5:         // gaussian/parabolic
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr/2.5) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = dr*dr*l*l; 
      }
    }
    break;    
  case 6:         // logistic
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 1.0 + 1.0/(1.0 + exp(-5.0*dr*(l-rmax/2.0))); 
        h[l][m] = 0.0; 
        g[l][m] = 0.0; 
      }
    }
    break;  
  case 7:         // flat metric with radial gaussian perturbation (sharper peak)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*6.0) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;      
  case 8:         // flat metric with radial gaussian perturbation (sharper peak)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*15.0) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;      
  case 9:         // flat metric with radial gaussian perturbation (sharper peak)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax*0.8)*(l-rmax*0.8)*dr*dr*6.0) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;  
  case 10:         // flat metric with radial gaussian perturbation (sharper peak)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax*0.2)*(l-rmax*0.2)*dr*dr*6.0) + 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 1.0; 
      }
    }
    break;  
  case 11:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 1.0; 
        h[l][m] = 0.5; 
        g[l][m] = 1.0; 
      }
    }
    break;
  case 12:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 1.0; 
        h[l][m] = 0.0; 
        g[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 1.0; 
      }
    }
    break;
  case 13:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 1.0;  
        g[l][m] = 2.0;
      }
    }
    break;
  case 14:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-50.0)*(m-50.0)*dtheta*dtheta*2.0)) + 1.0; 
        h[l][m] = 0.5; 
        g[l][m] = 1.0; 
      }
    }
    break;
  case 15:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 1.0; 
        h[l][m] = 0.5; 
        g[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-50.0)*(m-50.0)*dtheta*dtheta*2.0)) + 1.0; 
      }
    }
    break;
  case 16:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 2.0; 
        h[l][m] = 0.0;  
        g[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0) + 2.0;
      }
    }
    break;
  case 17:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 5.0*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*3.0) + 2.0; 
        h[l][m] = 0.0; 
        g[l][m] = 2.0*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*3.0) + 2.0; 
      }
    }
    break;    
  case 18:         // flat metric with radial gaussian perturbation
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.25*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*4.0)*(0.5*exp(-(m-50.0)*(m-50.0)*dtheta*dtheta*0.5)) + 1.0; 
        h[l][m] = 2.0;  
        g[l][m] = 0.25*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*4.0)*(0.5*exp(-(m-50.0)*(m-50.0)*dtheta*dtheta*0.5)) + 1.0;
      }
    }
    break;
  case 19:         // single gaussian
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.25*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*4.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.5)) + 2.0; 
        h[l][m] = 0.25*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*4.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.5)) + 1.0;  
        g[l][m] = 0.25*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*4.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.5)) + 2.0;
      }
    }
    break;
  case 20:         // single gaussian (small)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 1.0;  
        g[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 2.0;
      }
    }
    break;
  case 21:         // single gaussian (large)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 5.0*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(1.0*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(0.5*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 1.0;  
        g[l][m] = 2.0*exp(-(l-rmax/2.0)*(l-rmax/2.0)*dr*dr*2.0)*(1.0*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1)) + 2.0;
      }
    }
    break;
  case 22:         // double crescent
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax*0.25)*(l-rmax*0.25)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*0.1)) + 0.5*exp(-(l-rmax*0.75)*(l-rmax*0.75)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*0.1)) + 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax*0.25)*(l-rmax*0.25)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*0.1)) + 0.5*exp(-(l-rmax*0.75)*(l-rmax*0.75)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*0.1)) + 1.0;  
        g[l][m] = 0.5*exp(-(l-rmax*0.25)*(l-rmax*0.25)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*0.1)) + 0.5*exp(-(l-rmax*0.75)*(l-rmax*0.75)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*0.1)) + 2.0;
      }
    }
    break;  
  case 23:         // double gaussian (small perturbation)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 1.0;  
        g[l][m] = 0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 2.0;
      }
    }
    break;
  case 24:         // double gaussian (large perturbation)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        f[l][m] = 5.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(1.0*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 1.0*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 2.0; 
        h[l][m] = 0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(0.5*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 0.5*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 1.0;  
        g[l][m] = 2.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0)*(1.0*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*2.0) + 1.0*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*2.0)) + 2.0;
      }
    }
    break;
  case 25:         // four gaussians (large perturbation)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        metricParts = 1.0*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0); 
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*5.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        f[l][m] = metricParts;

        metricParts = 0.5*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 1.0;
        h[l][m] = metricParts;

        metricParts = 1.0*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 1.0*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*2.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        g[l][m] = metricParts;
      }
    }
    break;
  case 26:         // four gaussians (small perturbation)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        metricParts = 0.5*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0); 
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        f[l][m] = metricParts;

        metricParts = 0.5*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 1.0;
        h[l][m] = metricParts;

        metricParts = 0.5*exp(-(m-thetamax*(1.0/8.0))*(m-thetamax*(1.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(3.0/8.0))*(m-thetamax*(3.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(5.0/8.0))*(m-thetamax*(5.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts + 0.5*exp(-(m-thetamax*(7.0/8.0))*(m-thetamax*(7.0/8.0))*dtheta*dtheta*2.0);
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        g[l][m] = metricParts;
      }
    }
    break;
  case 27:         // eight gaussians (large)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 1.0*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*5.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        f[l][m] = metricParts;

        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 0.5*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 1.0;
        h[l][m] = metricParts;

        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 1.0*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*2.0*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        g[l][m] = metricParts;
      }
    }
    break;
  case 28:         // eight gaussians (small)
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 0.5*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        f[l][m] = metricParts;

        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 0.5*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 1.0;
        h[l][m] = metricParts;

        metricParts = 0.0;
        for (ntheta=1; ntheta<16; ntheta+=2) {
          metricParts = metricParts + 0.5*exp(-(m-ntheta*(thetamax/16.0))*(m-ntheta*(thetamax/16.0))*dtheta*dtheta*2.0);
        }
        metricParts = metricParts*0.5*exp(-(l-rmax*0.5)*(l-rmax*0.5)*dr*dr*2.0) + 2.0;
        g[l][m] = metricParts;
      }
    }
    break;
    default:    
    break;
  }

  for(l=0; l<rmax+1; l++) {
    for(m=0; m<thetamax+1; m++) {
      fn[l][m] = f[l][m];
      fp[l][m] = f[l][m];
      hn[l][m] = h[l][m];
      hp[l][m] = h[l][m];
      gn[l][m] = g[l][m];
      gp[l][m] = g[l][m];
    }
  }
  
  // initialize the velocity field using switches.
  switch (velocityProfile) {   
  case 0:      // null velocity field
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 0.0;
        v2[l][m] = 0.0;
      }
    }
    break;
  case 1:     // flat velocity field
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 1.0;
        v2[l][m] = 0.0;
      }
    }
    break;  
  case 2:     // Silk velocity profile
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        v1[l][m] = SilkVel(l);
        v2[l][m] = 0.0;
      }
    }
    break;
  case 3:     // gaussian peak
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 0.5*exp(-(l-rmax/2)*(l-rmax/2)*dr*dr*2.0);
        v2[l][m] = 0.0;
      }
    }
    break;
  case 8:   // 0.01*tanh centered, dr=0.1
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
      v1[l][m] = 0.01*(tanh(dr*(l-rmax/2)) - tanh(-rmax/2));
      v2[l][m] = 0.0;
     }
    }
    break;
  case 9:   // 0.01*tanh^20
    for(l=0; l<rmax+1; l++) {
      for(m=0; m<thetamax+1; m++) {
      v1[l][m] = 0.01*pow(tanh(dr*l),20);
      v2[l][m] = 0.0;
     }
    }
    break;
    case 10:   // logistic 
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 1.0/(1.0 + exp(-5.0*dr*(l-rmax/2)));
        v2[l][m] = 0.0;
      }
    }
    break;
    case 11:   // small sharp logistic 
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 0.01/(1.0 + exp(-5.0*dr*(l-rmax/2))) - 0.01/(1.0 + exp(-5.0*dr*(-rmax/2)));
        v2[l][m] = 0.0;
      }
    }
    break;
    case 12:   // small gradual logistic 
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 0.01/(1.0 + exp(-1.0*dr*(l-rmax/2))) - 0.01/(1.0 + exp(-1.0*dr*(-rmax/2)));
        v2[l][m] = 0.0;
      }
    }
    break;
    case 13:   // small sharpest logistic 
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        v1[l][m] = 0.01/(1.0 + exp(-10.0*dr*(l-rmax/2))) - 0.01/(1.0 + exp(-10.0*dr*(-rmax/2)));
        v2[l][m] = 0.0;
      }
    }
    break;
    case 14:   // small sharpest logistic with angular anisotropy
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        logisticPart = 0.01/(1.0 + exp(-10.0*dr*(l-rmax/2))) - 0.01/(1.0 + exp(-10.0*dr*(-rmax/2)));
        gaussianPart = 1.0*exp(-(m-thetamax/2.0)*(m-thetamax/2.0)*dtheta*dtheta*0.1);
        v1[l][m] = logisticPart*gaussianPart;
        v2[l][m] = 0.0;
      }
    }
    break;    
    case 15:   // small sharpest logistic with double angular anisotropy
      for(l=0; l<rmax+1; l++) {
        for(m=0; m<thetamax+1; m++) {
        logisticPart = 0.01/(1.0 + exp(-10.0*dr*(l-rmax/2))) - 0.01/(1.0 + exp(-10.0*dr*(-rmax/2)));
        gaussianPart = 1.0*exp(-(m-thetamax*0.25)*(m-thetamax*0.25)*dtheta*dtheta*0.1);
        gaussianPart = gaussianPart + 1.0*exp(-(m-thetamax*0.75)*(m-thetamax*0.75)*dtheta*dtheta*0.1);
        v1[l][m] = logisticPart*gaussianPart;
        v2[l][m] = 0.0;
      }
    }
    break;    
    default:
    break;
  }

  for(l=0; l<rmax+1; l++) {
    for(m=0; m<thetamax+1; m++) {
      v1n[l][m] = v1[l][m];
      v1p[l][m] = v1[l][m];
      v2n[l][m] = v2[l][m];
      v2p[l][m] = v2[l][m];
    }
  }
  
}

// finite differencing for spatial derivative operator
double GradR(double p[rmax+1][thetamax+1], int l, int m)
{
  if (l==0) return (p[l+1][m] - p[l][m])/dr; 
  else if (l==rmax) return (p[l][m] - p[l-1][m])/dr;
  else return (p[l+1][m] - p[l-1][m])/(2.0*dr);  //central differencing
}

// finite differencing for spatial derivative operator
double GradTheta(double p[rmax+1][thetamax+1], int l, int m)
{
  if (l==0) return 0.0;
  else if (m==0) return (p[l][m+1] - p[l][m])/dtheta; 
  else if (m==thetamax) return (p[l][m] - p[l][m-1])/dtheta;
  else return (p[l][m+1] - p[l][m-1])/(2.0*dtheta);  //central differencing
}

// laplacian function
double LapR(double p[rmax+1][thetamax+1], int l, int m)
{
  return (p[l+1][m] - 2.0*p[l][m] + p[l-1][m])/(dr*dr);
}

// laplacian function
double LapTheta(double p[rmax+1][thetamax+1], int l, int m)
{
  return (p[l][m+1] - 2.0*p[l][m] + p[l][m-1])/(dtheta*dtheta);
}

// mixed second derivative
double Mixed2ndDeriv(double p[rmax+1][thetamax+1], int l, int m)
{
  return (p[l+1][m+1] - p[l-1][m+1] - p[l+1][m-1] + p[l-1][m-1])/(4.0*dr*dtheta);
}

// Crank-Nicolson laplacian in radial direction
double CNR(double p[rmax+1][thetamax+1], double q[rmax+1][thetamax+1], int l, int m)
{
  double p1,p2;
  p1 = LapR(p,l,m);
  p2 = LapR(q,l,m);
  return 0.5*p1 + 0.5*p2;
}

// Crank-Nicolson in theta direction
double CNTheta(double p[rmax+1][thetamax+1], double q[rmax+1][thetamax+1], int l, int m)
{
  double p1,p2;
  p1 = LapTheta(p,l,m);
  p2 = LapTheta(q,l,m);
  return 0.5*p1 + 0.5*p2;
}

// Crank-Nicolson laplacian in radial direction divided by -det(metric)
double CNRdenom(double p[rmax+1][thetamax+1], double q[rmax+1][thetamax+1], int l, int m)
{
  double p1,p2,denom;

  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  p1 = LapR(p,l,m)/denom;

  denom = hp[l][m]*hp[l][m] - fp[l][m]*gp[l][m];
  p2 = LapR(q,l,m)/denom;

  return 0.5*p1 + 0.5*p2;
}

// Crank-Nicolson laplacian in radial direction divided by -det(metric)
double CNThetadenom(double p[rmax+1][thetamax+1], double q[rmax+1][thetamax+1], int l, int m)
{
  double p1,p2,denom;

  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  p1 = LapTheta(p,l,m)/denom;

  denom = hp[l][m]*hp[l][m] - fp[l][m]*gp[l][m];
  p2 = LapTheta(q,l,m)/denom;

  return 0.5*p1 + 0.5*p2;
}

// BoundaryR returns either a flat or const derivative boundary for radial direction
double BoundaryR(double p[rmax+1][thetamax+1], int l, int m, int boundarySwitch)
{
  double polar;
  int n;
  polar = 0.0;

  switch(boundarySwitch) {
  case 0:    // null derivative boundary
    if (l==0) return p[1][m];
    else return p[rmax-1][m];    
    break;
  case 1:    // constant derivative boundary
    if (l==0) return 2.0*p[1][m] - p[2][m];
    else return 2.0*p[rmax-1][m] - p[rmax-2][m];
    break;
  case 2:    // polar inner boundary
    if (l==0) {
      for (n=0; n<thetamax+1; n++) polar = polar + p[1][n];
      return polar/(101.0);
    }
    else return p[rmax-1][m];
    break;
  case 5:    // fixed at original value at origin, const derivative outer bndry
    if (l==0) return p[0][0]; 
    else return 2.0*p[rmax-1][m] - p[rmax-2][m];
    break;
 case 6:    // const 2nd deriv 
    if (l==0) return LapR(p,2,m) + 2.0*p[1][m] - p[2][m];
    else return LapR(p,rmax-2,m) + 2.0*p[rmax-1][m] - p[rmax-2][m];
    break;
 case 7:    // const 1st deriv 
    if (l==0) return p[2][m] - p[3][m] + p[1][m];
    else return p[rmax-2][m] - p[rmax-3][m] + p[rmax-1][m];
    break;
  default:
    break;
  }

  return 100.0;

}

// BoundaryTheta returns periodic boundary calculated by midpoint
double BoundaryTheta(double p[rmax+1][thetamax+1], int l, int m, int boundarySwitch)
{
  switch(boundarySwitch) {
  case 0: case 1: case 5: case 6: case 7:   
    return 0.5*(p[l][thetamax-1]+p[l][1]);
    break;
  default:
    break;
  }

  return 100.0;

}

// boundary conditions for dynamic variables
void DynamicVarBndries(int metricBoundary, int velBoundary)
{

  // 0 - null derivative boundary
  // 1 - constant derivative boundary
  // 5 - fixed value at origin (radial direction only)

  int l,m;
      
  for(m=1; m<thetamax; m++) {
    fn[0][m] = BoundaryR(fn,0,m,metricBoundary);
    fn[rmax][m] = BoundaryR(fn,rmax,m,metricBoundary);

    hn[0][m] = BoundaryR(hn,0,m,metricBoundary);
    hn[rmax][m] = BoundaryR(hn,rmax,m,metricBoundary);

    gn[0][m] = BoundaryR(gn,0,m,metricBoundary);
    gn[rmax][m] = BoundaryR(gn,rmax,m,metricBoundary);

    v1n[0][m] = BoundaryR(v1n,0,m,velBoundary);
    v1n[rmax][m] = BoundaryR(v1n,rmax,m,velBoundary);

    v2n[0][m] = BoundaryR(v2n,0,m,velBoundary);
    v2n[rmax][m] = BoundaryR(v2n,rmax,m,velBoundary);
  }

  for(l=1; l<rmax; l++) {
    fn[l][0] = BoundaryTheta(fn,l,0,metricBoundary);
    fn[l][thetamax] = BoundaryTheta(fn,l,thetamax,metricBoundary);

    hn[l][0] = BoundaryTheta(hn,l,0,metricBoundary);
    hn[l][thetamax] = BoundaryTheta(hn,l,thetamax,metricBoundary);

    gn[l][0] = BoundaryTheta(gn,l,0,metricBoundary);
    gn[l][thetamax] = BoundaryTheta(gn,l,thetamax,metricBoundary);

    v1n[l][0] = BoundaryTheta(v1n,l,0,velBoundary);
    v1n[l][thetamax] = BoundaryTheta(v1n,l,thetamax,velBoundary);

    v2n[l][0] = BoundaryTheta(v2n,l,0,velBoundary);
    v2n[l][thetamax] = BoundaryTheta(v2n,l,thetamax,velBoundary);
  }

  fn[0][0] = BoundaryR(fn,0,0,metricBoundary);
  fn[rmax][thetamax] = BoundaryTheta(fn,rmax,thetamax,metricBoundary);
  fn[0][thetamax] = BoundaryR(fn,0,0,metricBoundary);
  fn[rmax][0] = BoundaryTheta(fn,rmax,0,metricBoundary);

  hn[0][0] = BoundaryR(hn,0,0,metricBoundary);
  hn[rmax][thetamax] = BoundaryTheta(hn,rmax,thetamax,metricBoundary);
  hn[0][thetamax] = BoundaryR(hn,0,0,metricBoundary);
  hn[rmax][0] = BoundaryTheta(hn,rmax,0,metricBoundary);

  gn[0][0] = BoundaryR(gn,0,0,metricBoundary);
  gn[rmax][thetamax] = BoundaryTheta(gn,rmax,thetamax,metricBoundary);
  gn[0][thetamax] = BoundaryR(gn,0,0,metricBoundary);
  gn[rmax][0] = BoundaryTheta(gn,rmax,0,metricBoundary);

  v1n[0][0] = BoundaryR(v1n,0,0,velBoundary);
  v1n[rmax][thetamax] = BoundaryTheta(v1n,rmax,thetamax,velBoundary);
  v1n[0][thetamax] = BoundaryR(v1n,0,0,velBoundary);
  v1n[rmax][0] = BoundaryTheta(v1n,rmax,0,velBoundary);

  v2n[0][0] = BoundaryR(v2n,0,0,velBoundary);
  v2n[rmax][thetamax] = BoundaryTheta(v2n,rmax,thetamax,velBoundary);
  v2n[0][thetamax] = BoundaryR(v2n,0,0,velBoundary);
  v2n[rmax][0] = BoundaryTheta(v2n,rmax,0,velBoundary);
}

// condition on l=0 prevents divide by zero error
double Gamma111(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  if (l==0) return Gamma111(1,m);
  return ((2.0*h[l][m]*GradR(h,l,m) - GradR(f,l,m)*g[l][m])*l*dr + 2.0*h[l][m]*h[l][m] - GradTheta(f,l,m)*h[l][m])/(2.0*denom*l*dr);
}

double Gamma211(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  if (l==0) return Gamma211(1,m);
  return -((2.0*f[l][m]*GradR(h,l,m) - GradR(f,l,m)*h[l][m])*l*dr + 2.0*f[l][m]*h[l][m] - f[l][m]*GradTheta(f,l,m))/(2.0*denom*l*l*dr*dr);
}

double Gamma112(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  return (GradR(g,l,m)*h[l][m]*l*dr + 2.0*g[l][m]*h[l][m] - GradTheta(f,l,m)*g[l][m])/(2.0*denom);    
}

double Gamma212(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  if (l==0) return Gamma212(1,m);
  return -(f[l][m]*GradR(g,l,m)*l*dr - GradTheta(f,l,m)*h[l][m] + 2.0*f[l][m]*g[l][m])/(2.0*denom*l*dr);
}

double Gamma122(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  return (g[l][m]*GradR(g,l,m)*l*l*dr*dr + (-2.0*g[l][m]*GradTheta(h,l,m) + GradTheta(g,l,m)*h[l][m] + 2.0*g[l][m]*g[l][m])*l*dr)/(2.0*denom);
}

double Gamma222(int l, int m) {
  double denom;
  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
  return -(GradR(g,l,m)*h[l][m]*l*dr - 2.0*h[l][m]*GradTheta(h,l,m) + 2.0*g[l][m]*h[l][m] + f[l][m]*GradTheta(g,l,m))/(2.0*denom);
}

// returns negative of scalar curvature; calculation is based on the Riemann tensor of a 2D metric, of which the only
// unique component is R_1212
double ScalarCurvature(int l, int m) {
  double riemann,denom; 

  // if (l == 0) return 0.0;
  // else if (l == rmax) return 0.0;
  if (l == 0) return ScalarCurvature(1,m);
  else if (l == rmax) return ScalarCurvature(rmax-1,m);
  // else if (m == 0) return ScalarCurvature(l,1);
  // else if (m == thetamax) return ScalarCurvature(l,thetamax-1);
  else if (m == 0) return 0.5*(ScalarCurvature(l,1) + ScalarCurvature(l,thetamax-1));
  else if (m == thetamax) return 0.5*(ScalarCurvature(l,1) + ScalarCurvature(l,thetamax-1));
  else {
    denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];
    riemann = 0.0; 
    riemann = (2.0*GradR(g,l,m)*h[l][m]*GradR(h,l,m) - f[l][m]*GradR(g,l,m)*GradR(g,l,m) - GradR(f,l,m)*g[l][m]*GradR(g,l,m));
    riemann = riemann + (-2.0*h[l][m]*GradTheta(h,l,m) + 2.0*g[l][m]*h[l][m] + f[l][m]*GradTheta(g,l,m))*2.0*GradR(h,l,m)/(l*dr);
    riemann = riemann + (2.0*denom*Mixed2ndDeriv(h,l,m) + GradR(f,l,m)*g[l][m]*GradTheta(h,l,m) - 3.0*GradR(g,l,m)*h[l][m]*h[l][m])*2.0/(l*dr);
    riemann = riemann + (GradTheta(f,l,m)*GradR(g,l,m) - GradR(f,l,m)*GradTheta(g,l,m))*h[l][m]/(l*dr);
    riemann = riemann + (2.0*f[l][m]*g[l][m]*GradR(g,l,m) - GradR(f,l,m)*g[l][m]*g[l][m])*2.0/(l*dr); // **
    riemann = riemann + 2.0*(GradTheta(f,l,m)*h[l][m] - 2.0*f[l][m]*g[l][m])*GradTheta(h,l,m)/(l*l*dr*dr); 
    riemann = riemann + 2.0*(f[l][m]*GradTheta(g,l,m) + GradTheta(f,l,m)*g[l][m])*h[l][m]/(l*l*dr*dr);
    riemann = riemann - f[l][m]*GradTheta(f,l,m)*GradTheta(g,l,m)/(l*l*dr*dr) - GradTheta(f,l,m)*GradTheta(f,l,m)*g[l][m]/(l*l*dr*dr);
    riemann = riemann/(-4.0*denom*denom);
    riemann = riemann + CNRdenom(g,gp,l,m)/(2.0);  // **
    riemann = riemann + CNThetadenom(f,fp,l,m)/(2.0*l*l*dr*dr);
    return (-1.0)*riemann;
  }
}

double GrowthTensor11(int l, int m) {
  double gtParts;
  gtParts = 2.0*f[l][m]*(GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]);  
  gtParts = gtParts + 2.0*h[l][m]*l*dr*(GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]);
  return gtParts;
}

// divided by r
double GrowthTensor12(int l, int m) {
  double gtParts;
  if (l==0) return GrowthTensor12(1,m);
  gtParts = h[l][m]*(GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]);
  gtParts = gtParts + g[l][m]*l*dr*(GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]);
  gtParts = gtParts + f[l][m]*(GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m])/(l*dr);
  gtParts = gtParts + h[l][m]*(GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]);
  return gtParts;
}

// divided by r*r
double GrowthTensor22(int l, int m) {
  double gtParts;
  if (l==0) return GrowthTensor22(1,m);
  gtParts = 2.0*h[l][m]*(GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m])/(l*dr);
  gtParts = gtParts + 2.0*g[l][m]*(GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]);
  return gtParts;
}

// radial component of velocity field (v1)
double VelocityR(int l, int m, double curvatureCoupling, double advection, double velDif) {
  double vparts,denom;

  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];

  vparts = curvatureCoupling*(-gamma111[l][m]*v1[l][m]*v1[l][m] - 2.0*gamma112[l][m]*v1[l][m]*v2[l][m] - gamma122[l][m]*v2[l][m]*v2[l][m]);
  vparts = vparts - advection*v1[l][m]*(GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]);
  vparts = vparts - advection*v2[l][m]*(GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m]);

  vparts = vparts + velDif*(-g[l][m]/denom)*(LapR(v1,l,m) + GradR(gamma111,l,m)*v1[l][m] + GradR(gamma112,l,m)*v2[l][m]);
  vparts = vparts + velDif*(-g[l][m]/denom)*(gamma111[l][m]*GradR(v1,l,m) + 2.0*gamma112[l][m]*GradR(v2,l,m) - gamma211[l][m]*GradTheta(v1,l,m));
  vparts = vparts + velDif*(-g[l][m]/denom)*(gamma112[l][m]*gamma212[l][m] - gamma211[l][m]*gamma122[l][m])*v2[l][m];

  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(LapTheta(v1,l,m) + GradTheta(gamma112,l,m)*v1[l][m] + GradTheta(gamma122,l,m)*v2[l][m]);
  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(2.0*gamma112[l][m]*GradTheta(v1,l,m) + 2.0*gamma122[l][m]*GradTheta(v2,l,m) - gamma122[l][m]*GradR(v1,l,m) - gamma222[l][m]*GradTheta(v1,l,m));
  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(gamma112[l][m]*(gamma112[l][m] - gamma222[l][m]) + gamma122[l][m]*(gamma212[l][m] - gamma111[l][m]))*v1[l][m];

  return vparts;
}

// angular component of velocity field (v2)
double VelocityTheta(int l, int m, double curvatureCoupling, double advection, double velDif) {
  double vparts,denom;

  denom = h[l][m]*h[l][m] - f[l][m]*g[l][m];

  vparts = curvatureCoupling*(-gamma211[l][m]*v1[l][m]*v1[l][m] - 2.0*gamma212[l][m]*v1[l][m]*v2[l][m] - gamma222[l][m]*v2[l][m]*v2[l][m]);
  vparts = vparts - advection*v1[l][m]*(GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]);
  vparts = vparts - advection*v2[l][m]*(GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]);

  vparts = vparts + velDif*(-g[l][m]/denom)*(LapR(v2,l,m) + GradR(gamma211,l,m)*v1[l][m] + GradR(gamma212,l,m)*v2[l][m]);
  vparts = vparts + velDif*(-g[l][m]/denom)*(2.0*gamma211[l][m]*GradR(v1,l,m) + 2.0*gamma212[l][m]*GradR(v2,l,m) - gamma111[l][m]*GradR(v2,l,m) - gamma211[l][m]*GradTheta(v2,l,m));
  vparts = vparts + velDif*(-g[l][m]/denom)*(gamma211[l][m]*(gamma112[l][m] - gamma222[l][m]) + gamma212[l][m]*(gamma212[l][m] - gamma111[l][m]))*v2[l][m];

  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(LapTheta(v2,l,m) + GradTheta(gamma212,l,m)*v1[l][m] + GradTheta(gamma222,l,m)*v2[l][m]);
  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(2.0*gamma212[l][m]*GradTheta(v1,l,m) - gamma222[l][m]*GradTheta(v2,l,m) - gamma122[l][m]*GradR(v2,l,m));
  vparts = vparts + velDif*(-f[l][m]/(denom*l*l*dr*dr))*(gamma212[l][m]*gamma112[l][m] - gamma122[l][m]*gamma211[l][m])*v1[l][m];

  return vparts;
}

// derived quantities
void DerivedQuantities(int boundaryVar)
{
  // 0 - null derivative boundary
  // 1 - constant derivative boundary
  // 5 - fixed value at origin (radial direction only)

  int l,m;

  for(l=1; l<rmax; l++) {
    for(m=1; m<thetamax; m++) {
      gamma111[l][m] = Gamma111(l,m);
      gamma211[l][m] = Gamma211(l,m);
      gamma112[l][m] = Gamma112(l,m);    
      gamma212[l][m] = Gamma212(l,m);
      gamma122[l][m] = Gamma122(l,m);
      gamma222[l][m] = Gamma222(l,m);
      K[l][m] = ScalarCurvature(l,m);
      T11[l][m] = GrowthTensor11(l,m);
      T12[l][m] = GrowthTensor12(l,m);
      T22[l][m] = GrowthTensor22(l,m);
    }
  }

  for(m=1; m<thetamax; m++) { // R boundary
    l=0;
    gamma111[l][m] = Gamma111(l,m);
    gamma211[l][m] = Gamma211(l,m);
    gamma112[l][m] = Gamma112(l,m);    
    gamma212[l][m] = Gamma212(l,m);
    gamma122[l][m] = Gamma122(l,m);
    gamma222[l][m] = Gamma222(l,m);
    K[l][m] = ScalarCurvature(l,m);
    T11[l][m] = GrowthTensor11(l,m);
    T12[l][m] = GrowthTensor12(l,m);
    T22[l][m] = GrowthTensor22(l,m);

    l=rmax;
    gamma111[l][m] = Gamma111(l,m);
    gamma211[l][m] = Gamma211(l,m);
    gamma112[l][m] = Gamma112(l,m);    
    gamma212[l][m] = Gamma212(l,m);
    gamma122[l][m] = Gamma122(l,m);
    gamma222[l][m] = Gamma222(l,m);
    K[l][m] = ScalarCurvature(l,m);
    T11[l][m] = GrowthTensor11(l,m);
    T12[l][m] = GrowthTensor12(l,m);
    T22[l][m] = GrowthTensor22(l,m);
  }

  for(l=1; l<rmax; l++) { // Theta boundary
    m = 0;  
    gamma111[l][m] = Gamma111(l,m);
    gamma211[l][m] = Gamma211(l,m);
    gamma112[l][m] = Gamma112(l,m);    
    gamma212[l][m] = Gamma212(l,m);
    gamma122[l][m] = Gamma122(l,m);
    gamma222[l][m] = Gamma222(l,m);
    K[l][m] = ScalarCurvature(l,m);
    T11[l][m] = GrowthTensor11(l,m);
    T12[l][m] = GrowthTensor12(l,m);
    T22[l][m] = GrowthTensor22(l,m);

    m=thetamax;
    gamma111[l][m] = Gamma111(l,m);
    gamma211[l][m] = Gamma211(l,m);
    gamma112[l][m] = Gamma112(l,m);    
    gamma212[l][m] = Gamma212(l,m);
    gamma122[l][m] = Gamma122(l,m);
    gamma222[l][m] = Gamma222(l,m);
    K[l][m] = ScalarCurvature(l,m);
    T11[l][m] = GrowthTensor11(l,m);
    T12[l][m] = GrowthTensor12(l,m);
    T22[l][m] = GrowthTensor22(l,m);
  }

  gamma111[0][0] = Gamma111(0,0);
  gamma111[0][thetamax] = gamma111[0][0]; 
  gamma111[rmax][0] = Gamma111(rmax,0);
  gamma111[rmax][thetamax] = gamma111[rmax][0];

  gamma211[0][0] = Gamma211(0,0);
  gamma211[0][thetamax] = gamma211[0][0]; 
  gamma211[rmax][0] = Gamma211(rmax,0);
  gamma211[rmax][thetamax] = gamma211[rmax][0];

  gamma112[0][0] = Gamma112(0,0);
  gamma112[0][thetamax] = gamma112[0][0]; 
  gamma112[rmax][0] = Gamma112(rmax,0);
  gamma112[rmax][thetamax] = gamma112[rmax][0];

  gamma212[0][0] = Gamma212(0,0);
  gamma212[0][thetamax] = gamma212[0][0]; 
  gamma212[rmax][0] = Gamma212(rmax,0);
  gamma212[rmax][thetamax] = gamma212[rmax][0];

  gamma122[0][0] = Gamma122(0,0);
  gamma122[0][thetamax] = gamma122[0][0]; 
  gamma122[rmax][0] = Gamma122(rmax,0);
  gamma122[rmax][thetamax] = gamma122[rmax][0];

  gamma222[0][0] = Gamma222(0,0);
  gamma222[0][thetamax] = gamma222[0][0]; 
  gamma222[rmax][0] = Gamma222(rmax,0);
  gamma222[rmax][thetamax] = gamma222[rmax][0];

  K[0][0] = ScalarCurvature(0,0);
  K[0][thetamax] = K[0][0];  
  K[rmax][0] = ScalarCurvature(rmax,0);
  K[rmax][thetamax] = K[rmax][0];
 
  T11[0][0] = GrowthTensor11(0,0);
  T12[0][0] = GrowthTensor12(0,0);
  T22[0][0] = GrowthTensor22(0,0);

  T11[0][thetamax] = T11[0][0];
  T12[0][thetamax] = T12[0][0];
  T22[0][thetamax] = T22[0][0];

  T11[rmax][0] = GrowthTensor11(rmax,0);
  T12[rmax][0] = GrowthTensor12(rmax,0);
  T22[rmax][0] = GrowthTensor22(rmax,0);

  T11[rmax][thetamax] = T11[rmax][0];
  T12[rmax][thetamax] = T12[rmax][0];
  T22[rmax][thetamax] = T22[rmax][0];

}

double RadialDistance(int lmax, int m)
{
  int l;
  double partialLength;

  partialLength = 0.0;

  for (l=1;l<lmax+1;l++) {
    partialLength = partialLength + sqrt(f[l][m])*dr;
  }

  return partialLength;

}

double ThetaDistance(int l, int mmax)
{
  int m;
  double partialLength;

  partialLength = 0.0;

  for (m=1;m<mmax+1;m++) {
    partialLength = partialLength + sqrt(g[l][m])*dtheta*l*dr;
  }

  return partialLength;

}

// \nabla_1 v^1 = GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]
// \nabla_2 v^2 = GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]
// \nabla_1 v^2 = GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]
// \nabla_2 v^1 = GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m]

double ExpansionTensor(int l, int m) {
  double parts;
  parts = GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m] + GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m];  
  return parts;
}

double ShearTensor11(int l, int m) {
  double parts;
  parts = -f[l][m]*(GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]);  
  return parts;
}

double ShearTensor12(int l, int m) {
  double parts;
  parts = l*l*dr*dr*g[l][m]*(GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]);
  parts = parts + f[l][m]*(GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m]);
  parts = parts - l*dr*h[l][m]*ExpansionTensor(l,m);
  return 0.5*parts;
}

double ShearTensor22(int l, int m) {
  double parts;
  parts = -l*l*dr*dr*g[l][m]*(GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]);  
  return parts;
}

double RotationTensor12(int l, int m) {
  double parts;
  parts = l*dr*h[l][m]*(GradR(v1,l,m) + gamma111[l][m]*v1[l][m] + gamma112[l][m]*v2[l][m]);
  parts = parts - l*dr*h[l][m]*(GradTheta(v2,l,m) + gamma212[l][m]*v1[l][m] + gamma222[l][m]*v2[l][m]);
  parts = parts + l*l*dr*dr*g[l][m]*(GradR(v2,l,m) + gamma211[l][m]*v1[l][m] + gamma212[l][m]*v2[l][m]);
  parts = parts - f[l][m]*(GradTheta(v1,l,m) + gamma112[l][m]*v1[l][m] + gamma122[l][m]*v2[l][m]);
  return 0.5*parts;
}

double REGR(int l, int m, double kappa, double kappa1) {
  double parts;
  parts = kappa*K[l][m] + kappa1*ExpansionTensor(l,m);
  return parts;
}

void CollectData(int timestamp, int testnumber, double cCoupling, double gtCoupling)
{
  int l,m;
  int deltaR = 5;

  sprintf(filename, filetag, "pm",testnumber,"_thetaPi0_",timestamp,".txt");
  dynamicsData_thetaPi0 = fopen(filename,"w");
  sprintf(filename, filetag, "pm",testnumber,"_thetaPi4_",timestamp,".txt");
  dynamicsData_thetaPi4 = fopen(filename,"w");
  sprintf(filename, filetag, "pm",testnumber,"_thetaPi2_",timestamp,".txt");
  dynamicsData_thetaPi2 = fopen(filename,"w");
  sprintf(filename, filetag, "pm",testnumber,"_thetaPi_",timestamp,".txt");
  dynamicsData_thetaPi = fopen(filename,"w");  
  sprintf(filename,"%s%d%s","twoDdata",timestamp,".txt");
  twoDdata = fopen(filename,"w");
  sprintf(filename,"%s%d%s","modgrowthPi0_",timestamp,".txt");
  modgrowth_Pi0 = fopen(filename,"w");
  sprintf(filename,"%s%d%s","modgrowthPi4_",timestamp,".txt");
  modgrowth_Pi4 = fopen(filename,"w");
  sprintf(filename,"%s%d%s","modgrowthPi2_",timestamp,".txt");
  modgrowth_Pi2 = fopen(filename,"w");
  sprintf(filename,"%s%d%s","modgrowthPi_",timestamp,".txt");
  modgrowth_Pi = fopen(filename,"w");
  for(l=0; l<rmax+1; l++)
    {
      m = 0;
      fprintf(dynamicsData_thetaPi0,"%6d %6d %6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",timestamp,l,m,fn[l][m],hn[l][m],gn[l][m],v1n[l][m],v2n[l][m],-1.0*K[l][m],T11[l][m],T12[l][m],T22[l][m],ExpansionTensor(l,m),ShearTensor11(l,m),ShearTensor12(l,m),ShearTensor22(l,m),RotationTensor12(l,m),REGR(l,m,cCoupling,gtCoupling));
      if (l%deltaR == 0) fprintf(modgrowth_Pi0, "%6d %6e %6e %6e %6e %6e \n",timestamp,l*dr,RadialDistance(l,m),fn[l][m],hn[l][m],gn[l][m]);

      m = thetamax/8;
      fprintf(dynamicsData_thetaPi4,"%6d %6d %6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",timestamp,l,m,fn[l][m],hn[l][m],gn[l][m],v1n[l][m],v2n[l][m],-1.0*K[l][m],T11[l][m],T12[l][m],T22[l][m],ExpansionTensor(l,m),ShearTensor11(l,m),ShearTensor12(l,m),ShearTensor22(l,m),RotationTensor12(l,m),REGR(l,m,cCoupling,gtCoupling));
      if (l%deltaR == 0) fprintf(modgrowth_Pi4, "%6d %6e %6e %6e %6e %6e \n",timestamp,l*dr,RadialDistance(l,m),fn[l][m],hn[l][m],gn[l][m]);

      m = thetamax/4;
      fprintf(dynamicsData_thetaPi2,"%6d %6d %6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",timestamp,l,m,fn[l][m],hn[l][m],gn[l][m],v1n[l][m],v2n[l][m],-1.0*K[l][m],T11[l][m],T12[l][m],T22[l][m],ExpansionTensor(l,m),ShearTensor11(l,m),ShearTensor12(l,m),ShearTensor22(l,m),RotationTensor12(l,m),REGR(l,m,cCoupling,gtCoupling));
      if (l%deltaR == 0) fprintf(modgrowth_Pi2, "%6d %6e %6e %6e %6e %6e \n",timestamp,l*dr,RadialDistance(l,m),fn[l][m],hn[l][m],gn[l][m]);

      m = thetamax/2;
      fprintf(dynamicsData_thetaPi,"%6d %6d %6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",timestamp,l,m,fn[l][m],hn[l][m],gn[l][m],v1n[l][m],v2n[l][m],-1.0*K[l][m],T11[l][m],T12[l][m],T22[l][m],ExpansionTensor(l,m),ShearTensor11(l,m),ShearTensor12(l,m),ShearTensor22(l,m),RotationTensor12(l,m),REGR(l,m,cCoupling,gtCoupling));
      if (l%deltaR == 0) fprintf(modgrowth_Pi, "%6d %6e %6e %6e %6e %6e \n",timestamp,l*dr,RadialDistance(l,m),fn[l][m],hn[l][m],gn[l][m]);
    
      for(m=0; m<thetamax+1; m++) fprintf(twoDdata, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",l*(1.0/rmax),m*(2.0*M_PI/thetamax),RadialDistance(l,m),fn[l][m],hn[l][m],gn[l][m],v1n[l][m],v2n[l][m],-1.0*K[l][m],T11[l][m],T12[l][m],T22[l][m],ExpansionTensor(l,m),ShearTensor11(l,m),ShearTensor12(l,m),ShearTensor22(l,m),RotationTensor12(l,m),REGR(l,m,cCoupling,gtCoupling));
      fprintf(twoDdata, "%s\n", "");

    }  
  
}
/**********************************************************************************************/
/**********************************************************************************************/


// PlantModel(model #, curvature forces, v*Grad(v), v damping, v diffusion, ricci flow, v*v, growth tensor)
int PlantModel(int testnum, int initialMetric, int initialVelocity, int metricBndry, int velBndry, int bndryVarDerived, double c0,double c1,double c2,double c3,double k0,double k1,double k2)
{
  int i,j,t;
 
  sprintf(filename,"%s%d%s","plantLength",testnum,".txt");
  plantLength = fopen(filename,"w");

  // print values of coupling constants to the console.
  printf("%6f %6f %6f %6f %6f %6f %6f \n",c0,c1,c2,c3,k0,k1,k2);

  // wipe arrays
  for(i=0; i<rmax+1; i++) {
    for(j=0; j<thetamax+1; j++) {
      f[i][j] = 0.0;
      h[i][j] = 0.0;
      g[i][j] = 0.0;
      v1[i][j] = 0.0;
      v2[i][j] = 0.0;
      fn[i][j] = 0.0;
      hn[i][j] = 0.0;
      gn[i][j] = 0.0;
      v1n[i][j] = 0.0;
      v2n[i][j] = 0.0;
      fp[i][j] = 0.0;
      hp[i][j] = 0.0;
      gp[i][j] = 0.0;
      v1p[i][j] = 0.0;
      v2p[i][j] = 0.0;
      gamma111[i][j] = 0.0;
      gamma211[i][j] = 0.0;
      gamma112[i][j] = 0.0;
      gamma212[i][j] = 0.0;
      gamma122[i][j] = 0.0;
      gamma222[i][j] = 0.0;
      T11[i][j] = 0.0;
      T12[i][j] = 0.0;
      T22[i][j] = 0.0;
      K[i][j] = 0.0;
    }
  }

  // Initialize the metric and velocity profile based on an analytic function (t=0)
  t = 0;
  InitialConditions(initialMetric,initialVelocity);
  DynamicVarBndries(metricBndry,velBndry);

  // calculate derived quantities
  DerivedQuantities(bndryVarDerived);
  
  // write initial data to a file  
  CollectData(t,testnum,k0,k2);
  fprintf(plantLength,"%6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",t,RadialDistance(rmax,0), RadialDistance(rmax,thetamax/8), RadialDistance(rmax,thetamax/4), RadialDistance(rmax,thetamax/2), ThetaDistance(1,thetamax), ThetaDistance(rmax/4,thetamax), ThetaDistance(rmax/2,thetamax), ThetaDistance(rmax,thetamax)); 

  // calculate metrics at successive time steps
  for(t=1; t<tmax; t++) {

    for(i=1; i<rmax; i++) {
      for(j=1; j<thetamax; j++) {

        fn[i][j] = f[i][j] + dt*(k0*K[i][j]*f[i][j] + k2*T11[i][j]);  // Holy shit, calculus works!
        hn[i][j] = h[i][j] + dt*(k0*K[i][j]*h[i][j] + k2*T12[i][j]);
        gn[i][j] = g[i][j] + dt*(k0*K[i][j]*g[i][j] + k2*T22[i][j]);

        v1n[i][j] = v1[i][j] + dt*VelocityR(i,j,c0,c1,c3);
        v2n[i][j] = v2[i][j] + dt*VelocityTheta(i,j,c0,c1,c3);
	
      }
    }
  
      // boundary conditions for dynamical variables
      DynamicVarBndries(metricBndry,velBndry);
      
      // update the dynamical variables for next iteration
      for(i=0; i<rmax+1; i++) {
        for(j=0; j<thetamax+1; j++) {
      	  fp[i][j] = f[i][j];
      	  f[i][j] = fn[i][j];
  	      hp[i][j] = h[i][j];
      	  h[i][j] = hn[i][j];
      	  gp[i][j] = g[i][j];
      	  g[i][j] = gn[i][j];
      	  v1p[i][j] = v1[i][j];
      	  v1[i][j] = v1n[i][j];
      	  v2p[i][j] = v2[i][j];
      	  v2[i][j] = v2n[i][j];
      	}
      }

      // calculate derived quantities
      DerivedQuantities(bndryVarDerived);

      // radius, circumference
      if (t%10 == 0) {
        fprintf(plantLength,"%6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",t,RadialDistance(rmax,0), RadialDistance(rmax,thetamax/8), RadialDistance(rmax,thetamax/4), RadialDistance(rmax,thetamax/2), ThetaDistance(1,thetamax), ThetaDistance(rmax/4,thetamax), ThetaDistance(rmax/2,thetamax), ThetaDistance(rmax,thetamax)); 
      }

      switch (t) {
      case 10: case 50: case 100: case 500: case 1000: case 5000: case 10000: case 50000: case 100000: case 500000: case 750000:
        CollectData(t,testnum,k0,k2);
        break;
      default:
        break;
      }

  }

  return 2;
}

/**********************************************************************************************/
/**********************************************************************************************/

// PlantModel(model #, initial metric, initial velocity, dynamic vars boundary condition, derived vars boundary condition, curvature forces, v*Grad(v), v damping, v diffusion, ricci flow, v*v, growth tensor)
int main(int argc, char *argv[])
{

  int modelnumber;
  double testParam;
  int modelmax = 1;

  /* parameters that indicate the intial metric and velocity functions. Passed to InitialConditions function.
   metricTest: 0 - flat 
               1 - flat metric with gaussian perturbation in g11
               2 - parabolic g11
               3 - silk type g11
               4 - centered gaussian g11
               5 - gaussian g11, parabolic g22
               6 - logistic g11
               7 - sharper gaussian g11
               8 - sharpest gaussian g11
               9 - sharper gaussian g11, shifted right
               10 - sharper gaussian g11, shifted left
               11 - gaussian g11
               12 - gaussian g22
               13 - gaussian g12
               14 - gaussian in (r,theta) g11
               15 - gaussian in (r,theta) g22
               16 - gaussian in (r) g11 and g22, small
               17 - gaussian in (r) g11 and g22, large
               18 - gaussian in (r,theta) g11 and g22
               19 - tiny single gaussian in (r,theta) g11, g12 and g22
               20 - small single gaussian in (r,theta) g11, g12 and g22
               21 - large single gaussian in (r,theta) g11, g12 and g22
               22 - double cresent
               23 - small double gaussian in g11, g12, g22
               24 - large double gaussian in g11, g12, g22
               25 - large 4-gaussian in g11, g12, g22
               26 - small 4-gaussian in g11, g12, g22
               27 - large 8-gaussian in g11, g12, g22
               28 - small 8-gaussian in g11, g12, g22
   velocityTest: 0 - null
                 1 - flat v1
                 2 - Silk velocity profile, v1
                 3 - gaussian v1
                 8 - 0.01*tanh centered, dr=0.1
                 9 - 0.01*tanh^20 off center
                 10 - logistic
                 11 - small sharp logistic, zero @inner boundary (used in most tests, based on Silk+Erickson data)
                 12 - small gradual logistic, zero @inner boundary
                 13 - small sharpest logistic, zero @inner boundary
   bndry[Metric/Vel/Derived]Var:  0 - null derivative
                     1 - constant derivative
                     2 - polar inner boundary
                     5 - fixed to initial value at inner boundary, const derivative at outer bndry
                     6 - constant 2nd derivative
                     7 - constant 1st derivative using 3 points
  */                


  int testCase = atoi(argv[1]);
  int metricTest = atoi(argv[2]);
  int velocityTest = atoi(argv[3]);
  // int bndryMetricVar = 7;
  int bndryMetricVar = 0;
  int bndryVelVar = 5;
  int bndryDerivedVar = 0;

  for (modelnumber=0;modelnumber<modelmax+1;modelnumber+=2)
    {
      testParam = (double) modelnumber;

      switch (testCase) {

      case 0: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,0.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // null test
      break; 
      case 1: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,20.0,0.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // vel diff 
      break;
      case 2: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,5.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // Ricci flow
      break;
      case 3: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,-100.0,-100.0,0.0,0.0,0.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // advection 
      break;
      case 4: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,0.0,0.0,1.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // Growth Tensor
      break;
      case 5: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.008,5.0,0.0,3.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 6: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,10.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // Ricci Flow
      break;
      case 7: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,-5.0,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // Ricci Flow
      break;
      case 8: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.02,3.0,0.0,3.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 9: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,1.0,0.0,1.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // RF and GT
      break;
      case 10: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,5.0,0.0,1.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // RF and GT
      break;
      case 11: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,1.0,0.0,10.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // RF and GT
      break;
      case 12: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,1.0,0.0,(testParam+1.0)) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   //
      break;
      case 13: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.1,6.0,0.0,6.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 14: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2*(testParam+1.0),0.0,0.1) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // negative RF
      break;
      case 15: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2*(testParam+1.0),0.0,-0.1) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // negative RF and GT
      break;
      case 16: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2,0.0,-0.1*(testParam+1.0)) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // negative RF and GT
      break;
      case 17: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5,0.0,5.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 18: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.75,0.0,0.75) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 19: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.005,0.75,0.0,0.75) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 20: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,0.75,0.0,0.75) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // most tested paramter set
      break;
      case 21: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5,0.0,3.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 22: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5,0.0,(testParam+1.0)) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 24: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,0.5,0.0,0.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 25: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,0.75,0.0,0.75) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   // 
      break;
      case 26: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.75,0.75,0.0,0.75) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   //
      break;
      case 27: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.005,1.0,0.0,1.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   //
      break;
      case 28: if (PlantModel(modelnumber,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.005,2.0,0.0,2.0) == 2) printf("%s %d %s","model",modelnumber,"ran successfully \n");   //
      break;
      default:
      break;

      }
      
    }

}

/*************************************************************************************/
