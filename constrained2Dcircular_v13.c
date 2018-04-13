/*
Covariant metric structure:

g11 = f[i]
g12 = g21 = h[i]*r = h[i]*i*du
g22 = g[i]*r*r = g[i]*i*i*du*du

Contravariant metric structure:
g11upper = -g[i]/den
g12upper = g21upper = h[i]/(den*r) = h[i]/(den*i*du)
g22upper = -f[i]/(den*r*r) = -f[i]/(den*i*i*du*du)

   where den = h[i]*h[i] - f[i]*g[i]



January 2014

For intial conformally flat metric and null velocity:
f = 1; h = 0; g = 1; v1 = 0; v2 = 0;
- no dynamics even for very long times (tmax = 10,000)
- when v1_initial = 1, instabilities occur in growth tensor 
and velocity diffusion tests. Could be inner boundary condition.


Feb
- Trying various combinations of null vs. constant derivative at boundaries
(for dynamical variables only).
- center weighted differencing appears to work better than forward differencing
- radial velocity drops at origin, affects growth tensor (or vice versa??)

JustVelDiffusion folder: [intial metric]_[dynamic var boundary]_[derived val boundary]
JustRicciFlow folder: [intial metric]_[dynamic var boundary]_[derived val boundary]

March 24 2014
Long-time simulations (10k time steps) have an instability near the origin.
- lowering the amplitude of the tanh function helps a bit
- constant derivative boundary conditions for the metric do not help
- best combination seems to be a fixed inner boundary for metric and velocity,
  and a constant derivative inner bndry for the derived values like GT, connection
  coefficients, etc.

April 24: (RF15 tests)
- for large velocities, the ficticious boundary model outperforms the true boundary model
  in any one test case, and in RF15 to 10k time steps
- for big vel, the ficticious and true boundary simulations fail in similar ways at ~25k time steps.
- for small velocities, the true boundary model is more stable for an initially flat metric
  and the centered gaussian metric

April 25: Simple tests
True Bndry, small vel, centered gaussian:
- RF has non-zero velocity values at origin (10^-24)
- vel diff develops inner bndry instability for ~0.2, fails around 1000t
- GT develops inner bndry instability around 1.0, fails around 1700t even for null coupling
- Same for Lagrangian coupling
- RF16-2 breaks around velDiff=0.2 - inner boundary instability in Velocity field.
- value of vel at inner boundary is several orders of magnitude higher in true boundary 
  model than in ficticious boundary 

May 2: long time sims, centered guass, small vel
- RF16-1 fails at ~550k for highest vel diffusion coupling
- RF17 fails at ~125k for lowest RF coupling
- RF18 fails at ~100k for highest GT coupling

May 7: implemented direct calculation of derived values at the inner boundary
- allowed simulations to run to 10^6 time steps

May 8: knock-out tests @dt=10^-4
- GT coupling: lower bound is 0; upper bound ~2.0 (test vel diff -3)
- RF coupling: lower bound is ~0.1; upper bound ~1.0 (test RF-2; test zero GT, vel diff)
- vel diffusion: lower bound is 0; upper bound ~0.01 (test GT-2)

December 2015:
- Dufort-Frankel method on second derivs of f and g in curvature and velocity eqns 
not as effective as in 1D.
- Crank-Nicolson second deriv of f for calculating curvature stabilizes results.

Feb 2015
- use null boundary condition on metric when starting flat; 3-point (case7) condition 
when starting non-flat.
 */

/******************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************/

const int imax = 100;        //max radial coordinate
const int thetaMax = 100;        //max radial coordinate
//const int tmax = 1001;
//const int tmax = 10001;
// const int tmax = 100001;
// const int tmax = 750001;
// const int tmax = 1000001;
const int tmax = 2000001;
// const int tmax = 4000001;
// const int tmax = 10000001;
//const int tmax = 100000001;
// const double dt = 0.000401;     // 
// const double dt = 0.000201;     // 
const double dt = 0.000101;     // regular dt (10^-4)
// const double dt = 0.0000501;    //small dt (5*10^-5)
//const double dt = 0.00000101;  //extra small dt (10^-6)
//const double dt = 0.000000501;  //xxsmall dt (5*10^-7)    
// const double du = 0.0503;       
const double du = 0.103;       

double fn[imax+1];       //metric component being calculated (t+1) - covariant 
double f[imax+1];        //current metric component (t) - covariant
double fp[imax+1];       //past metric (t-1) - covariant
double hn[imax+1];       //metric component being calculated (t+1) - covariant 
double h[imax+1];        //current metric component (t) - covariant
double hp[imax+1];       //past metric (t-1) - covariant
double gn[imax+1];       //metric component being calculated - covariant
double g[imax+1];        //current metric component - covariant
double gp[imax+1];       //past metric (t-1) - covariant
double v1n[imax+1];        //new radial velocity - contravariant
double v1[imax+1];         //current radial velocity - contravariant
double v1p[imax+1];        //past radial velocity - contravariant
double v2n[imax+1];        //new theta velocity - contravariant
double v2[imax+1];         //current theta velocity - contravariant
double v2p[imax+1];        //past theta velocity - contravariant
double gamma111[imax+1];  //connection coefficients
double gamma211[imax+1];  //connection coefficients
double gamma112[imax+1];  //connection coefficients
double gamma212[imax+1];  //connection coefficients
double gamma122[imax+1];  //connection coefficients
double gamma222[imax+1];  //connection coefficients
double T11[imax+1];        //Growth tensor component 
double T12[imax+1];        //Growth tensor component 
double T22[imax+1];        //Growth tensor component 
double K[imax+1];          //Scalar curvature

int i,t,theta;
double length;

FILE *oneDdata,*modLengths,*plantLength,*info,*twoDdata;
char filename[30];
char filetag[10] = "%s%d%s%d%s";

/******************************************************************/
/******************************************************************/

// The InitialConditions function exploits the fact that the metric and 
// velocity fields are global variables. The switch is set in main(), in
// the function call to PlantModel.
void InitialConditions(int metricFunction, int velocityProfile)
{
  // initialize the metric function using switches.
  switch (metricFunction) {
  case 0:        // flat metric
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0;
      h[i] = 0.0;  
      g[i] = 1.0;
    }
    break;
  case 1:         // flat metric with gaussian perturbation
    for(i=0; i<imax+1; i++) {
      f[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du/2.5) + 1.0; 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;    
  case 2:         // parabolic 
    for(i=0; i<imax+1; i++) {
      f[i] = du*du*i*i + 1.0; 
      h[i] = 0.0;  
      g[i] = 1.0;
    }
    break;  
  case 4:         // centered gaussian 
    for(i=0; i<imax+1; i++) {
      f[i] = 0.5*exp(-i*i*du*du/2.5) + 1.0; 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break; 
  case 5:         // flat metric with gaussian perturbation
    for(i=0; i<imax+1; i++) {
      f[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 1.0; 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;      
  case 12:         // flat metric with radial gaussian perturbation
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0; 
      h[i] = 0.0; 
      g[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 1.0; 
    }
    break; 
  case 17:         // flat metric with gaussian perturbation
    for(i=0; i<imax+1; i++) {
      f[i] = 0.2*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*4.0) + 1.0; 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;       
  case 18:         // limit of many symmetric gaussians (small)
    for(i=0; i<imax+1; i++) {
      f[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 2.0; 
      h[i] = 0.0; 
      g[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 2.0; 
    }
    break;       
  case 19:         // limit of many symmetric gaussians (large)
    for(i=0; i<imax+1; i++) {
      f[i] = 5.0*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*3.0) + 2.0; 
      h[i] = 0.0; 
      g[i] = 2.0*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*3.0) + 2.0; 
    }
    break;       
  case 20:         // initially + curvature
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (2.0-1.0)/(1+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;       
  case 21:         // initially + curvature
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (1.1-1.0)/(1+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break; 
  case 22:         // aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (1.25-1.0)/(1+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;        
  case 23:         // constant curvature simple
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0/(1-0.004*i*i*du*du); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 24:         // curvature mimicing numerical results
    for(i=0; i<imax+1; i++) {
      f[i] = 0.5*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 0.9; // + 0.1/(1.0 + exp(-5.0*du*(i-imax/2.0))); 
      h[i] = 0.0; 
      g[i] = 0.5 + 0.1/(1.0 + exp(-5.0*du*(i-imax/2.0))) - 0.1/(1.0 + exp(-5.0*du*(-imax/2.0))); 
    }
    break;  
  case 25:         // curvature mimicing numerical results
    for(i=0; i<imax+1; i++) {
      f[i] = 0.9*exp(-(i-imax/2.0)*(i-imax/2.0)*du*du*2.0) + 0.9 + (1.25-1.0)/(1+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 0.8 + 0.1/(1.0 + exp(-5.0*du*(i-imax/2.0))); 
    }
    break;  
  case 26:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 27:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + 0.02*i*i*du*du*exp(-i*du) + 0.25/(1.0+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break; 
 case 28:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + 0.25/(1.0+exp(-0.75*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;          
 case 29:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + 0.25/(1.0+exp(-1.2*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;           
 case 30:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + 0.25/(1.0+5.0*exp(-1.0*(i*du - 4.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break; 
  case 31:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.02*i*du)/(1.0+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;      
  case 32:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*du)/(1.0+exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 33:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.05*exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 34:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+5.0*exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 35:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.0001*exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;              
  case 36:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.00001*exp(-1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 37:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+exp(1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;              
  case 38:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.00001*exp(1.0*(i*du - 5.0))); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;             
  case 39:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.00001*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                             
  case 40:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.002*i*i*du*du)/(1.0+0.0001*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;          
  case 41:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.01*i*i*du*du)/(1.0+0.00007*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                   
  case 42:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.01*i*i*du*du)/(1.0+0.0001*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                  
  case 43:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.01*i*i*du*du)/(1.0+0.00002*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                 
  case 44:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.01*i*i*du*du)/(1.0+0.00001*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;  
  case 45:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.004*i*i*du*du)/(1.0+0.00002*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                   
  case 46:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.006*i*i*du*du)/(1.0+0.00004*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;
  case 47:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.008*i*i*du*du)/(1.0+0.00006*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break; 
  case 48:         // modified aceta thesis results
    for(i=0; i<imax+1; i++) {
      f[i] = 1.0 + (0.003*i*i*du*du)/(1.0+0.00001*pow(i*du,4.0)); 
      h[i] = 0.0; 
      g[i] = 1.0; 
    }
    break;                                                                                                             
    default:
    break;
  }

  for(i=0; i<imax+1; i++) {
    fn[i] = f[i];
    fp[i] = f[i];
    hn[i] = h[i];
    hp[i] = h[i];
    gn[i] = g[i];
    gp[i] = g[i];
  }
  
  // initialize the velocity field using switches.
  switch (velocityProfile) {   
  case 0:      // null velocity field
    for(i=0; i<imax+1; i++) {
      v1[i] = 0.0;
      v2[i] = 0.0;
    }
    break;
  case 1:     // flat velocity field
    for(i=0; i<imax+1; i++) {
      v1[i] = 1.0;
      v2[i] = 0.0;
    }
    break;  
  case 3:   // tanh^2
    for(i=0; i<imax+1; i++) {
    v1[i] = tanh(du*i)*tanh(du*i);
    v2[i] = 0.0;
  }
  break;
  case 4:   // tanh^4
    for(i=0; i<imax+1; i++) {
    v1[i] = pow(tanh(du*i),4);
    v2[i] = 0.0;
  }
  break;
  case 5:   // tanh^6
    for(i=0; i<imax+1; i++) {
    v1[i] = pow(tanh(du*i),6);
    v2[i] = 0.0;
  }
  break;
  case 6:   // tanh^8
    for(i=0; i<imax+1; i++) {
    v1[i] = pow(tanh(du*i),8);
    v2[i] = 0.0;
  }
  break;
  case 7:   // tanh^12
    for(i=0; i<imax+1; i++) {
    v1[i] = pow(tanh(du*i),12);
    v2[i] = 0.0;
  }
  break;
  case 8:   // tanh^20
    for(i=0; i<imax+1; i++) {
    v1[i] = pow(tanh(du*i),20);
    v2[i] = 0.0;
  }
  break;
  case 9:   // 0.01*tanh^20
    for(i=0; i<imax+1; i++) {
    v1[i] = 0.01*pow(tanh(du*i),20);
    v2[i] = 0.0;
  }
  break;
  case 10:   // logistic 
    for(i=0; i<imax+1; i++) {
    v1[i] = 1.0/(1.0 + exp(-5.0*du*(i-imax/2.0)));
    v2[i] = 0.0;
  }
  break;
  case 11:   // small sharp logistic; v(0) = 0
    for(i=0; i<imax+1; i++) {
      v1[i] = 0.01/(1.0 + exp(-5.0*du*(i-imax/2.0))) - 0.01/(1.0 + exp(-5.0*du*(-imax/2.0)));
      v2[i] = 0.0;
  }
  break;
  case 12:   // small sharp logistic; v(0) = 0
    for(i=0; i<imax+1; i++) {
      v1[i] = 0.01/(1.0 + exp(-5.0*du*(i-imax*0.75))) - 0.01/(1.0 + exp(-5.0*du*(-imax*0.75)));
      v2[i] = 0.0;
  }
  break;  
  case 13:   // small sharp logistic; v(0) = 0
    for(i=0; i<imax+1; i++) {
      v1[i] = 0.01/(1.0 + exp(-5.0*du*(i-imax*0.25))) - 0.01/(1.0 + exp(-5.0*du*(-imax*0.25)));
      v2[i] = 0.0;
  }
  break;
  default:
    break;
  }

  for(i=0; i<imax+1; i++) {
    v1n[i] = v1[i];
    v1p[i] = v1[i];
    v2n[i] = v2[i];
    v2p[i] = v2[i];
  }
  
}

// finite differencing for spatial derivative operator
double Grad(double p[imax], int i)
{
  if (i==0) return (p[i+1] - p[i])/du; 
  else if (i==imax) return (p[i] - p[i-1])/du;
  else return (p[i+1] - p[i-1])/(2.0*du);  //central differencing
}

// laplacian function
double Lap(double p[imax], int i)
{
  return (p[i+1] - 2.0*p[i] + p[i-1])/(du*du);
}

// Laplacian operator multiplied by 1/f
double Lapf(double p[imax], int i)
{
  double p1,p2;
  p1 = (f[i] + f[i-1])/2.0;
  p2 = (f[i] + f[i+1])/2.0;
  return ( (p[i-1] - p[i])/p1 + (p[i+1] - p[i])/p2 )/(du*du);
}

// Laplacian operator multiplied by 1/g
double Lapg(double p[imax], int i)
{
  double p1,p2;
  p1 = (g[i] + g[i-1])/2.0;
  p2 = (g[i] + g[i+1])/2.0;
  return ( (p[i-1] - p[i])/p1 + (p[i+1] - p[i])/p2 )/(du*du);
}

// Laplacian operator multiplied by 1/f/g
double Lapfg(double p[imax], int i)
{
  double p1,p2;
  p1 = (g[i] + g[i-1])/2.0;
  p1 = p1*(f[i] + f[i-1])/2.0;
  p2 = (g[i] + g[i+1])/2.0;
  p2 = p2*(f[i] + f[i+1])/2.0;
  return ( (p[i-1] - p[i])/p1 + (p[i+1] - p[i])/p2 )/(du*du);
}


// Laplacian operator multiplied by 1/f/f
double Lapf2(double p[imax], int i)
{
  double p1,p2;
  p1 = (f[i] + f[i-1])/2.0;
  p2 = (f[i] + f[i+1])/2.0;
  return ( (p[i-1] - p[i])/p1/p1 + (p[i+1] - p[i])/p2/p2 )/(du*du);
}

// Dufort-Frankel laplacian on f(r) (tested in GradR(gamma111))
double DFf(int i)
{
  double frr;

  frr = (f[i+1] - fp[i] - fn[i] + f[i-1])/(du*du); 
  return 0.5*(frr/f[i] - Grad(f,i)*Grad(f,i)/f[i]/f[i]);
}

// Dufort-Frankel laplacian on g(r)
double DFg(int i)
{
  return (g[i+1] - gp[i] - gn[i] + g[i-1])/(du*du); 
}

// Crank-Nicolson laplacian on f(r); equivalent to GradR(gamma111,i)
double CNf(int i)
{
  double p1,p2,frr;
  p1 = Lap(f,i)/f[i];
  p2 = Lap(fp,i)/fp[i];
  frr = 0.5*p1 + 0.5*p2;
  return 0.5*(frr - Grad(f,i)*Grad(f,i)/f[i]/f[i]);
}

// Crank-Nicolson laplacian/g on g(r)
double CNg(int i)
{
  double p1,p2;
  p1 = Lap(g,i)/g[i];
  p2 = Lap(gp,i)/gp[i];
  return 0.5*p1 + 0.5*p2;
}

// Crank-Nicolson laplacian on g(r)
double CN(int i)
{
  double p1,p2;
  p1 = Lap(g,i);
  p2 = Lap(gp,i);
  return 0.5*p1 + 0.5*p2;
}

// BoundaryType returns either a flat or const derivative boundary
double BoundaryType(double p[imax], int i, int boundarySwitch)
{
  switch(boundarySwitch) {
  case 0:    // null derivative boundary
    return p[i];
    break;
  case 1:    // constant derivative boundary
    if (i == 1) return 2.0*p[i] - p[i+1];
    else return 2.0*p[i] - p[i-1];
    break;
  case 5:    // fixed at original value at origin, const derivative outer bndry
    if (i==1) return p[0];
    else return 2.0*p[i] - p[i-1];;
    break;
  case 6:    // fixed at original value at origin, null derivative
    if (i==1) {
      p[1] = p[0];
      return p[0];
    }
    else return 2.0*p[i] - p[i-1];;
    break;
  case 7:    // constant 1st deriv using 3 points
    if (i==1) return p[2] - p[3] + p[1];
    else return p[imax-2] - p[imax-3] + p[imax-1];
    break;
  default:
    break;
  }

  return 100.0;

}

// boundary conditions for dynamic variables
void DynamicVarBndries(int metricBoundary, int velBoundary)
{
 
  double nonLinTerm,den; 
  
  fn[0] = BoundaryType(fn,1,metricBoundary);
  fn[imax] = BoundaryType(fn,imax-1,metricBoundary);

  hn[0] = BoundaryType(hn,1,metricBoundary);
  hn[imax] = BoundaryType(hn,imax-1,metricBoundary);

  gn[0] = BoundaryType(gn,1,metricBoundary);
  gn[imax] = BoundaryType(gn,imax-1,metricBoundary);

  v1n[0] = BoundaryType(v1n,1,velBoundary);
  v1n[imax] = BoundaryType(v1n,imax-1,velBoundary);

  v2n[0] = BoundaryType(v2n,1,velBoundary);
  v2n[imax] = BoundaryType(v2n,imax-1,velBoundary);
}

// derived quantities
void DerivedQuantities(int boundaryVar)
{
  // 0 - null derivative boundary
  // 1 - constant derivative boundary

  double denom;
  length = 0.0;

  for(i=1; i<imax; i++)
    {
      denom = -f[i]*g[i];

      gamma111[i] = Grad(f,i)/(2.0*f[i]);
      gamma211[i] = 0.0;
      gamma112[i] = 0.0;
      gamma212[i] = (Grad(g,i)*i*du + 2.0*g[i])/(2.0*g[i]*i*du);
      gamma122[i] = -(Grad(g,i)*i*i*du*du + 2.0*g[i]*i*du)/(2.0*f[i]);
      gamma222[i] = 0.0;
      
      // calculates negative of scalar curvature = -R1212/det(g)
      K[i] = 0.5*CNg(i)/f[i] - (f[i]*Grad(g,i)*Grad(g,i) + Grad(f,i)*g[i]*Grad(g,i))/(4.0*denom*denom) + (4.0*f[i]*g[i]*Grad(g,i) - 2.0*Grad(f,i)*g[i]*g[i])/(4.0*denom*denom*i*du);

      T11[i] = 2.0*f[i]*(Grad(v1,i) + gamma111[i]*v1[i]);
      T12[i] = 0.0;
      T22[i] = 2.0*g[i]*(gamma212[i]*v1[i]); //divided by r*r

      length = length + sqrt(f[i])*du;
    }

  denom = -f[0]*g[0];  

  gamma111[0] = Grad(f,0)/(2.0*f[0]);
  gamma111[imax] = Grad(f,imax)/(2.0*f[imax]);

  gamma211[0] = 0.0;
  gamma211[imax] = 0.0;

  gamma112[0] = 0.0;
  gamma112[imax] = 0.0;

  gamma212[0] = BoundaryType(gamma212,1,boundaryVar);
  gamma212[imax] = Grad(g,imax)/(2.0*g[imax]) + 1.0/imax/du;

  gamma122[0] = 0.0;
  gamma122[imax] =  -(imax*du/f[imax])*(Grad(g,imax)*imax*du*0.5 + g[imax]);

  gamma222[0] = 0.0;
  gamma222[imax] = 0.0;

  K[0] = BoundaryType(K,1,boundaryVar);
  K[imax] = BoundaryType(K,imax-1,boundaryVar);

  // T11[0] = 2.0*f[0]*(Grad(v1,0) + gamma111[0]*v1[0]);
  // T12[0] = 0.0;
  // T22[0] = 2.0*g[0]*gamma212[0]*v1[0];

  T11[0] = T11[1];
  T12[0] = 0.0;
  T22[0] = T22[1];

  T11[imax] = 2.0*f[imax]*(Grad(v1,imax) + gamma111[imax]*v1[imax]);
  T12[imax] = 0.0;
  T22[imax] = 2.0*g[imax]*gamma212[imax]*v1[imax];

  length = length + sqrt(f[imax])*du;
}

// metric evolution equation without g_rr term
double NonLinMetric(int l, double RFcoupling, double GTcoupling)
{
  double denom, nonLinTerm;
  denom = -f[l]*g[l];
  nonLinTerm = -(f[l]*Grad(g,l)*Grad(g,l) + Grad(f,l)*g[l]*Grad(g,l))/(4.0*denom*denom) + (4.0*f[l]*g[l]*Grad(g,l) - 2.0*Grad(f,l)*g[l]*g[l])/(4.0*denom*denom*l*du);
  nonLinTerm = RFcoupling*g[l]*nonLinTerm;
  nonLinTerm = nonLinTerm + GTcoupling*T22[i]; // T22 already divided by r*r
  return nonLinTerm;
}

// velocity equation without v1_rr term
double NonLinVel(int l, double ccCoupling, double advection, double vdiff)
{
  double vparts;
  vparts = ccCoupling*(-gamma111[l]*v1[l]*v1[l]); //using contravariant velocity
  vparts = vparts - advection*v1[l]*(Grad(v1,l) + gamma111[l]*v1[l]);
  vparts = vparts + (vdiff/f[l])*(Grad(gamma111,l)*v1[l] + gamma111[l]*Grad(v1,l));
  vparts = vparts + (vdiff/(g[l]*l*l*du*du))*(-gamma122[l]*Grad(v1,l) + gamma122[l]*(gamma212[l] - gamma111[l])*v1[l]);
  return vparts;
}

double LocalLength(int lmax)
{
  int l;
  double partialLength;

  partialLength = 0.0;

  for (l=1;l<lmax+1;l++) {
    partialLength = partialLength + sqrt(f[l])*du;
  }

  return partialLength;

}

void CollectData(int timestamp, int modelnumber)
{
  int l,di;
  di = 5;

  sprintf(filename,filetag, "pm",modelnumber,"_",timestamp,".txt");
  oneDdata = fopen(filename,"w");
  sprintf(filename,"%s%d%s", "twoDdata",timestamp,".txt");
  twoDdata = fopen(filename,"w");
  sprintf(filename,"%s%d%s", "modLength",timestamp,".txt");
  modLengths = fopen(filename,"w");
  for(l=0; l<imax+1; l++) {
    fprintf(oneDdata,"%6d %6d %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e \n",timestamp,l,fn[l],hn[l],gn[l],v1n[l],v2n[l],-1.0*K[l],T11[l],T12[l],T22[l],gamma111[l],gamma211[l],gamma112[l],gamma212[l],gamma122[l],gamma222[l],LocalLength(l));
    for(theta=0;theta<thetaMax+1;theta++) fprintf(twoDdata, "%e %e %e %e %e %e \n",l*(1.0/imax),theta*(2.0*M_PI/thetaMax),LocalLength(l),-1.0*K[l],sqrt(fn[l]),sqrt(gn[l])); 
    fprintf(twoDdata, "%s\n", "");
    if (l%di == 0) fprintf(modLengths, "%6d %6f %6f %6f %6f %6f \n",timestamp,l*du,LocalLength(l),fn[l],hn[l],gn[l]);
  }

}
/**********************************************************************************************/
/**********************************************************************************************/


// PlantModel(model #, curvature forces, v*Grad(v), v damping, v diffusion, ricci flow, v*v, growth tensor)
int PlantModel(int testnum, int initialMetric, int initialVelocity, int metricBndry, int velBndry, int bndryVarDerived, double c0,double c1,double c2,double c3,double k0,double k1,double k2)
{

  double alpha,finalLength,den;
 
  sprintf(filename,"%s%d%s","plantLength",testnum,".txt");
  plantLength = fopen(filename,"w");


  // print values of coupling constants to the console.
  printf("%6f %6f %6f %6f %6f %6f %6f \n",c0,c1,c2,c3,k0,k1,k2);

  // wipe arrays
  for(i=0; i<imax+1; i++)
    {
      f[i] = 0.0;
      h[i] = 0.0;
      g[i] = 0.0;
      v1[i] = 0.0;
      v2[i] = 0.0;
      fn[i] = 0.0;
      hn[i] = 0.0;
      gn[i] = 0.0;
      v1n[i] = 0.0;
      v2n[i] = 0.0;
      fp[i] = 0.0;
      hp[i] = 0.0;
      gp[i] = 0.0;
      v1p[i] = 0.0;
      v2p[i] = 0.0;
      gamma111[i] = 0.0;
      gamma211[i] = 0.0;
      gamma112[i] = 0.0;
      gamma212[i] = 0.0;
      gamma122[i] = 0.0;
      gamma222[i] = 0.0;
      T11[i] = 0.0;
      T12[i] = 0.0;
      T22[i] = 0.0;
      K[i] = 0.0;
    }

  // Initialize the metric and velocity profile based on analytic functions (t=0)
  t = 0;
  InitialConditions(initialMetric,initialVelocity);

  // calculate derived quantities
  DerivedQuantities(bndryVarDerived);
  
  // write initial data to a file
  CollectData(t,testnum);
  fprintf(plantLength,"%6d %.9e %.9e %.9e \n",t,length,2.0*M_PI*(imax/2)*du*sqrt(gn[imax/2]),2.0*M_PI*imax*du*sqrt(gn[imax])); 
  
  t = 1;
  for(i=1; i<imax; i++)
    {
      fn[i] = f[i] + dt*(k0*K[i]*f[i] + k2*T11[i]);
      hn[i] = 0.0;
      gn[i] = g[i] + dt*(k0*K[i]*g[i] + k2*T22[i]);

      den = -f[i]*g[i];

      v1n[i] = v1[i] + dt*(NonLinVel(i,c0,c1,c3) + (c3/f[i])*Lap(v1,i));

    }

  // boundary conditions for dynamic variables
  DynamicVarBndries(metricBndry,velBndry);
  
  // update dynamical variable arrays
  for(i=0; i<imax+1; i++)
    {
      fp[i] = f[i];
      f[i] = fn[i];
      hp[i] = h[i];
      h[i] = hn[i];
      gp[i] = g[i];
      g[i] = gn[i];
      v1p[i] = v1[i];
      v1[i] = v1n[i];
      v2p[i] = v2[i];
      v2[i] = v2n[i];
    }

  // calculate derived quantities
  DerivedQuantities(bndryVarDerived);
  
  // calculate metrics at successive time steps using Dufort-Frankel method
  // ...Holy shit! Calculus works!
  for(t=2; t<tmax; t++)
  {
    for(i=1; i<imax; i++)
      {
        den = -f[i]*g[i];

        fn[i] = f[i] + dt*(k0*K[i]*f[i] + k2*T11[i]); 
        hn[i] = 0.0;   

	      alpha = (k0/f[i])*(dt/du/du);
        gn[i] = gp[i]*(1.0-alpha)/(1.0+alpha) + (g[i+1]+g[i-1])*alpha/(1.0+alpha) + 2.0*dt*NonLinMetric(i,k0,k2)/(1.0+alpha);

	      alpha = (c3/f[i])*(2.0*dt/du/du);	  
	      v1n[i] = v1p[i]*(1.0-alpha)/(1.0+alpha) + (v1[i+1]+v1[i-1])*alpha/(1.0+alpha) + 2.0*dt*NonLinVel(i,c0,c1,c3)/(1.0+alpha);	
      }
  
      // boundary conditions for dynamical variables
      DynamicVarBndries(metricBndry,velBndry);

      // update the dynamical variables for next iteration
      for(i=0; i<imax+1; i++) 
    	{
    	  fp[i] = f[i];
    	  f[i] = fn[i];
	      hp[i] = h[i];
    	  h[i] = hn[i];
    	  gp[i] = g[i];
    	  g[i] = gn[i];
    	  v1p[i] = v1[i];
    	  v1[i] = v1n[i];
    	  v2p[i] = v2[i];
    	  v2[i] = v2n[i];
    	}
      
      // calculate derived quantities
      DerivedQuantities(bndryVarDerived);

      if (t%1000 == 0) {
	      fprintf(plantLength,"%6d %.9e %.9e %.9e \n",t,length,2.0*M_PI*(imax/2)*du*sqrt(gn[imax/2]),2.0*M_PI*imax*du*sqrt(gn[imax])); 
      }
      
      switch (t) {
      case 1000: case 10000: case 50000: case 100000: case 500000: case 750000: case 1000000: case 1500000: case 2000000: case 2500000: case 3000000: case 3500000: case 4000000:
        CollectData(t,testnum);
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

  int i;
  double di;
  int imax = 1;
//  int imax = 4;

  /* parameters that indicate the intiail metric and velocity functions. Passed to InitialConditions function.
   metricTest: 0 - flat 
               1 - flat metric with gaussian perturbation
               2 - parabolic
               4 - centered gaussian 
               5 - flat metric with gaussian perturbation in g11
               12 - flat metric with radial gaussian perturbation in g22
               17 - flat metric with gaussian perturbation in g11 (small)
               18 - limit of many symmetric gaussians (small)
               19 - limit of many symmetric gaussians (large)
               22 - aceta thesis results
               23 - initially constant curvature
   velocityTest: 0 - null
                 1 - flat
                 3 - tanh^2
                 4 - tanh^4
                 5 - tanh^6
                 6 - tanh^8
                 7 - tanh^12
                 8 - tanh^20
                 9 - 0.01*tanh^20
                 10 - logistic
                 11 - small sharp logistic; v(0) = 0 (most widely tested, based on Silk+Erickson data)
   bndry[Metric/Vel/Derived]Var:  0 - null derivative
                     1 - constant derivative
                     5 - fixed to initial value at inner boundary, const derivative at outer bndry
                     6 - fixed to IV at inner bndry and null first derivative
                     7 - constant 1st derivative using 3 points
  */                

  int testCase = atoi(argv[1]);     // simulation paramaters can be set on the console or in a script
  int metricTest = atoi(argv[2]);
  int velocityTest = atoi(argv[3]);
  int bndryMetricVar = 0; 
  // int bndryMetricVar = 7; 
  int bndryVelVar = 5;
  int bndryDerivedVar = 0;

  for (i=0;i<imax+1;i+=2)
    {
      di = (double) i;

      switch (testCase) {
 
      case 1: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.075,1.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero vel diff -1
      break;
      case 2: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.075,1.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero vel diff -2
      break;
      case 3: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.07,1.0,0.0,2.5) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero vel diff -3 
      break;
      case 4: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.075,2.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -1
      break;
      case 5: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.008,5.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      case 6: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.007,1.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero GT -1
      break;
      case 7: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.005,1.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero GT -2 (max-ish value of vel diffusion)
      break;
      case 8: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.02,3.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero GT -3
      break;
      case 9: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.03,3.0,0.0,5.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero GT, vel diff (max-ish RF ~1.0)
      break;
      case 10: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,5.0,0.0,7.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // #16-3
      break;
      case 11: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,7.0,0.0,7.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // #17-1
      break;
      case 12: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,7.0,0.0,5.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // #18-1
      break;
      case 13: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.1,6.0,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // negative GT
      break;
      case 14: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2*(di+1.0),0.0,0.1) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // negative RF
      break;
      case 15: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2*(di+1.0),0.0,-0.1) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // negative RF and GT
      break;
      case 16: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,-0.2,0.0,-0.1*(di+1.0)) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // negative RF and GT
      break;
      case 17: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5,0.0,5.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 18: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.75,0.0,0.75) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 19: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5*(di+1.0),0.0,0.1) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max RF
      break;
      case 20: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,0.75,0.0,0.75) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 21: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,0.5,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 22: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,1.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 23: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.001,1.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // v8
      break;
      case 24: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,0.0,0.0,0.0,0.0,0.5,0.0,0.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 25: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,1.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 26: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.1,1.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 27: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.05,1.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 28: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.1,1.0,0.0,5.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 29: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.1,1.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // max GT
      break;
      case 30: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,5.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      case 31: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,4.0,0.0,2.5) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      case 32: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.015,2.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      case 33: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.015,3.0,0.0,3.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      case 34: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,0.5,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 35: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,0.1,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 36: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.015,0.1,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // all max (same as case 12)
      break;
      case 37: if (PlantModel(i,metricTest,velocityTest,bndryMetricVar,bndryVelVar,bndryDerivedVar,1.0,1.0,0.0,0.01,1.5,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n");   // zero RF -2 (minimal value of RF coupling for 1M time steps)
      break;
      default:
      break;

      }
      
    }

}

/*************************************************************************************/
