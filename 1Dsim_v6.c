// 1D plant growth simulation
// author: Julia Pulwicki, Feb 2016



/******************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************/

const int imax = 1000;
const int jmax = 1;
const int tmax = 100001;
const int sectionStep = 10;
const int timestep = 1000;
//const double dt = 0.00501;   
//const double dt = 0.000101;    
const double dt = 0.0000101;    
//const double dt = 0.00000101;    
const double du = 0.0503;       
const double dv = 0.0503;

double hp[imax+1][1];    //future metric function
double ho[imax+1][1];    //present metric function
double hn[imax+1][1];    //past metric function
double vp[imax+1][1];    //future velocity field
double vo[imax+1][1];    //present velocity field
double vn[imax+1][1];    //past velocity field
double  R[imax+1][1];    //Ricci tensor 
double GT[imax+1][1];    //growth tensor 
double  s[imax+1][1];    //modular lengths
double  g[imax+1][1];    //first time derivative of scale factor
double cc[imax+1][1];    //connection coefficients
double kpz[imax+1][1];   //kpz in dufort frankel
double v1[imax+1][1];    //first spatial deriv of vel
double v2[imax+1][1];    //second spatial deriv of vel

int i,j,t;
double A,b;

FILE *plantGrowthData,*sectionalGrowthData,*plantLength,*info;
char filename[50];
char filetag[50] = "%s%d%s%d%s";

// function h: the analytic function used to initialize the metric
// and set boundary conditions.
// Select one by uncommenting
double h(double dblx)
{
  //return 5.0*exp(-(dblx)*(dblx)*du*du/500.0) + 1.0;
  //return (1.4*cos(dblx*du*0.1)*cos(dblx*du*0.1) + 1.0)*exp(-dblx/500.0);
  //return 0.5;
  //return du*dblx/50.3 + 1.0; //linear
  return 1.0;
  //return pow(0.3*(1.0 + exp(0.02*dblx-5.0)),1.0/15.0);
  //return 0.5*tanh(0.007*dblx - 5.0) + 1.0;
  //return exp(-0.1*du*du*(dblx-800.0)*(dblx-800.0)) + 1.0;
  //return exp(-0.1*du*du*(dblx-500.0)*(dblx-500.0)) + 1.0;
  //return 1.0/(1.0 + exp(-0.2*du*(dblx-500.00))) + 1.0;


  // double q;
  // q = 1.0 + exp(0.01*dblx-5.0);
  // //return 0.2*(1.0/q + 1.0);
  // //return 0.45*(1.0/q + 1.0);
  // //return 0.6*(1.0/q + 1.0);
  // return 0.9*(1.0/q + 1.0);

}

// function v: the analytic function used to initialize the velocity field
// Select one by uncommenting
double v(double dblx)
{
  //return 0.0005;
  //return 0.0;
  //return 1.0;
 
  //return pow(tanh(0.05*du*dblx),20);
  //return 0.5*tanh(du*(dblx - imax*0.8)) + 0.5;
  //return tanh(du*(dblx - imax*0.8)) + 1.0;
  //return tanh(du*(0.2*dblx - imax*0.2)) + 1.0;
  //return 0.45*tanh(du*(dblx - imax*0.8)) + 0.45; //steep velocity
  //return 0.25*tanh(du*0.25*(dblx - imax*0.8)) + 0.25; //shallow, small velocity
  //return 0.45*tanh(du*0.4*(dblx - imax*0.7)) + 0.45; // silk vel 
  //return 1.0*tanh(du*0.4*(dblx - imax*0.8)) + 1.0; // silk vel large
  //return 0.5*tanh(du*0.4*(dblx - imax*0.7)) + 0.5; // silk vel med
  //return 0.5*tanh(du*(dblx - imax*0.8)) + 0.5; // steep vel, big
  //return 1.0*tanh(du*(dblx - imax*0.8)) + 1.0; // steep vel, max 2

  //return (atan(du*(dblx - imax*0.8)) + atan(du*imax*0.8))/atan(du*999.0)/2.0; 
  //return (atan(0.5*du*(dblx - imax*0.8)) + atan(0.5*du*imax*0.8))/atan(0.5*du*999.0)/2.0; 
  //return (atan(0.2*du*(dblx - imax*0.8)) + atan(0.2*du*imax*0.8))/atan(0.2*du*999.0)/2.0; 

  // logistic functions 
  //return 1.0/(1.0 + exp(-1.0*du*(dblx-imax*0.8)));
  //return 1.0/(1.0 + exp(-0.5*du*(dblx-imax*0.8)));
  //return 1.0/(1.0 + exp(-0.2*du*(dblx-imax*0.8)));
  //return 1.0/(1.0 + exp(-1.0*du*(dblx-imax*0.5)));
  //return 1.0/(1.0 + exp(-0.5*du*(dblx-imax*0.5)));
  return 1.0/(1.0 + exp(-0.2*du*(dblx-imax*0.5)));

  //return du*dblx/50.3; //linear

  // // shallow stepped linear
  // if (dblx<imax*0.2) return 0.0;
  // else if (dblx<800.0) return 1.5*du*(dblx-imax*0.2)/45.27;
  // else return 1.0;

  // // shallower stepped linear
  // if (dblx<imax*0.2) return 0.0;
  // else if (dblx<800.0) return 0.75*du*(dblx-imax*0.2)/45.27;
  // else return 0.5;  

  // // stepped linear
  // if (dblx<500.0) return 0.0;
  // else if (dblx<800.0) return 3.0*du*(dblx-imax*0.5)/45.27;
  // else return 1.0;


/*
  // Silk velocity profile (decreasing):
  double q,n;
  n = 5.0;
  q = 1.0 + exp(0.02*dblx-5.0);
  //return 0.9*pow(q,-1.0/n);
  return 0.9*pow(q,-1.0/n) + 0.6;
  */
}

// finite differencing for spatial derivative operator
double grad(double p[imax][jmax], int i, int j)
{
  if (i==0) return (p[i+1][j] - p[i][j])/du; 
  else if (i==imax-1) return (p[i][j] - p[i-1][j])/du;
  return (p[i+1][j] - p[i-1][j])/(2.0*du);  //central differencing
  //else return (p[i+1][j] - p[i][j])/du;            //upwind
  //return (p[i][j] - p[i-1][j])/du;          //downwind
}

// laplacian operator 
double lap(double p[imax][jmax], int i, int j)
{
  return (p[i+1][j] - 2.0*p[i][j] + p[i-1][j])/(du*du);  
}

// laplacian operator multiplied by 1/f 
double lapf(double p[imax][jmax], int i, int j)
{
  double p1,p2;
  p1 = (ho[i][j] + ho[i-1][j])/2.0;
  p2 = (ho[i][j] + ho[i+1][j])/2.0;
  return ( (p[i-1][j] - p[i][j])/p1 + (p[i+1][j] - p[i][j])/p2 )/(du*du);
  //return (p[i+1][j] - 2.0*p[i][j] + p[i-1][j])/(du*du);
}

// laplacian operator multiplied by 1/f/f
double lapf2(double p[imax][jmax], int i, int j)
{
  double p1,p2;
  p1 = (ho[i][j] + ho[i-1][j])/2.0;
  p2 = (ho[i][j] + ho[i+1][j])/2.0;
  return ( (p[i-1][j] - p[i][j])/p1/p1 + (p[i+1][j] - p[i][j])/p2/p2 )/(du*du);
  //return (p[i+1][j] - 2.0*p[i][j] + p[i-1][j])/(du*du);
}

// crank-nicholson type laplacian/f/f
double crankNichols(int i, int j)
{
  double p1,p2;
  p1 = lap(ho,i,j)/ho[i][j]/ho[i][j];
  p2 = lap(hp,i,j)/hp[i][j]/hp[i][j];
  return 0.5*p1 + 0.5*p2;
}

// dufort-frankel type laplacian/f/f
double dufortFrankel(int i, int j)
{
  return (ho[i+1][j] - hp[i][j] - hn[i][j] + ho[i-1][j])/(ho[i][j]*ho[i][j]*du*du);
}


/**********************************************************************************************/

// plantModel(model #, curvature forces, v*grad(v), v damping, v diffusion, ricci flow, v*v, growth tensor)
int plantModel(int testnum, double c,double c1,double c2,double c3,double k,double k1,double k2)
{

  double alpha,hDiff,hRest,vDiff,vRest,length,h0,h1,h2,h3,h4;

  info = fopen("info.txt", "w");
  fprintf(info, "%6f %6f %6f %6f %6f %6f %6f \n",c,c1,c2,c3,k,k1,k2);
  sprintf(filename,"%s%d%s","plantLength",testnum,".txt");
  plantLength = fopen(filename,"w");
  plantGrowthData = fopen("plantGrowthData.txt","w");

  // Initialize the metric and velocity profile based on an analytic function (t=0)

  t = 0;
  j = 0;
  for(i=0; i<imax; i++)
    // for(j=0; j<jmax; j++)
      {
      	hn[i][j] = h(i);
      	ho[i][j] = h(i);
        hp[i][j] = h(i);
      	vn[i][j] = v(i);
      	vo[i][j] = v(i);
        vp[i][j] = v(i);
      }

// calculate derived quantities
  for(i=1; i<imax-1; i++)
    {
      // first derivative of scale factor
      R[i][j] = grad(ho,i,j);
      // Growth tensor
      GT[i][j] = (2.0*ho[i][j]*grad(vo,i,j) + grad(ho,i,j)*vo[i][j]);
      // first time derivative of scale factor
      g[i][j] = (sqrt(hp[i][j]) - sqrt(hn[i][j]))/(2.0*dt);
      // connection coefficient
      cc[i][j] = grad(ho,i,j)*vo[i][j]*vo[i][j]/(2.0*ho[i][j]);
      // KPZ in Dufort Frankel
      kpz[i][j] = -0.5*k*(grad(ho,i,j)*grad(ho,i,j)/ho[i][j]/ho[i][j] - dufortFrankel(i,j)*ho[i][j]);
      // first vel deriv
      v1[i][j] = grad(vo,i,j);
      // second vel deriv
      v2[i][j] = lap(vo,i,j);
    }

    R[0][0] = grad(ho,0,j);
    R[imax-1][0] = grad(ho,imax-1,j);
    GT[0][0] = (2.0*ho[0][j]*grad(vo,0,j) + grad(ho,0,j)*vo[0][j]);
    GT[imax-1][0] = (2.0*ho[imax-1][j]*grad(vo,imax-1,j) + grad(ho,imax-1,j)*vo[imax-1][j]);
    cc[0][0] = grad(ho,0,j)/(2.0*ho[0][j]);
    cc[imax-1][0] = grad(ho,imax-1,j)/(2.0*ho[imax-1][j]);
    kpz[0][0] = kpz[1][0];
    kpz[imax-1][0] = kpz[imax-2][0];
    v1[0][0] = grad(vo,0,j);
    v1[imax-1][0] = grad(vo,imax-1,j);
    v2[0][0] = v2[1][0];
    v2[imax-1][0] = v2[imax-2][0];

  // length of leaf    
  length = 0.0;  
  for(i=0; i<imax; i++)
    {
      length = length + sqrt(ho[i][j])*du;
      s[i][j] = length;   //length of plant up to n*di*du (di defined as a global constant)
    }
  fprintf(plantLength,"%6d %6f \n",t,length); 

  sprintf(filename,filetag, "modLength",testnum,"_",t,".txt");
  sectionalGrowthData = fopen(filename,"w");
  for(i=0; i<imax; i++)
    {
      if (i%sectionStep == 0) fprintf(plantGrowthData,"%6d %6d %6d %e %e %e %e %e %e %e %e %e %e \n",t,i,j,ho[i][j],vo[i][j],R[i][j],GT[i][j],g[i][j],cc[i][j],kpz[i][j],v1[i][j],v2[i][j],s[i][j]);
      if (i%sectionStep == 0) fprintf(sectionalGrowthData, "%6d %6d %6f \n",t,i,s[i][j]);
    }
  fprintf(plantGrowthData, "%s\n", "");

  // calculate the metric for the next time step, go (t=1), using simple finite differencing.
  t = 1;
  for(i=1; i<imax-1; i++)
      {
      	hRest = 0.5*k*(lap(hn,i,j)/hn[i][j] - grad(hn,i,j)*grad(hn,i,j)/hn[i][j]/hn[i][j]);
      	hRest = hRest + k1*vn[i][j]*vn[i][j]*hn[i][j]*hn[i][j];
      	hRest = hRest + k2*(2.0*hn[i][j]*grad(vn,i,j) + grad(hn,i,j)*vn[i][j]);
      	ho[i][j] = hn[i][j] + dt*hRest;

      	vRest = -c*grad(hn,i,j)*vn[i][j]*vn[i][j]/(2.0*hn[i][j]);
      	vRest = vRest - c1*vn[i][j]*(grad(vn,i,j) + grad(hn,i,j)*vn[i][j]/(2.0*hn[i][j]));
      	vRest = vRest - c2*vn[i][j];
      	vRest = vRest + c3*lap(vn,i,j)/hn[i][j]; 
        vRest = vRest + c3*0.5*lap(hn,i,j)*vn[i][j]/hn[i][j]/hn[i][j];
        vRest = vRest + c3*0.5*grad(hn,i,j)*grad(vn,i,j)/hn[i][j]/hn[i][j];
        vRest = vRest - c3*0.5*grad(hn,i,j)*grad(hn,i,j)*vn[i][j]/hn[i][j]/hn[i][j]/hn[i][j];
      	vo[i][j] = vn[i][j] + dt*vRest;

      }

  // set boundary conditions
  // null boundary
  ho[0][0] = ho[1][0]; 
  ho[imax-1][0] = ho[imax-2][0];
  // constant inner boundary
  //ho[0][0] = h(0);
  // null first derivative boundary
  //ho[0][0] = 2.0*ho[1][0] - ho[2][0];
  //ho[imax-1][0] = 2.0*ho[imax-2][0] - ho[imax-3][0];
  // null boundary
  //vo[0][0] = vo[1][0];
  //vo[imax-1][0] = vo[imax-2][0];
  // constant inner boundary, constant first deriv outer boundary
  vo[0][0] = v(0);
  vo[imax-1][0] = 2.0*vo[imax-2][0] - vo[imax-3][0];

  // calculate the Ricci tensor, growth tensor
  for(i=1; i<imax-1; i++)
    {
      // first derivative of scale factor
      R[i][j] = grad(ho,i,j);
      // Growth tensor
      GT[i][j] = (2.0*ho[i][j]*grad(vo,i,j) + grad(ho,i,j)*vo[i][j]);
      // first time derivative of scale factor
      g[i][j] = (sqrt(hp[i][j]) - sqrt(hn[i][j]))/(2.0*dt);
      // connection coefficients
      cc[i][j] = grad(ho,i,j)/(2.0*ho[i][j]);
      // KPZ in Dufort Frankel
      kpz[i][j] = -0.5*k*(grad(ho,i,j)*grad(ho,i,j)/ho[i][j]/ho[i][j] - dufortFrankel(i,j)*ho[i][j]);
      // first vel deriv
      v1[i][j] = grad(vo,i,j);
      // second vel deriv 
      v2[i][j] = lap(vo,i,j);
    }
  
    R[0][0] = grad(ho,0,j);
    R[imax-1][0] = grad(ho,imax-1,j);
    GT[0][0] = (2.0*ho[0][j]*grad(vo,0,j) + grad(ho,0,j)*vo[0][j]);
    GT[imax-1][0] = (2.0*ho[imax-1][j]*grad(vo,imax-1,j) + grad(ho,imax-1,j)*vo[imax-1][j]);
    cc[0][0] = grad(ho,0,j)/(2.0*ho[0][j]);
    cc[imax-1][0] = grad(ho,imax-1,j)/(2.0*ho[imax-1][j]);
    kpz[0][0] = kpz[1][0];
    kpz[imax-1][0] = kpz[imax-2][0];
    v1[0][0] = grad(vo,0,j);
    v1[imax-1][0] = grad(vo,imax-1,j);
    v2[0][0] = v2[1][0];
    v2[imax-1][0] = v2[imax-2][0];



  // calculate metrics at successive time steps using Dufort-Frankel method
  for(t=2; t<tmax; t++)
    {
      for(i=1; i<imax-1; i++)
    	{   
    	  alpha = k*dt/(ho[i][j]*du*du);
    	  hDiff = hn[i][j]*(1.0 - alpha)/(1.0 + alpha) + (ho[i+1][j] + ho[i-1][j])*alpha/(1.0 + alpha);
    	  hRest = -k*grad(ho,i,j)*grad(ho,i,j)/2.0/ho[i][j]/ho[i][j];
    	  hRest = hRest + k1*ho[i][j]*ho[i][j]*vo[i][j]*vo[i][j];
    	  hRest = hRest + k2*(2.0*ho[i][j]*grad(vo,i,j) + grad(ho,i,j)*vo[i][j]);
    	  hRest = hRest*2.0*dt/(1.0 + alpha);
    	  hp[i][j] = hDiff + hRest;
    	  
    	  alpha = c3*2.0*dt/(ho[i][j]*du*du);
    	  vDiff = vn[i][j]*(1.0 - alpha)/(1.0 + alpha) + (vo[i+1][j] + vo[i-1][j])*alpha/(1.0 + alpha);
    	  vRest = -c*grad(ho,i,j)*vo[i][j]*vo[i][j]/(2.0*ho[i][j]);
    	  vRest = vRest - c1*vo[i][j]*(grad(vo,i,j) + grad(ho,i,j)*vo[i][j]/(2.0*ho[i][j]));
    	  vRest = vRest - c2*vo[i][j];
    	  h2 = ho[i][j]*ho[i][j];
    	  h3 = h2*ho[i][j];
    	  vRest = vRest + c3*0.5*dufortFrankel(i,j)*vo[i][j] + c3*0.5*grad(vo,i,j)*grad(ho,i,j)/h2 - c3*0.5*grad(ho,i,j)*grad(ho,i,j)*vo[i][j]/h3;
    	  vRest = vRest*2.0*dt/(1.0 + alpha);
    	  vp[i][j] = vDiff + vRest;
    	}
    	
      // set boundary conditions
      // null boundary
      hp[0][0] = hp[1][0];
      hp[imax-1][0] = hp[imax-2][0];
      // constant inner boundary
      //hp[0][0] = h(i);
      // const first derivative boundary
      //hp[0][0] = 2.0*hp[1][0] - hp[2][0];
      //hp[imax-1][0] = 2.0*hp[imax-2][0] - hp[imax-3][0];
      // constant inner boundary
      vp[0][0] = v(0);
      // null boundary
      //vp[imax-1][0] = vp[imax-2][0];
      //vp[0][0] = vp[1][0];
      // const first derivative boundary
      //vp[0][0] = 2.0*vp[1][0] - vp[2][0];
      vp[imax-1][0] = 2.0*vp[imax-2][0] - vp[imax-3][0];
            
      // calculate derived quantities
      g[0][j] = 0.0;
      for(i=1; i<imax-1; i++)
      	{
          // first derivative of scale factor
      	  R[i][j] = grad(ho,i,j);
      	  // Growth tensor
      	  GT[i][j] = (2.0*ho[i][j]*grad(vo,i,j) + grad(ho,i,j)*vo[i][j]);
      	  // first time derivative of scale factor
      	  g[i][j] = (sqrt(hp[i][j]) - sqrt(hn[i][j]))/(2.0*dt);
      	  // connection coefficients
      	  cc[i][j] = grad(ho,i,j)/(2.0*ho[i][j]);
      	  // KPZ in Dufort Frankel
          kpz[i][j] = 0.5*k*(dufortFrankel(i,j)*ho[i][j] - grad(ho,i,j)*grad(ho,i,j)/ho[i][j]/ho[i][j]);
      	  // first vel deriv
      	  v1[i][j] = grad(vo,i,j);
      	  // second vel deriv
      	  v2[i][j] = (vo[i+1][j] - vn[i][j] - vp[i][j] + vo[i-1][j])/du/du;
      	}
      R[0][0] = grad(ho,0,j);
      R[imax-1][0] = grad(ho,imax-1,j);
      GT[0][0] = (2.0*ho[0][j]*grad(vo,0,j) + grad(ho,0,j)*vo[0][j]);
      GT[imax-1][0] = (2.0*ho[imax-1][j]*grad(vo,imax-1,j) + grad(ho,imax-1,j)*vo[imax-1][j]);
      cc[0][0] = grad(ho,0,j)/(2.0*ho[0][j]);
      cc[imax-1][0] = grad(ho,imax-1,j)/(2.0*ho[imax-1][j]);
      kpz[0][0] = kpz[1][0];
      kpz[imax-1][0] = kpz[imax-2][0];
      v1[0][0] = grad(vo,0,j);
      v1[imax-1][0] = grad(vo,imax-1,j);
      v2[0][0] = v2[1][0];
      v2[imax-1][0] = v2[imax-2][0];


      // length of leaf  
      length = 0.0;
      for(i=0; i<imax; i++) {
          length = length + sqrt(ho[i][j])*du;
          s[i][j] = length;   //length of plant up to n*di*du (di defined as a global constant)
        }

      if (t%timestep == 0) {
      	fprintf(plantLength,"%6d %6f \n",t,length); 

        sprintf(filename,filetag, "modLength",testnum,"_",t,".txt");
        sectionalGrowthData = fopen(filename,"w");

        for(i=0; i<imax; i++) {
            if (i%sectionStep == 0) fprintf(plantGrowthData,"%6d %6d %6d %e %e %e %e %e %e %e %e %e %e \n",t,i,j,ho[i][j],vo[i][j],R[i][j],GT[i][j],g[i][j],cc[i][j],kpz[i][j],v1[i][j],v2[i][j],s[i][j]);
            if (i%sectionStep == 0) fprintf(sectionalGrowthData, "%6d %6d %6f \n",t,i,s[i][j]);
          }
        fprintf(plantGrowthData, "%s\n", "");
      }
      
      
      // switch (t) {
      // case 10: case 1000: case 10000: case 20000: case 30000: case 40000: case 50000: case 60000: case 80000: case 100000: case 350000: case 500000: case 850000: case 1000000: case 1500000:
      //   sprintf(filename,filetag, "pm",testnum,"_",t,".txt");
      //   currentfile = fopen(filename,"w");
      //   sprintf(filename,filetag, "modLength",testnum,"_",t,".txt");
      //   sectionalGrowthData = fopen(filename,"w");
      //   for(i=0; i<imax; i++) {
      //     fprintf(currentfile,"%6d %6d %6d %e %e %e %e %e %e %e %e %e %e \n",t,i,j,ho[i][j],vo[i][j],R[i][j],GT[i][j],g[i][j],cc[i][j],kpz[i][j],v1[i][j],v2[i][j],s[i][j]);
      //     if (i%sectionStep == 0) fprintf(sectionalGrowthData, "%6d %6d %6f \n",t,i,s[i][j]);
      //   }
      //   break;      
      // default:
      // 	break;
      // }
      
      // update the metrics for next iteration
      for(i=0; i<imax; i++)
    	{
    	  hn[i][j] = ho[i][j];
    	  vn[i][j] = vo[i][j];
    	  ho[i][j] = hp[i][j];
    	  vo[i][j] = vp[i][j];
    	}
      
    }
  
  return 2;
  
}

/**********************************************************************************************/
/**********************************************************************************************/


int main()
{

  int i;
  int imax = 0;
  double di = (double) i;

  // plantModel(model #, c: curvature forces, c1: v*grad(v), c2: v damping, c3: v diffusion, k: ricci flow, k1: v*v, k2: growth tensor)

  for (i=0;i<imax+1;i++)
    {
      di = (double) i;

      // test cases; select one by uncommenting one line.

      // simulations with increasing sigmoidal velocity field (tanh)
      //if (plantModel(i,1.0,1.0,0.0,2.0,1.0,0.0,0.005) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // single test
      //if (plantModel(i,1.0,1.0,0.0,0.5,1.0,0.0,0.5) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // single test
      //if (plantModel(i,1.0,1.0,0.0,0.1,10.0,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); //  long time test (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.0,10.0,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); //  long time test null vel dif (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.0,0.0,0.0,0.005) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDiff, RF null
      // velocity still 'diffuses'... due to growth in metric? Just gets stretched out.

      //if (plantModel(i,1.0,1.0,0.0,1.0*di,1.0,0.0,0.005) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest1 (10^-4 dt)
      // steep vel: 0,1 are good values; sim fails mid way for 2,3.
      // shallow vel: 0 is good; 1,2,3 oscillate
      //if (plantModel(i,1.0,1.0,0.0,0.5*di,1.0,0.0,0.05) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest2 (10^-5 dt)
      // steep vel: highest coupling has ocsillations
      //if (plantModel(i,1.0,1.0,0.0,0.5*di,1.0,0.0,0.5) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest3 (10^-6 dt)
      // steep vel works very well
      ///if (plantModel(i,1.0,1.0,0.0,2.0*di,1.0,0.0,2.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest4 (10^-6 dt) 
      // steep vel: oscillatory for high coupling
      //if (plantModel(i,1.0,1.0,0.0,0.1*di,1.0,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest5 (10^-6 dt)
      // steep vel, silk vel: 0.2,0.3 too high
      //if (plantModel(i,1.0,1.0,0.0,0.2*di,10.0,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest6 (10^-6 dt)
      // silk vel: 0.4,0.6 too high
      //if (plantModel(i,1.0,1.0,0.0,0.1*di,40.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest7 (10^-6 dt)
      // silk vel: 0.2,0.3 too high
      // steep vel: 0.3 too high
      //if (plantModel(i,1.0,1.0,0.0,0.05*di,20.0,0.0,40.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest8 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,30.0*di,20.0,0.0,0.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest9 (10^-6 dt) 
      // works for all values up to VD=90.
      //if (plantModel(i,1.0,1.0,0.0,0.1*di,2.0,0.0,10.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest10 (10^-6 dt), linear vel
      //if (plantModel(i,1.0,1.0,0.0,0.1*di,2.0,0.0,1.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest11 (10^-6 dt), linear vel
      // small 'blips' in derivatives of metric
      //if (plantModel(i,1.0,1.0,0.0,0.2*di,40.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest12 (10^-6 dt) with crank-nicholson
      //if (plantModel(i,1.0,1.0,0.0,0.5*di,40.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest13 (10^-6 dt) with dufort-frankel
      //if (plantModel(i,1.0,1.0,0.0,2.0*di,40.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest14 (10^-6 dt) with Dufort-Frankel
      //if (plantModel(i,1.0,1.0,0.0,10.0*di,40.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest15 (10^-5 dt) with Dufort-Frankel
      //if (plantModel(i,1.0,1.0,0.0,10.0*di,10.0,0.0,10.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // velDifTest16 (10^-5 dt) with Dufort-Frankel


      //if (plantModel(i,1.0,1.0,0.0,0.1*di,1.0*di,0.0,0.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // VDRFtest1 (10^-6 dt), linear vel


      //if (plantModel(i,1.0,1.0,0.0,1.0,3.0*di,0.0,0.005) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest1 (10^-4 dt)
      // steep and shallow vel: zero RF fails part way through simulation; higher RF makes scale factor more uniform; length almost identical for all RF values.
      //if (plantModel(i,1.0,1.0,0.0,0.5,2.0*di,0.0,0.05) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest2 
      //if (plantModel(i,1.0,1.0,0.0,0.5,2.0*di,0.0,0.5) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest3; GT coupling too high with silk vel
      //if (plantModel(i,1.0,1.0,0.0,0.1,2.0*di,0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest4 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.1,5.0*(di+1.0),0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest5 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.1,10.0*(di+1.0),0.0,6.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest6 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,30.0,10.0*di,0.0,0.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // RFTest7 (10^-6 dt) Gaussian metric


      //if (plantModel(i,1.0,1.0,0.0,1.0,1.0,0.0,0.0015*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest1
      //if (plantModel(i,1.0,1.0,0.0,0.0,1.0,0.0,0.0015*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest2
      // steep vel: works great!
      // shallow vel: higest GT coupling gets oscillations 
      //if (plantModel(i,1.0,1.0,0.0,0.5,1.0,0.0,0.5*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest3 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.0,1.0,0.0,2.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest4 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.1,10.0,0.0,2.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest5 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.0,10.0,0.0,10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest6 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.1,20.0,0.0,20.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest7 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.1,2.0,0.0,10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest8 (10^-6 dt), linear vel
      //if (plantModel(i,1.0,1.0,0.0,0.1,20.0,0.0,10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest9 (10^-6 dt), linear vel
      //if (plantModel(i,1.0,1.0,0.0,10.0,20.0,0.0,10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest10 (10^-5 dt) with Dufort-Frankel
      if (plantModel(i,1.0,1.0,0.0,10.0,20.0,0.0,30.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTTest10 (10^-5 dt) with Dufort-Frankel


      //if (plantModel(i,1.0,1.0,0.0,0.0,20.0*(di+1.0),0.0,20.0*(di+1.0)) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTRFtest1 (10^-6 dt)
      //if (plantModel(i,1.0,1.0,0.0,0.2,20.0*(di+1.0),0.0,20.0*(di+1.0)) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTRFtest2 (10^-6 dt)
      // silk vel fails at long times
      //if (plantModel(i,1.0,1.0,0.0,0.0,5.0*(di+1.0),0.0,5.0*(di+1.0)) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // GTRFtest3 (10^-6 dt)


      //if (plantModel(i,0.0,0.0,0.0,0.0,10.0*di,0.0,10.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // SFadvection1 (10^-5 dt) with Dufort-Frankel, Gaussian metric
      //if (plantModel(i,0.0,0.0,0.0,0.0,10.0,0.0,10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // SFadvection2 (10^-5 dt) with Dufort-Frankel, Gaussian metric
      //if (plantModel(i,0.0,0.0,0.0,0.0,10.0,0.0,-10.0*di) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // SFadvection3 (10^-5 dt) with Dufort-Frankel, Gaussian metric

      //if (plantModel(i,1.0,1.0,0.0,10.0*di,0.0,0.0,0.0) == 2) printf("%s %d %s","model",i,"ran successfully \n"); // burgerTest (10^-5 dt) with Dufort-Frankel


    }

}

/*************************************************************************************/

