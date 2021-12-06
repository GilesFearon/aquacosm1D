#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#ifndef M_PI
  #define M_PI (3.14159265358979323846)
#endif

void diffuse(double *Particles, double *SEM,
	     int Npts, int Nscalars,
	     double scalefct, double denom, double radius){
  int Pstride = Nscalars+2;
  int i, j, ns;
  double d, q, MassExch;
  double *aux;

  //Zeroes the elements of the 'Sum of Exchanged Masses' array.
  for (aux=SEM; aux<SEM+Npts*Nscalars; aux++){
    *aux = 0.;
  }
  //Exchanges the scalars among the i-th particle and all other
  //particles closer than 'radius'.
  for (i=0; i<Npts-1; i++){
    j = i+1;
    while (true){
      if (j==Npts) break;
      d = Particles[j*Pstride+1] - Particles[i*Pstride+1];
      if  (d>=radius) break;
      q = scalefct*exp(-d*d/denom);
      for (ns=0; ns<Nscalars; ns++){
  	MassExch = q*(Particles[j*Pstride+ns+2] - Particles[i*Pstride+ns+2]);
  	SEM[i*Nscalars+ns] += MassExch;
  	SEM[j*Nscalars+ns] -= MassExch;
      }
      j += 1;
    }
  }

  //Updates the concentration of scalars of each particle
  for (i=0; i<Npts; i++){
    for (ns=0; ns<Nscalars; ns++){
      Particles[i*Pstride+ns+2] += SEM[i*Nscalars+ns];
    }
  }
}

/*-------------------------------------------------------------*/

void diffuse_varD(double *Particles, double *SEM,
		  int Npts, int Nscalars,
		  double p, double radius, double dt,
		  double *EddyDiff
		  ){
  int Pstride = Nscalars+2;
  int i, j, ns;
  double d, q, MassExch, denom;
  double keddy;
  double *aux;

  //Zeroes the elements of the 'Sum of Exchanged Masses' array.
  for (aux=SEM; aux<SEM+Npts*Nscalars; aux++){
    *aux = 0.;
  }
  //Exchanges the scalars among the i-th particle and all other
  //particles closer than 'radius'.
  for (i=0; i<Npts-1; i++){
    j = i+1;
    while (true){
      if (j==Npts) break;
      d = Particles[j*Pstride+1] - Particles[i*Pstride+1];
      if (d>=radius) break;
      keddy = fmin(EddyDiff[j], EddyDiff[i]);
      denom = 4*keddy*dt;
      q = p*exp(-d*d/denom)/sqrt(M_PI*denom);
      for (ns=0; ns<Nscalars; ns++){
  	MassExch = q*(Particles[j*Pstride+ns+2] - Particles[i*Pstride+ns+2]);
  	SEM[i*Nscalars+ns] += MassExch;
  	SEM[j*Nscalars+ns] -= MassExch;
      }
      j += 1;
    }
  }

  //Updates the concentration of scalars of each particle
  for (i=0; i<Npts; i++){
    for (ns=0; ns<Nscalars; ns++){
      Particles[i*Pstride+ns+2] += SEM[i*Nscalars+ns];
    }
  }
}
