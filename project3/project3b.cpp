#include <iostream>
#include <cmath>
#include "gauss-laguerre.h"

using namespace std;

int alpha=2;
double const  pi = 3.14159265359;

void LinearSpacedArray(double *x, double a, double b, int n){
  double h = (b-a)/n;
  x[0] = a;
  x[-1] = b;
  for (int i=1;i<n-1;i++){
    x[i] = x[i-1] + h;
  }
}

double int_function(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
  double numerator = sin(theta1)*sin(theta2);
  double denominator = 1024.*sqrt(r1*r1 + r2*r2 - 2.*r1*r2*(cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2)));
  if (denominator == 0.) {
    return 0;
  }
  else return numerator/denominator;
}
//     Here we define various functions called by the main program
//     this function defines the function to integrate
int main()
{
     int n;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *wr1 = new double [n+1];
     double *wr2 = new double [n+1];
     double *wtheta1 = new double [n+1];
     double *wtheta2 = new double [n+1];
     double *wphi1 = new double [n+1];
     double *wphi2 = new double [n+1];

     double *r1 = new double [n+1];
     double *r2 = new double [n+1];
     double *theta1 = new double [n+1];
     double *theta2 = new double [n+1];
     double *phi1 = new double [n+1];
     double *phi2 = new double [n+1];
//   set up the mesh points and weights
     LinearSpacedArray(r1, 0, 3, n+1);
     LinearSpacedArray(r2, 0, 3, n+1);
     LinearSpacedArray(theta1, 0, pi, n+1);
     LinearSpacedArray(theta2, 0, pi, n+1);
     LinearSpacedArray(phi1, 0, 2*pi, n+1);
     LinearSpacedArray(phi2, 0, 2*pi, n+1);

     gauss_laguerre(r1, wr1, n, alpha);
     gauss_laguerre(r2, wr2, n, alpha);
     gauss_laguerre(theta1, wtheta1, n, alpha);
     gauss_laguerre(theta2, wtheta2, n, alpha);
     gauss_laguerre(phi1, wphi1, n, alpha);
     gauss_laguerre(phi2, wphi2, n, alpha);
//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum
     double int_gauss = 0.;
     for ( int i1 = 1;  i1 < n; i1++){
       for ( int i2 = 1;  i2 < n; i2++){
         for ( int i3 = 1;  i3 < n; i3++){
           for ( int i4 = 1;  i4 < n; i4++){
             for ( int i5 = 1;  i5 < n; i5++){
               for ( int i6 = 1;  i6 < n; i6++){
                  int_gauss+=wr1[i1]*wr2[i2]*wtheta1[i3]*wtheta2[i4]*wphi1[i5]*wphi2[i6]*int_function(r1[i1],r2[i2],theta1[i3],theta2[i4],phi1[i5],phi2[i6]);
               }
             }
           }
         }
       }
     }
     double e1 = int_gauss - 5*pi*pi/(16*16);
     double error = sqrt(e1*e1);
//    final output
      cout << "Gaussian quad = " << int_gauss << " Error = " << error << endl;
      delete [] r1;
      delete [] wr1;
      delete [] r2;
      delete [] wr2;
      delete [] theta1;
      delete [] wtheta1;
      delete [] theta2;
      delete [] wtheta2;
      delete [] phi1;
      delete [] wphi1;
      delete [] phi2;
      delete [] wphi2;
      return 0;
}  // end of main program
//  this function defines the function to integrate // end of function to evaluate
