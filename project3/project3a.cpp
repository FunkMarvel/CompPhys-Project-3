#include <iostream>
#include <cmath>
#include "lib.h"

using namespace std;

int alpha=2;
double const  pi = 3.14159265359;

double int_function(double x1, double x2, double y1, double y2, double z1, double z2)
{
  double r1, r2, r12;
  r1 = sqrt(x1*x1+y1*y1+z1*z1);
  r2 = sqrt(x2*x2+y2*y2+z2*z2);
  r12 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  if (r12==0) {
    return 0;
  } else
  {
  return exp(-2*alpha*(r1+r2))*1/r12;
  }
}
//     Here we define various functions called by the main program
//     this function defines the function to integrate
int main()
{
     int n;
     double a, b;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *x = new double [n];
     double *w = new double [n];
//   set up the mesh points and weights
     gauleg(a, b, x, w, n);
//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum
     double int_gauss = 0.;
     for ( int i1 = 0;  i1 < n-1; i1++){
       for ( int i2 = 0;  i2 < n-1; i2++){
         for ( int i3 = 0;  i3 < n-1; i3++){
           for ( int i4 = 0;  i4 < n-1; i4++){
             for ( int i5 = 0;  i5 < n-1; i5++){
               for ( int i6 = 0;  i6 < n-1; i6++){
                  int_gauss+=w[i1]*w[i2]*w[i3]*w[i4]*w[i5]*w[i6]*int_function(x[i1],x[i2],x[i3],x[i4],x[i5],x[i6]);
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
      delete [] x;
      delete [] w;
      return 0;
}  // end of main program
//  this function defines the function to integrate // end of function to evaluate
