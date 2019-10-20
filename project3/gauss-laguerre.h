/*
 * The definition module
 *                      gauss-laguerre.h
 * for the library function common for all C programs.
 */

 // Standard ANSI-C++ include files


#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
#define   UL         unsigned long
#define   EPS 3.0e-14
#define   MAXIT 10


// Function declarations

double gammln(double);
void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);
