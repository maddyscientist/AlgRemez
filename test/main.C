/*
  Mike Clark - 25th May 2005

  Quick test code for calculating optimal rational approximations for
  the functions x^(y/z) with appropriate bounds.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"alg_remez.h"

int main (int argc, char* argv[]) {

  if (argc != 8) {
    printf("./test <numerator> <denominator> <deg numerator> <deg denominator> <lambda_low> <lambda_high> <precision>\n");
    exit(0);
  }  
  
  int i=0;
  int n; // The degree of the numerator polynomial
  int d; // The degree of the denominator polynomial
  int y; // The numerator of the exponent
  int z; // The denominator of the exponent
  int precision; // The precision that gmp uses
  double lambda_low, lambda_high; // The bounds of the approximation

  // Set the exponent
  sscanf(argv[++i],"%d",&y);
  sscanf(argv[++i],"%d",&z);  
  if(y <= 0 || z <= 0) {
    printf("Both the numerator y=%d and the denominator z=%d must be positive.\n", y, z);
    exit(0);
  }
  
  // Set the required degree of approximation
  sscanf(argv[++i],"%d",&n);
  sscanf(argv[++i],"%d",&d);
  if(n <= 0 || d <= 0) {
    printf("Both the numerator degree n=%d and the denominator degree d=%d must be positive.\n", n, d);
    exit(0);
  }

  // Set the approximation bounds
  sscanf(argv[++i],"%le",&lambda_low);
  sscanf(argv[++i],"%le",&lambda_high);
  if(lambda_low <= 0 || lambda_high <= 0) {
    printf("Both the lower bound lambda_low=%e and the upper bound lambda_high=%e must be positive.\n", lambda_low, lambda_high);
    exit(0);
  }

  // Set the precision of the arithmetic
  sscanf(argv[++i],"%d",&precision);
  if(precision <= 0) {
    printf("The bit length precision=%d must be positive.\n", precision);
    exit(0);
  }
  
  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res = new double[n];
  double *pole = new double[d];

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));
  char FORCE_FILE[100], ENERGY_FILE[100];
  sprintf(FORCE_FILE, "force_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);
  sprintf(ENERGY_FILE, "energy_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);

  // Instantiate the Remez class
  AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez.generateApprox(n,d,y,z);

  FILE *output = fopen("approx.dat", "w");

  fprintf(output, "Approximation to f(x) = x^(%d/%d)\n\n", y, z);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez.getPFE(res,pole,&norm);
  
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	    i+1, res[i], i+1, pole[i]);
  }

  // Find pfe of inverse function
  remez.getIPFE(res,pole,&norm);
  fprintf(output, "\nApproximation to f(x) = x^(-%d/%d)\n\n", y, z);
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	    i+1, res[i], i+1, pole[i]);
  }

  fclose(output);

  FILE *error_file = fopen("error.dat", "w");
  for (double x=lambda_low; x<lambda_high; x*=1.01) {
    double f = remez.evaluateFunc(x);
    double r = remez.evaluateApprox(x);
    fprintf(error_file,"%e %e\n", x,  (r - f)/f);
  }
  fclose(error_file);

  delete[] res;
  delete[] pole;

  exit(0);

}     
