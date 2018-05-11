#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



int main(int argc, char **argv)
{
	/*mpz_t ordre;
	mpz_init(ordre);

	mpz_set_ui(ordre,8);
	

	mpz_set(ordre, mpz_sizeinbase(ordre,2));

	gmp_printf("%Zd ", ordre);*/

	double n = 1217.0;
	double logn = log(n);
	double demi = 0.5;

	double res0 =  log( logn );
	double res1 =  res0 * logn;
	double res2 = pow(res1, demi);
	double res3 = (sqrt(2)/2) + 1; //OK
	double res4 = res3 * res2;
	double res5 = exp(res4);

	double res = exp( ( (sqrt(2)/2) + 1 )  *  pow( (  log( log(n) ) * log(n) ), 0.5 ) );

	printf("\nres = %lf ; log(%lf) = %lf \n", res, n, logn);
	printf("\nres0 = %lf res1 = %lf ; res2 = %lf ; res3 = %lf ; res4 = %lf ; res5 = %lf ; \n", res0 , res1, res2, res3, res4, res5);


	return 0;
}