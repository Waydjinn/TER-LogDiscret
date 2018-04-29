#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void index_calculus(/*mpz_t resultat,*/ mpz_t ordre, mpz_t generateur, mpz_t elt/*, int* factorBase*/){
	
	mpz_t res_puissance;
	mpz_init(res_puissance);
	
	int i;
	for(i=1;i<100;i++){
		mpz_powm_ui(res_puissance,generateur,i,ordre);
		gmp_printf("%Zd ^(%d) = %Zd \n",generateur,i,res_puissance);
	}
}

int main(int argc, char **argv)
{
	mpz_t generateur,ordre,elt;
	mpz_init(generateur);mpz_init(ordre);mpz_init(elt);
	
	mpz_set_ui(generateur,3);
	mpz_set_ui(ordre,1217);
	mpz_set_ui(elt,37);
	
	index_calculus(ordre,generateur,elt);
	
	return 0;
}

