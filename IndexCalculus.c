#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//on cherche la plus grande puissance d'un nombre premier composant le nombre elt
void prime_pow(mpz_t resultat,int prime, int* factorBase,int tailleFB, mpz_t elt){
	mpz_t res;
	mpz_init(res);
	
	int i;
	int cmp=1;
	
	for(i=0;i<tailleFB;i++){ //on stock dans cmp la multiplication de tous les nombres de factorBase
		
		if(factorBase[i]!=prime){ //sauf pour le nbr dont on veut trouver la puissance max
			cmp=cmp*factorBase[i];
		}
	
	}
	
	mpz_powm_ui(res,cmp,1,ordre); //je met cmp à la puissance 1 et je lui applique le modulo ordre
	
	//donc dans le mpz_t res j'ai récupérer un nombre prime à la puissance x (le x que je cherche)
	
	
}

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

