#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>




void baby_step_giant_step(mpz_t resultat,mpz_t ordre,mpz_t generateur, mpz_t eltB){
	

	
	
	int x;
	
	mpz_t resultatF;
	mpz_init(resultatF);
	
	mpz_t m,entier;
	mpz_init(m);mpz_init(entier);
	

	
	
	mpz_sqrt(m,ordre); //m = (racine carré(ordre))
	mpz_set_ui(entier,mpz_get_ui(m)+1); //entier = partie entière (racine carré(ordre)) + 1
	
	
	
	mpz_t* t = malloc(mpz_get_ui(entier)*sizeof(mpz_t)); //table de hachage : allouer de la mémoire de m entiers
	
	//peut être une boucle
	for(x=0;x<mpz_get_ui(entier);x++){
		mpz_init(t[x]);
	}
	
	
	//==============================================================================================
	
	
	int j;

	for(j=0;j<mpz_get_ui(entier);j++){ //Baby Step
		
		
		mpz_t res,tmp;
		mpz_init(res);mpz_init(tmp);
		mpz_set_ui(tmp,j);
		
	
		mpz_powm(res,generateur,tmp,ordre); //res = generateur^j mod ordre 			
		
		
	
		
		mpz_set(t[j],res); //stocker dans une table de hachage
	
	}
	
	
	mpz_t res1,g_inv,mpz_entier,gamma; //Calculer alpha^(-m) mod n
	mpz_init(res1);mpz_init(g_inv);mpz_init(mpz_entier),mpz_init(gamma);
	

	
	int inv = mpz_invert(g_inv,generateur,ordre);
	mpz_powm(res1,g_inv,entier,ordre);
	
	
	mpz_set(gamma,eltB);

	int i;
	
	
	for(i=0;i<mpz_get_ui(entier)-1;i++){ //Giant Step
		
		if(mpz_cmp(t[i],gamma)==0){
			free(t);
			mpz_mul_ui(resultatF,entier,i); // resF = entier * i
			mpz_add_ui(resultatF,resultatF,i); // resF = (entier * i) +1
			mpz_set(resultat,resultatF); //affectation 
			
			
		}else{
		
			mpz_mul(gamma,gamma,res1);
			mpz_mod(gamma,gamma,ordre);
			
		}
	}
	
		
	
	
	}
	
int main (int argc, char *argv[]){
	
	mpz_t a,ordre,elt;
	mpz_init(a);mpz_init(ordre);mpz_init(elt);
	mpz_t resultat;
	mpz_init(resultat);
	
	mpz_set_ui(a,11);
	mpz_set_ui(ordre,6);
	mpz_set_ui(elt,10);
	
	baby_step_giant_step(resultat,a,ordre,elt);
	gmp_printf("resultat=%Zd \n",resultat); //a=5
	
	
	
	return 1;
}
