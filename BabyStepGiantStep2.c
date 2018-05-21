#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int baby_step_giant_step(mpz_t resultat,mpz_t ordre,mpz_t generateur, mpz_t eltB){	
	
	int x;
	
	mpz_t resultatF;
	mpz_init(resultatF);
	
	mpz_t m,entier;
	mpz_init(m);mpz_init(entier);
	
	mpz_sqrt(m,ordre); //m = (racine carré(ordre))
	mpz_set_ui(entier,mpz_get_ui(m)+1); //entier = partie entière (racine carré(ordre)) + 1
	
	mpz_clear(m);
	
	//Table : allouer de la mémoire pour m entiers
	mpz_t* t = malloc(mpz_get_ui(entier)*sizeof(mpz_t)); 
	
	//On l'initialise
	for(x=0;x<mpz_get_ui(entier);x++){
		mpz_init(t[x]);
	}
	
	//======================================Baby Step========================================================
	
	
	int j;

	for(j=0;j<mpz_get_ui(entier);j++){ //Baby Step
		
		
		mpz_t res,tmp;
		mpz_init(res);mpz_init(tmp);
		mpz_set_ui(tmp,j);
		
	
		mpz_powm(res,generateur,tmp,ordre); //res = generateur^j mod ordre 			
		
		mpz_set(t[j],res); //On stocke dans la table 
		
		mpz_clear(tmp);
		mpz_clear(res);
	
	}
	
	//====================================Calcul de g⁽⁻m⁾ : g_inv⁽m⁾ ===================================================
	
	mpz_t res1,g_inv,gamma; //Calculer alpha^(-m) mod n
	mpz_init(res1);mpz_init(g_inv);mpz_init(gamma);
	
	int inv = mpz_invert(g_inv,generateur,ordre);//calcul de l'inverse de g
	
	mpz_powm(res1,g_inv,entier,ordre); //g_inv à la puissance m
	
	mpz_clear(g_inv);
	
	//===============================Giant Step==========================================================
	
	mpz_set(gamma,eltB);

	int i;
	
	for(i=0;i<mpz_get_ui(entier);i++){ //Giant Step
		
		for(j=0;j<mpz_get_ui(entier);j++){
		
			if(mpz_cmp(t[j],gamma)==0){
				printf("Solution trouvée ! \n");
			
				free(t);
			
				mpz_mul_ui(resultatF,entier,i); // resF = entier * i
				mpz_add_ui(resultatF,resultatF,j); // resF = (entier * i) + j
				mpz_set(resultat,resultatF); //affectation 
			
				mpz_clear(resultatF);
				mpz_clear(entier);
				mpz_clear(gamma);
				mpz_clear(res1);
				gmp_printf("Le résultat = %Zd \n",resultat);
				
				return 0;
			}	
		}
		
			mpz_mul(gamma,gamma,res1);
			mpz_mod(gamma,gamma,ordre);	
	}
}
	
int main (int argc, char *argv[]){
	
	mpz_t generateur, ordre, elt, reponse, resultat;
	mpz_init(generateur); mpz_init(ordre); mpz_init(elt); mpz_init(reponse), mpz_init(resultat);
	
	/////////////////////////////////////////
	//Plusieurs valeurs possibles à essayer :
	/////////////////////////////////////////

	/*mpz_set_ui(generateur,872052);
	mpz_set_ui(ordre,1000667);
	mpz_set_ui(elt,588654);*/

	/*
	mpz_set_ui(generateur,78277094);//g=23
	mpz_set_ui(ordre,100000127);//n=53
	mpz_set_ui(elt,42);//h=3 //ligne 4 */

	/*
	mpz_set_ui(generateur,78277094);//g=23
	mpz_set_ui(ordre,100000127);//n=53
	mpz_set_ui(elt,7777);//h=3 //ligne 4 */

	mpz_set_ui(generateur,5430827534);//g=23
	mpz_set_ui(ordre,10000000259);//n=53
	mpz_set_ui(elt,2090591416);

	/*
	mpz_set_ui(generateur,5430827534);//g=23
	mpz_set_ui(ordre,10000000259);//n=53
	mpz_set_ui(elt,4278787);
	 */

	/*
	mpz_set_ui(generateur,3);//g=23
	mpz_set_ui(ordre,10);//n=53
	mpz_set_ui(elt,4);
	 */
	
	/*
	mpz_set_ui(generateur,654438348730);//g=23
	mpz_set_ui(ordre,1000000000547);//n=53
	mpz_set_ui(elt,24333);
	 */

	int g = baby_step_giant_step(resultat,ordre,generateur,elt);
	
	return 1;
}
