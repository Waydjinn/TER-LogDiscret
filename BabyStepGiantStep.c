#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>




int baby_step_giant_step(mpz_t resultat,mpz_t ordre,mpz_t generateur, mpz_t eltB){	
	
	int x;
	
	printf("Début du programme. \n");
	
	mpz_t resultatF;
	mpz_init(resultatF);
	
	mpz_t m,entier;
	mpz_init(m);mpz_init(entier);
	

	
	
	mpz_sqrt(m,ordre); //m = (racine carré(ordre))
	mpz_set_ui(entier,mpz_get_ui(m)+1); //entier = partie entière (racine carré(ordre)) + 1
	
	mpz_clear(m);
	
	mpz_t* t = malloc(mpz_get_ui(entier)*sizeof(mpz_t)); //table de hachage : allouer de la mémoire de m entiers
	
	//peut être une boucle
	for(x=0;x<mpz_get_ui(entier);x++){
		mpz_init(t[x]);
	}
	
	
	//======================================Baby Step========================================================
	
	
	int j;
	printf("Baby Step \n");
	for(j=0;j<mpz_get_ui(entier);j++){ //Baby Step
		
		
		mpz_t res,tmp;
		mpz_init(res);mpz_init(tmp);
		mpz_set_ui(tmp,j);
		
	
		mpz_powm(res,generateur,tmp,ordre); //res = generateur^j mod ordre 			
		
		
	
		
		mpz_set(t[j],res); //stocker dans une table de hachage
		
		mpz_clear(tmp);
		mpz_clear(res);
	
	}
	
	printf("On affiche la table de hashage. \n");
	for(x=0;x<mpz_get_ui(entier);x++){
		gmp_printf("t[%d]=%Zd \n",x,t[x]);
	}
	
	//====================================Calcul de g⁽⁻m⁾ : g_inv⁽m⁾ ===================================================
	
	mpz_t res1,g_inv,gamma; //Calculer alpha^(-m) mod n
	mpz_init(res1);mpz_init(g_inv);mpz_init(gamma);
	

	
	int inv = mpz_invert(g_inv,generateur,ordre);//calcul de l'inverse de g
	gmp_printf("g_inverse : generateur = %Zd  , ordre = %Zd \n",generateur,ordre);
	
	mpz_powm(res1,g_inv,entier,ordre); //g_inv à la puissance m
	
	gmp_printf("%Zd à la puissance %Zd = %Zd (modulo %Zd) \n",g_inv,entier,res1,ordre);
	
	mpz_clear(g_inv);
	
	//===============================Giant Step==========================================================
	
	mpz_set(gamma,eltB);



	int i;
	
	printf("Giant Step \n");
	for(i=0;i<mpz_get_ui(entier);i++){ //Giant Step
		printf("Pour i=%d \n",i);
		//gmp_printf("t[%d]=%Zd compare gamma = %Zd \n",i,t[i],gamma);
		
		for(j=0;j<mpz_get_ui(entier);j++){
			
			gmp_printf("Pour j=%d et t[j]=%Zd \n",j,t[j]);
		
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
		
			gmp_printf("gamma = %Zd * %Zd \n",gamma,res1);
			mpz_mul(gamma,gamma,res1);
			
			mpz_mod(gamma,gamma,ordre);
			gmp_printf("Je modifie gamma !\n",gamma);
			
		
	}
	
		
	
	
	}
	
int main (int argc, char *argv[]){
	
	mpz_t generateur,ordre,elt;
	mpz_init(generateur);mpz_init(ordre);mpz_init(elt);
	mpz_t resultat;
	mpz_init(resultat);
	
	//
	mpz_set_ui(generateur,2);//g=23
	mpz_set_ui(ordre,1019);//n=53
	mpz_set_ui(elt,5);//h=3    6¹⁰ = 60466176 mod 11 = 1
	
	int g = baby_step_giant_step(resultat,ordre,generateur,elt);
	//gmp_printf("resultat=%Zd \n",resultat); //a=5
	
	
	
	return 1;
}
