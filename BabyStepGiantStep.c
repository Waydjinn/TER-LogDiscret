#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>




int baby_step_giant_step(int ordre,int generateur, int eltB){
	
	
	mpz_t m,n,g;
	mpz_init(m);mpz_init(n);mpz_init(g); //initialisation : utile ? 
	
	mpz_set_ui(n,ordre); //ordre mpz_t
	mpz_set_ui(g,generateur); //générateur mpz_t
	mpz_pow_ui(m,n,0.5); //initialisation de m = (racine carré(m)) + 1
	int entier = mpz_get_ui(m)+1; // pour avoir la partie entière de m
	
	int* t = malloc(mpz_get_ui(m)*sizeof(int)); //table de hachage : allouer de la mémoire de m entiers
	
	//==============================================================================================
	
	int j;
	for(j=0;j<mpz_get_ui(m);j++){ //Baby Step
		
		mpz_t res,tmp;
		mpz_set_ui(tmp,j);
		mpz_powm(res,g,tmp,n); //alpha^j mod n : n = ordre du groupe
		
		int resultat=mpz_get_ui(res);
		t[j]=resultat; //stocker dans une table de hachage
	
	}
	mpz_t res1; //Calculer alpha^(-m) mod n
	mpz_powm(res1,g,m,n);
	
	int gamma=eltB;
	int i;
	
	for(i=0;i<(mpz_get_ui(m)-1);i++){ //Giant Step
		if(t[i]==gamma){
			return i*entier + i;
		}else{
			gamma=(gamma*mpz_get_ui(res1))%ordre;
		}
	}
	
		
	free(t);
	return 0; 
	}
