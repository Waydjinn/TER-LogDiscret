#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>




int baby_step_giant_step(int ordre,int generateur, int eltB){
	

	
	mpz_t m,n,g;
	mpz_init(m);mpz_init(n);mpz_init(g); //initialisation : utile ? 
	
	
	

	mpz_set_ui(n,ordre); //ordre mpz_t	
	mpz_set_ui(g,generateur); //générateur mpz_t

	
	
	mpz_sqrt(m,n); //initialisation de m = (racine carré(n)) + 1
	int entier = mpz_get_ui(m)+1; // pour avoir la partie entière de m
	

	
	int* t = malloc(entier*sizeof(int)); //table de hachage : allouer de la mémoire de m entiers
	
	
	
	//==============================================================================================
	
	
	int j;
	for(j=0;j<entier;j++){ //Baby Step
		
		
		mpz_t res,tmp;
		mpz_init(res);mpz_init(tmp);
		mpz_set_ui(tmp,j);
		
	
		mpz_powm(res,g,tmp,n); //res = alpha^j mod n 			: n = ordre du groupe et alpha = g
		
		
		int resultat=(int)mpz_get_ui(res);
		
		t[j]=resultat; //stocker dans une table de hachage
	
	}
	
	
	mpz_t res1,g_inv,mpz_entier; //Calculer alpha^(-m) mod n
	mpz_init(res1);mpz_init(g_inv);mpz_init(mpz_entier);
	
	mpz_set_ui(mpz_entier,entier);
	
	int inv = mpz_invert(g_inv,g,n);
	
	mpz_powm(res1,g_inv,mpz_entier,n);
	
	int gamma=eltB;
	int i;
	
	
	for(i=0;i<entier-1;i++){ //Giant Step
		if(t[i]==gamma){
			free(t);
			return i*entier + i;
			
		}else{
			gamma=(gamma*mpz_get_ui(res1))%ordre;
		}
	}
	
		
	
	return 0; 
	}
	
int main (int argc, char *argv[]){
	int a = baby_step_giant_step(11,6,10);
	printf("a=%d \n",a); //a=5
	return 1;
}
