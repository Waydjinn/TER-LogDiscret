#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void naif_methode(mpz_t res, mpz_t ordre, mpz_t generateur, mpz_t elt){

	mpz_t k;
	mpz_init(k);
	mpz_set_ui(k,0);

	do{
		//On teste toutes les puissances possibles
		mpz_powm(res,generateur,k,ordre);

	//Quand on trouve le résultat
	if(mpz_cmp(res,elt)==0){
		gmp_printf("Le résultat est : k = %Zd \n", k);
		mpz_set(res,k);
		mpz_clear(k);
		return;
	}

	mpz_add_ui(k,k,1);
	}while(1);
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

	//////////////////////////////////////////////////////////////////////////
	//On peut aussi utiliser mpz_set_str pour entrer des arguments plus grands
	//////////////////////////////////////////////////////////////////////////
	/*
	mpz_set_str(generateur,"2", 10);//g=23
	mpz_set_str(ordre,"1000000000121",10);//n=53
	mpz_set_str(elt,"99999999", 10);*/
	
	naif_methode(resultat,ordre,generateur,elt);

	return 1;
}
