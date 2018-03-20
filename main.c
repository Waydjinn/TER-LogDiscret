#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct {
    mpz_t p;        // nombre premier
    mpz_t alpha;    // élément primitif de Zp*
    mpz_t beta;     // égal à alpha^d mod p, avec d clé privée
}PubKey;

mpz_t encrypt (mpz_t message, mpz_t key, PubKey bob) {

    mpz_t cleEphemere;
    mpz_init (cleEphemere);
    mpz_powm (cleEphemere, bob->alpha, key, bob->p);

    mpz_t cleMasquage;
    mpz_init (cleMasquage);
    mpz_powm (cleMasquage, bob->beta, key, bob->p);

    mpz_t tmp;
    mpz_init (tmp);
    mpz_mul (tmp, message, cleMasquage);
    mpz_t encMessage;
    mpz_init (encMessage);
    mpz_mod (encMessage, tmp, bob->p);

    return encMessage;
}

mpz_t decrypt (mpz_t encMessage, mpz_t cleEphemere, mpz_t key, mpz_t p) {

    mpz_t cleMasquage;
    mpz_init (cleMasquage);
    mpz_powm (cleMasquage, cleEphemere, key, p);

    //Calcul de l'inverse de cleMasquage : Euclide étendu + exponentiation rapide (Square And Multiply)
    int n;
    n = mpz_invert (cleMasquage, cleMasquage, p);
    // fin

    mpz_t tmp;
    mpz_init (tmp);
    mpz_mul (tmp, encMessage, cleMasquage);

    mpz_t message;
    mpz_init (message);
    mpz_mod (message, tmp, p);

    return message;

}

int main(){

    return 0;
}
