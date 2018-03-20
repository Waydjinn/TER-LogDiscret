#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct {
    mpz_t p;        // nombre premier
    mpz_t alpha;    // élément primitif de Zp*
    mpz_t beta;     // égal à alpha^d mod p, avec d clé privée
}PubKey;

void encrypt (mpz_t message, mpz_t key, PubKey bob) {

    mpz_t cleEphemere;
    mpz_init (cleEphemere);
    mpz_powm (cleEphemere, bob.alpha, key, bob.p);

    mpz_t cleMasquage;
    mpz_init (cleMasquage);
    mpz_powm (cleMasquage, bob.beta, key, bob.p);

    mpz_t tmp;
    mpz_init (tmp);
    mpz_mul (tmp, message, cleMasquage);

    mpz_mod (message, tmp, bob.p);

}

void decrypt (mpz_t encMessage, mpz_t cleEphemere, mpz_t key, mpz_t p) {

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

    mpz_mod (encMessage, tmp, p);

}

int main(){

    PubKey bob;
    mpz_init(bob.p);
    mpz_set_ui(bob.p, 8999);
    mpz_init(bob.alpha);
    mpz_set_ui(bob.alpha, 6426);
    mpz_init(bob.beta);
    mpz_t x;
    mpz_init(x);
    mpz_set_ui(x, 3659);
    mpz_powm(bob.beta, bob.alpha, x, bob.p);

    mpz_t message;
    mpz_init(message);
    mpz_set_ui(message, 123456789);
    gmp_printf ("%s %Zd\n", "clair = ", message);

    mpz_t y;
    mpz_init(y);
    mpz_set_ui(y, 4563);
    encrypt(message, y, bob);
    gmp_printf ("%s %Zd\n", "chiffre = ", message);

    mpz_t cleE;
    mpz_init (cleE);
    mpz_powm(cleE, bob.alpha, y, bob.p);

    decrypt(message, cleE, x, bob.p);
    gmp_printf ("%s %Zd\n", "chiffre = ", message);

    return 0;
}
