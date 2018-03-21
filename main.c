#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

    int n;
    n = mpz_invert (cleMasquage, cleMasquage, p);
    if (n==0) exit(0);

    mpz_t tmp;
    mpz_init (tmp);
    mpz_mul (tmp, encMessage, cleMasquage);

    mpz_mod (encMessage, tmp, p);

}

void keyGen(mpz_t privateKey, mpz_t p) {

    unsigned long int seed = time(NULL);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm (privateKey, state, p);
    gmp_printf ("%s %Zd\n", "pvKey :", privateKey);
    gmp_randclear(state);

}

int main(int argc, char ** args){

    // Génération de la clé publique de Bob
    PubKey bob;
    mpz_init(bob.p);
    mpz_set_ui(bob.p, 8999);
    mpz_init(bob.alpha);
    mpz_set_ui(bob.alpha, 6426);
    mpz_init(bob.beta);

    mpz_t x;
    mpz_init(x);
    keyGen(x, bob.p);
    mpz_powm(bob.beta, bob.alpha, x, bob.p);

    // Chiffrement d'un message
    mpz_t message;
    mpz_init(message);
    mpz_set_ui(message, 502);
    gmp_printf ("%s %Zd\n", "clair = ", message);

    mpz_t y;    // clé privée d'Alice
    mpz_init(y);
    keyGen(y, bob.p);
    encrypt(message, y, bob);
    gmp_printf ("%s %Zd\n", "chiffre = ", message);

    // Dechiffrement du message chiffré
    mpz_t cleE;
    mpz_init (cleE);
    mpz_powm(cleE, bob.alpha, y, bob.p);

    decrypt(message, cleE, x, bob.p);
    gmp_printf ("%s %Zd\n", "clair = ", message);

    mpz_clear(bob.p);
    mpz_clear(bob.alpha);
    mpz_clear(bob.beta);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(message);
    mpz_clear(cleE);

    return 0;

}
