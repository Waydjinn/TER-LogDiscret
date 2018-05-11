#ifndef POHLIGHELLMAN_H_INCLUDED
#define POHLIGHELLMAN_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
/*#include "../PollardsRho/pollardrho.c"
#include "../BabyStepGiantStep/BabyStepGiantStep.c"*/

typedef struct PrimeFact
{
    mpz_t prime;
    mpz_t pow;
    mpz_t solution;
}PrimeFact;

void initTab(PrimeFact* P);
void affiche(PrimeFact* P);
void primeFactDecomp(PrimeFact* P, mpz_t n);
void chineseRemainder(mpz_t res, PrimeFact* P, mpz_t n);

#endif // POHLIGHELLMAN_H_INCLUDED
