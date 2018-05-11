#include "pohligHellman.h"

void initTab(PrimeFact* P)
{
    int i;
    for (i=0; i<1000; i++)
    {
        mpz_init(P[i].prime);
        mpz_init(P[i].pow);
        mpz_init(P[i].solution);
        mpz_set_ui(P[i].prime, 1);
        mpz_set_ui(P[i].pow, 0);
    }
}

void affiche(PrimeFact* P)
{
    int i=0;
    while(mpz_cmp_ui(P[i].prime, 1) != 0)
    {
        gmp_printf("Facteur : %Zd\nPuissance : %Zd\n\n", P[i].prime, P[i].pow);
        i++;
    }
}

void primeFactDecomp(PrimeFact* P, mpz_t n)      // décomposition de n en facteurs premiers
{
    mpz_t copyn;
    mpz_init_set(copyn, n);

    mpz_t rootn;
    mpz_init(rootn);
    mpz_sqrt(rootn, copyn);
    mpz_add_ui(rootn, rootn, 1);

    mpz_t i;
    mpz_init_set_ui(i, 2);

    mpz_t puis;
    mpz_init_set_ui(puis, 0);

    int indice = 0;

    while (mpz_cmp(i, rootn) <= 0)          // Tant que les facteurs testés sont inférieurs ou égaux à la racine de n
    {
        if(mpz_probab_prime_p(copyn, 40) > 0)
        {                   // Si n est premier ou supposé premier après 40 tests on le renvoie
            mpz_init_set(P[indice].prime, copyn);
            mpz_init_set_ui(P[indice].pow, 1);
            indice++;

            break;
        }

        mpz_t test;
        mpz_init(test);
        do
        {
            mpz_mod(test, copyn, i);
            if(mpz_cmp_ui(test, 0) == 0)        // Si i divise n
            {
                mpz_div(copyn, copyn, i);
                mpz_add_ui(puis, puis, 1);
            }
        }while (mpz_cmp_ui(test, 0) == 0);

        if (mpz_cmp_ui(puis, 0) != 0)       // On n'affiche pas les puissance à 0.
        {
            mpz_set(P[indice].prime, i);
            mpz_set(P[indice].pow, puis);
            indice++;
        }
        mpz_nextprime(i, i);
        mpz_set_ui(puis, 0);

        mpz_clear(test);
    }
}

void chineseRemainder(mpz_t res, PrimeFact* P, mpz_t n)     //Appliquer le théorème des restes chinois pour trouver la solution finale
{
    printf("\n\n****************** CHINESE REMAINDER ***********************\n\n");
    int i = 0;


    while(mpz_cmp_ui(P[i].prime, 1) != 0)
    {
        mpz_t M;
        mpz_init(M);

        mpz_t inv;
        mpz_init(inv);

        mpz_t tmp;      // tmp = pi^ui mod n
        mpz_init(tmp);
        mpz_powm(tmp, P[i].prime, P[i].pow, n);
        mpz_div(M, n, tmp);     // M = n/(pi^ui)
        gmp_printf("M = %Zd\n", M);

        mpz_invert(inv, M, tmp);
        mpz_mod(inv, inv, tmp);
        gmp_printf("l'inverse de M mod %Zd est : %Zd\n", tmp, inv);

        mpz_t xtmp;        // xtmp = xi * M^(-1) * M
        mpz_init(xtmp);
        mpz_mul(xtmp, P[i].solution, inv);
        mpz_mul(xtmp, xtmp, M);
        gmp_printf("xtmp = %Zd\n", xtmp);

        mpz_add(res, res, xtmp);
        gmp_printf("res = %Zd\n\n", res);
        i++;

        mpz_clear(M);
        mpz_clear(inv);
        mpz_clear(tmp);
        mpz_clear(xtmp);
    }

    mpz_mod(res, res, n);
    gmp_printf("x = %Zd\n\n", res);
}
