#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "../PollardsRho/pollardrho.c"
#include "../BabyStepGiantStep/BabyStepGiantStep.c"
#include "pohligHellman.c"

/** Description rapide de l'algo **
    On se place dans le cadre d'un groupe G de taille n : G = Z/nZ. A noter que n n'est pas forcément premier ici.
    Etant donné h appartient à G, le but est de retrouver g^x = h mod n, avec g un générateur de g et x appartient à [0, n-1]

    Pour cela, nous commençons par décomposer n en un produit de puissances de facteurs premiers. n = p0^(u0) * ... * pk^(uk)
    Ces facteurs étant premiers, on peut appliquer un algo tel que Pollard's Rho ou BSGS pour trouver des "sous-solutions" xi telles que
gi^xi = hi mod(pi^ui)

    gi = g^(n/(pi^ui))      hi = h^(n/(pi^ui))

    Une fois tous les xi obtenus, on applique le théorème des restes chinois.
        x = x0 mod(p0^u0)
        x = x1 mod(p1^u1)
        ...
        x = xk mod(pk^uk)

    On trouve alors une solution x étant la solution du problème initial.
**/


int main()
{
    double duree_ecoulee2;
    clock_t debut_programme2;
    debut_programme2 = clock();

    mpz_t n;
    mpz_init_set_str(n, "8263", 10);

    mpz_t h;
    mpz_init_set_str(h, "542", 10);

    mpz_t g;
    mpz_init(g);
    //findGenerator(g, n);
    mpz_set_str(g, "1693", 10);
    //gmp_printf("g = %Zd\n", g);

    mpz_t ordreg;
    mpz_init(ordreg);
    calculOrdre(ordreg, g, n);
    gmp_printf("ordreg = %Zd\n", ordreg);

    PrimeFact P[1000];
    initTab(P);
    affiche(P);

    primeFactDecomp(P, ordreg);
    affiche(P);

    int i;
    for(i=0; mpz_cmp_ui(P[i].prime, 1) != 0; i++)
    {
        mpz_t beta, alpha;
        mpz_init(beta);
        mpz_init(alpha);


        mpz_t fact;
        mpz_init_set(fact, P[i].prime);         // fact_i = (p_i)^(e_i)
        mpz_powm(fact, fact, P[i].pow, n);

        mpz_t m;
        mpz_init(m);            // m = n/(p_i)^(e_i)
        mpz_div(m, ordreg, fact);



        mpz_powm(beta, h, m, n);
        gmp_printf("\nbeta = %Zd\n\n", beta);      // beta = 542^1009

        mpz_powm(alpha, g, m, n);
        gmp_printf("\nalpha = %Zd\n\n", alpha);

        //pollardsRho(P[i].solution, n, beta, alpha);
        baby_step_giant_step(P[i].solution, n, alpha, beta);


        gmp_printf("sol = %Zd\n", P[i].solution);
    }

    mpz_t x;
    mpz_init_set_ui(x, 0);
    chineseRemainder(x, P, ordreg);

    mpz_clear(n);

    for(i=0; i<1000; i++)
    {
        mpz_clear(P[i].prime);
        mpz_clear(P[i].pow);
        mpz_clear(P[i].solution);
    }

    duree_ecoulee2 = clock() - debut_programme2;
    printf ( "\n\nTemps d'execution 2 : %f\n", duree_ecoulee2 );

    return 0;
}
