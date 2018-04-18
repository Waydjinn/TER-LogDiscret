#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/** On se place dans le cadre d'un groupe G fini d'ordre n, avec n+1 un nombre premier.
Le but de cet algo est de trouver un gamma tel que alpha^gamma = beta.
Pour cela, on va dire que chaque nombre peut s'écrire sous la forme de alpha^a * beta^b.
On va utiliser l'algo de détection de cycle de Floyd, c'est à dire qu'on va parcourir deux fois la boucle, dont une deux
fois plus vite que l'autre. Ces exécutions se font en parallèles.
l'algo va faire x = f(x) à chaque étape (ou X = f(f(x)) pour le plus rapide) et si à un moment on a x = X, alors on a un cycle et on a
un nombre égal à x et X soit:

    x = X => alpha^a * beta^b = alpha^A * beta^B
            => alpha^(a-A) = beta^(B-b)

Comme on cherche un gamma tel que alpha^gamma = beta, on obtient l'équation suivante :

    alpha^(a-A) = (alpha^gamma)^(B-b)
                = alpha^((B-b)*gamma)

Il nous suffit alors de résoudre a-A = (B-b)*gamma  <=>  gamma = (a-A)/(B-b)

**/

mpz_t N;

mpz_t ordreAlpha;

mpz_t alpha, beta;

void calculOrdre(mpz_t res, mpz_t base, mpz_t mod)
{
    int n;
    n = 1;
    mpz_mod(res, base, mod);
    while(mpz_cmp_ui(res, 1) != 0)
    {
        mpz_mul(res, res, base);
        mpz_mod(res, res, mod);
        n = n+1;
    }
    mpz_set_ui(res, n);
}

void findGenerator(mpz_t gen, mpz_t n)
{
    mpz_t test;             // test sera le nombre dont on va tester si il est générateur ou non
    mpz_init_set_ui(test, 2);

    while(mpz_cmp(test, n) != 0)        // tant que test != n
    {
        mpz_t compt;                    // compt est un compteur qui s'incrémentera à chaque vérification
        mpz_init_set_ui(compt, 1);

        mpz_t slow, fast;           // On applique l'algo de détection de Floyd, slow est la valeur prise par l'algo "lent" et fast par le "rapide"
        mpz_init_set_ui(slow, 1);
        mpz_init_set_ui(fast, 1);

        do
        {
            mpz_mul(slow, slow, test);          // slow = (slow * test) mod [n]
            mpz_mod(slow, slow, n);

            mpz_mul(fast, fast, test);          // fast = ((fast * test)* test) mod[n]
            mpz_mul(fast, fast, test);
            mpz_mod(fast, fast, n);

            mpz_add_ui(compt, compt, 1);        // compt++

        }while(mpz_cmp(slow, fast) != 0);      // Tant qu'un cycle n'a pas été trouvé

        if(mpz_cmp(compt, n) == 0)      // si n éléments ont été parcouru avant d'atteindre un cycle ( compt = n )
        {
            mpz_set(gen, test);
            gmp_printf("%Zd est un generateur de Z/%ZdZ\n", test, n);
            mpz_set(test, n);       // On casse la boucle
        }
        else        // Si on a atteint un cycle sans obtenir tous les éléments de Z/nZ
        {
            mpz_add_ui(test,test,1);    // test++
        }

        mpz_clear(compt);
        mpz_clear(slow);
        mpz_clear(fast);
    }

    mpz_clear(test);
}

void applyFunction(mpz_t x, mpz_t a, mpz_t b)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod_ui(tmp, x, 3);
                            ////////////////////////////////////////////////////
    if (mpz_cmp_ui(tmp, 0) == 0)    // x = 0 mod[3]
    {
        mpz_pow_ui(x, x, 2);        // x = x*x mod[N]
        mpz_mod(x, x, N);

        mpz_mul_ui(a, a, 2);        // a= a*2 mod[n]
        mpz_mod(a, a, ordreAlpha);

        mpz_mul_ui(b, b, 2);        // b = b*2 mod[n]
        mpz_mod(b, b, ordreAlpha);
    }                       ////////////////////////////////////////////////////

    if (mpz_cmp_ui(tmp, 1) == 0)    // sinon si x = 1 mod[3]
    {
        mpz_mul(x, x, alpha);       // x = x*alpha mod[N]
        mpz_mod(x, x, N);

        mpz_add_ui(a, a, 1);        // a = a+1 mod[n]
        mpz_mod(a, a, ordreAlpha);
    }                          ////////////////////////////////////////////////

    if (mpz_cmp_ui(tmp, 2) == 0)    // sinon, si x = 2 mod[3]
    {
        mpz_mul(x, x, beta);        // x= x*beta mod[N]
        mpz_mod(x, x, N);

        mpz_add_ui(b, b, 1);        // b = b+1 mod[n]
        mpz_mod(b, b, ordreAlpha);
    }

    mpz_clear(tmp);
}

int main()
{
    mpz_init_set_str(N, "11", 10);        // N est le nombre premier qui forme Z/NZ
    mpz_init_set_str(beta, "", 10);       // beta est un élément de G (celui dont on veut résoudre le PLD)

    mpz_init(alpha);        // alpha est un élément générateur de Z/nZ
    mpz_init(ordreAlpha);
    findGenerator(alpha, N);
    //mpz_set_ui(alpha, 6);
    calculOrdre(ordreAlpha, alpha, N);       // n = N-1

    mpz_t x, X, a, A, b, B;         // Les variables majuscules vont évouluer 2X plus vite que celles en minuscules
    mpz_init_set_ui(a, 0);
    mpz_init_set_ui(b, 0);
    mpz_init_set_ui(A, 0);
    mpz_init_set_ui(B, 0);
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(X, 1);

    mpz_t i;                    // i = 1
    mpz_init_set_ui(i, 1);
    while(mpz_cmp(i, N) < 0)        // tant que i < N
    {
        applyFunction(x, a, b);         // f(x)
        gmp_printf("x = %Zd,  ", x);

        applyFunction(X, A, B);         // f(f(x))
        applyFunction(X, A, B);
        gmp_printf("X = %Zd\n", X);

        if (mpz_cmp(x, X) == 0)     // Si on a trouvé un cycle (x = X)
        {
            gmp_printf("x : %Zd\na : %Zd\nb : %Zd\nX : %Zd\nA : %Zd\nB : %Zd\n ", x, a, b, X, A, B);

            mpz_t tmp1, tmp2;
            mpz_init(tmp1);
            mpz_init(tmp2);

            mpz_sub(tmp1, B, b);            // tmp1 = (B-b) mod[N]
            mpz_mod(tmp1, tmp1, ordreAlpha);

            if(mpz_cmp_ui(tmp1, 0) == 0)    // exception div par 0
            {
                printf("********** FAILURE **********\n");
                printf("\ndivision par 0 impossible\n\n");
                break;
            }

            mpz_sub(tmp2, a, A);            // a = (a-A) mod[N]
            mpz_mod(tmp2, tmp2, ordreAlpha);

            gmp_printf("%Zd, %Zd\n", tmp1, tmp2);

            mpz_t gamma;
            mpz_init(gamma);

            mpz_invert(tmp1, tmp1, N);
            gmp_printf("inverse tmp1 : %Zd \n", tmp1);
            mpz_mul(gamma, tmp2, tmp1);

            mpz_mod(gamma, gamma, N);

            gmp_printf("\nalpha = %Zd\n\n", alpha);
            gmp_printf("\ngamma = %Zd\n\n", gamma);

            mpz_clear(gamma);

            break;
        }

        mpz_add_ui(i, i, 1);
    }

    mpz_clear(x);
    mpz_clear(X);
    mpz_clear(a);
    mpz_clear(A);
    mpz_clear(b);
    mpz_clear(B);
    mpz_clear(ordreAlpha);
    mpz_clear(N);
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_clear(i);

    printf("\n\n");

    return 0;
}
