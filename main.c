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

mpz_t n;

mpz_t alpha, beta;

int calcul_mod(int base, int mod)
{
    int res, n;
    n = 1;
    res = base%mod;
    while (res != 1 && res != mod-1)
    {
        res = (res*base)%mod;
        n = n+1;
    }
    printf("n = %d\nres = %d", n, res);
    return n;
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
                                                // compt++
            mpz_add_ui(compt, compt, 1);
        }while (mpz_cmp(slow, fast) != 0);      // Tant qu'un cycle n'a pas été trouvé

        if(mpz_cmp(compt, n) == 0)      // si n éléments ont été parcouru avant d'atteindre un cycle ( compt = n )
        {
            mpz_set(gen, test);
            gmp_printf("%Zd est un générateur de Z/%ZdZ\n", test, n);
            mpz_set(test, n);
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

    if (mpz_cmp_ui(tmp, 0) == 0)
    {
        mpz_pow_ui(x, x, 2);
        mpz_mod(x, x, N);

        mpz_mul_ui(a, a, 2);
        mpz_mod(a, a, n);

        mpz_mul_ui(b, b, 2);
        mpz_mod(b, b, n);
    }
    if (mpz_cmp_ui(tmp, 1) == 0)
    {
        mpz_mul(x, x, alpha);
        mpz_mod(x, x, N);

        mpz_add_ui(a, a, 1);
        mpz_mod(a, a, n);
    }
    if (mpz_cmp_ui(tmp, 2) == 0)
    {
        mpz_mul(x, x, beta);
        mpz_mod(x, x, N);

        mpz_add_ui(b, b, 1);
        mpz_mod(b, b, n);
    }

    mpz_clear(tmp);
}

int main()
{
    mpz_init_set_str(N, "1009", 10);
    gmp_printf("N = %Zd\n", N);

    mpz_init(n);
    mpz_sub_ui(n, N, 1);

    //mpz_init_set_ui(alpha, 2);      // alpha est est générateur de G
    mpz_init(alpha);
    findGenerator(alpha, N);
    mpz_init_set_ui(beta, 3);       // beta est un élément de G

    mpz_t x, X, a, A, b, B;         // Les variables majuscules vont évouluer 2X plus vite que celles en minuscules
    mpz_init_set_ui(a, 0);
    mpz_init_set_ui(b, 0);
    mpz_init_set_ui(A, 0);
    mpz_init_set_ui(B, 0);
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(X, 1);

    mpz_t i;
    mpz_init_set_ui(i, 1);
    while(mpz_cmp(i, N) < 0)
    {
        applyFunction(x, a, b);
        gmp_printf("x = %Zd\n", x);

        applyFunction(X, A, B);
        applyFunction(X, A, B);
        gmp_printf("X = %Zd\n", X);

        if (mpz_cmp(x, X) == 0)
        {
            gmp_printf("x : %Zd\na : %Zd\nb : %Zd\nX : %Zd\nA : %Zd\nB : %Zd\n ", x, a, b, X, A, B);

            mpz_t gamma, tmp1, tmp2;
            mpz_init(gamma);
            mpz_init(tmp1);
            mpz_init(tmp2);

            mpz_sub(tmp1, B, b);
            mpz_mod(tmp1, tmp1, N);
            //mpz_abs(tmp1, tmp1);
            mpz_sub(tmp2, a, A);
            mpz_mod(tmp2, tmp2, N);
            //mpz_abs(tmp2, tmp2);
            if(mpz_cmp_ui(tmp1, 0) == 0)
            {
                printf("\n\ndivision par 0 impossible\n\n");
                break;
            }
            gmp_printf("%Zd, %Zd\n", tmp2, tmp1);
            mpz_cdiv_q(gamma, tmp2, tmp1);
            mpz_mod(gamma, gamma, n);

            gmp_printf("\nalpha = %Zd\n\n", alpha);
            gmp_printf("\ngamma = %Zd\n\n", gamma);

            mpz_clear(gamma);
            mpz_clear(tmp1);
            mpz_clear(tmp2);

            exit(0);
        }

        mpz_add_ui(i, i, 1);
    }
    printf("********** FAILURE **********\n");

    mpz_clear(x);
    mpz_clear(X);
    mpz_clear(a);
    mpz_clear(A);
    mpz_clear(b);
    mpz_clear(B);
    mpz_clear(n);
    mpz_clear(N);
    mpz_clear(alpha);
    mpz_clear(beta);
    mpz_clear(i);

    /*printf("\n\n\n\n");
    int n;
    n = calcul_mod(11, 1009);*/

    return 0;
}
