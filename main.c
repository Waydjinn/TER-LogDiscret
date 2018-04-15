#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/** On se place dans le cadre d'un groupe G fini d'ordre n, avec n+1 un nombre premier.
Puisque n+1 est premier, un élément générateur du groupe G est 2.

**/

mpz_t N;

mpz_t n;

mpz_t alpha, beta;

void findGenerator(mpz_t gen, mpz_t n)
{
    mpz_t test;
    mpz_init_set_ui(test, 2);

    while(mpz_cmp(test, n) != 0)
    {
        mpz_t compt;
        mpz_init_set_ui(compt, 1);

        mpz_t slow, fast;
        mpz_init_set_ui(slow, 1);
        mpz_init_set_ui(fast, 1);

        do
        {
            mpz_mul(slow, slow, test);
            mpz_mod(slow, slow, n);

            mpz_mul(fast, fast, test);
            mpz_mul(fast, fast, test);
            mpz_mod(fast, fast, n);

            mpz_add_ui(compt, compt, 1);
        }while (mpz_cmp(slow, fast) != 0);

        if(mpz_cmp(compt, n) == 0)
        {
            mpz_set(gen, test);
            gmp_printf("%Zd est un générateur de Z/%ZdZ\n", test, n);
            mpz_set(test, n);
        }
        else
        {
            mpz_add_ui(test,test,1);
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
    mpz_init_set_str(N, "1019", 10);
    gmp_printf("N = %Zd\n", N);

    mpz_init(n);
    mpz_sub_ui(n, N, 1);

    //mpz_init_set_ui(alpha, 2);      // alpha est est générateur de G
    findGenerator(alpha, N);
    mpz_init_set_ui(beta, 1014);       // beta est un élément de G

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
            //mpz_abs(tmp1,tmp1);
            mpz_sub(tmp2, a, A);
            //mpz_abs(tmp2, tmp2);
            if(mpz_cmp_ui(tmp1, 0) == 0)
            {
                printf("\n\ndivision par 0 impossible\n\n");
                break;
            }
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

    return 0;
}
