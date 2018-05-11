#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

void initVars(mpz_t x, mpz_t X, mpz_t a, mpz_t A, mpz_t b, mpz_t B) {
    mpz_init(x);
    mpz_init(X);
    mpz_init(a);
    mpz_init(A);
    mpz_init(b);
    mpz_init(B);
}

void resetVars(mpz_t x, mpz_t X, mpz_t a, mpz_t A, mpz_t b, mpz_t B) {
    mpz_set_str(x, "1", 10);
    mpz_set_str(X, "1", 10);
    mpz_set_str(a, "0", 10);
    mpz_set_str(A, "0", 10);
    mpz_set_str(b, "0", 10);
    mpz_set_str(B, "0", 10);
}

void applyFunction(mpz_t x1, mpz_t a1, mpz_t b1, mpz_t n1, mpz_t g1, mpz_t h1) {
    mpz_t tmp, ordre;
    mpz_init(tmp);
    mpz_init(ordre);
    mpz_mod_ui(tmp, x1, 3);
    mpz_sub_ui(ordre, n1, 1);


    if (mpz_cmp_ui(tmp, 0) == 0) {
        mpz_powm_ui(x1, x1, 2, n1);
        mpz_mul_ui(a1, a1, 2);
        mpz_mod(a1, a1, ordre);
        mpz_mul_ui(b1, b1, 2);
        mpz_mod(b1, b1, ordre);
    }

    if (mpz_cmp_ui(tmp, 1) == 0) {
        mpz_mul(x1, x1, g1);
        mpz_mod(x1, x1, n1);
        mpz_add_ui(a1, a1, 1);
        mpz_mod(a1, a1, ordre);
    }

    if (mpz_cmp_ui(tmp, 2) == 0) {
        mpz_mul(x1, x1, h1);
        mpz_mod(x1, x1, n1);
        mpz_add_ui(b1, b1, 1);
        mpz_mod(b1, b1, ordre);
    }

    mpz_clear(tmp);
    mpz_clear(ordre);
}

void clearVars(mpz_t x, mpz_t X, mpz_t a, mpz_t A, mpz_t b, mpz_t B, mpz_t n, mpz_t g, mpz_t h, mpz_t i) {
    mpz_clear(x);
    mpz_clear(X);
    mpz_clear(a);
    mpz_clear(A);
    mpz_clear(b);
    mpz_clear(B);
    mpz_clear(n);
    mpz_clear(g);
    mpz_clear(h);
    mpz_clear(i);
}

void verif(mpz_t g, mpz_t res, mpz_t h, mpz_t n) {
    mpz_t verif;
    mpz_init(verif);

    mpz_powm(verif, g, res, n);
    if (mpz_cmp(verif, h) == 0) {
        printf("\nCe resultat est bien verifie\n\n");
    }

    mpz_clear(verif);
}

int main()
{
    mpz_t n, g, h;
    mpz_init_set_str(n, "1009", 10);
    mpz_init_set_str(g, "669", 10);
    mpz_init_set_str(h, "542", 10);

    mpz_t x, X, a, A, b, B;
    initVars(x, X, a, A, b, B);
    resetVars(x, X, a, A, b, B);

    mpz_t i;
    mpz_init_set_ui(i, 0);
    while (mpz_cmp(i, n) != 0) {
        applyFunction(x, a, b, n, g, h);

        applyFunction(X, A, B, n, g, h);
        applyFunction(X, A, B, n, g, h);

        gmp_printf("%Zd : %Zd          %Zd          %Zd          %Zd          %Zd          %Zd\n", i, x, a, b, X, A, B);

        if (mpz_cmp(x, X) == 0) {
            break;
        }

        mpz_add_ui(i, i, 1);
    }

    mpz_t BMoinsb, aMoinsA, res, ordre;
    mpz_init(BMoinsb);
    mpz_init(aMoinsA);
    mpz_init(res);
    mpz_init(ordre);
    mpz_sub(BMoinsb, B, b);
    mpz_sub(aMoinsA, a, A);
    mpz_sub_ui(ordre, n, 1);

    if (mpz_cmp_ui(BMoinsb, 0) == 0) {
        printf("\nImpossible de trouver une solution. (Division par 0)\n\n");
    }
    else {
        mpz_invert(BMoinsb, BMoinsb, n);
        mpz_mul(res, BMoinsb, aMoinsA);
        mpz_mod(res, res, n);

        gmp_printf("\nUne solution trouvée : %Zd\n\n", res);

        verif(g, res, h, n);
    }

    clearVars(x, X, a, A, b, B, n, g, h, i);
    mpz_clear(BMoinsb);
    mpz_clear(aMoinsA);
    mpz_clear(res);
    mpz_clear(ordre);

    return 0;
}
