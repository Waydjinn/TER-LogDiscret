#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void applyFunction(mpz_t res, mpz_t x, mpz_t n)
{
    mpz_pow_ui(res, x, 2);
    mpz_add_ui(res, res, 1);
    mpz_mod(res, res, n);
}

int main()
{
    mpz_t rho, x, y, n, tmp;
    mpz_init_set_ui(x, 2);
    mpz_init_set_ui(y, 2);
    mpz_init_set_str(n, "10 023 859 281 455 311 421", 10);
    mpz_init(tmp);
    mpz_init(rho);

    do
    {
        applyFunction(x, x, n);
        gmp_printf("x = %Zd\n", x);

        applyFunction(y, y, n);
        applyFunction(y, y, n);
        gmp_printf("y = %Zd\n", y);


        mpz_sub(tmp, y, x);
        mpz_mod(tmp, tmp, n);
        mpz_abs(tmp, tmp);
        gmp_printf("tmp = %Zd\n", tmp);


        mpz_gcd(rho, tmp, n);
        mpz_mod(rho, rho, n);
        gmp_printf("rho = %Zd\n", rho);
        if (mpz_cmp_ui(rho, 1) > 0)
        {
            gmp_printf("facteur trouvé : %Zd\n", rho);
            exit(0);
        }
    }while (mpz_cmp(y, x) != 0);

    printf("*************** FAILURE ***************");

    return 0;
}
