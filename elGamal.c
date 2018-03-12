#include <gmp.h>
#include <stdio.h>

typedef struct
{
  mpz_t p;	    /* nb premier */
  mpz_t g;	    /* generateur du groupe */
  mpz_t y;	    /* g^x mod p */
} pk;


typedef struct
{
  mpz_t p;	    /* nb premier */
  mpz_t g;	    /* generateur du groupe */
  mpz_t y;	    /* g^x mod p */
  mpz_t x;	    /* exposant secret */
} sk;

void key_gen()
{
  gmp_randstate_t state;
  mpz_t n;
  mpz_init(n);
  mpz_t res;
  mpz_init(res);
  mpz_set_ui(n, 10);
  gmp_randinit_default (state);
  mpz_urandomm (res, state, n);

  mpz_clear(n);


}

int main()
{
  return 0;
}
