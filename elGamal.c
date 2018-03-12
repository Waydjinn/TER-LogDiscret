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

void key_gen(pk *pk;)
{
  gmp_printf ("%s is an mpz\n", "here");

  //pk *pk;

  gmp_printf ("%s is an mpz\n", "here");

  mpz_set_ui(pk pk.p, 6883);
  mpz_set_ui(pk->g, 4344);

  gmp_printf ("%s is an mpz\n", "here");

  gmp_randstate_t state;
  mpz_t n;
  mpz_init(n);
  mpz_t res;
  mpz_init(res);
  mpz_set(n, pk->p);
  gmp_randinit_default (state);
  mpz_urandomm (res, state, n);


  mpz_powm_sec(pk->y, pk->g, res, pk->p);
  gmp_printf ("%s is an mpz %Zd\n", "here", pk->y);


  mpz_clear(n);
  mpz_clear(pk->g);
  mpz_clear(res);
  mpz_clear(pk->p);
  mpz_clear(pk->y);


}

int main()
{
  key_gen();

  return 0;
}
