#include <gmp.h>
#include <stdio.h>
#include <time.h>

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

void key_gen(pk *pk)
{
/////////initialisation////////
  mpz_t p;
  mpz_init (pk->p);
  mpz_t g;
  mpz_init (pk->g);
  mpz_set_ui(pk->p, 6883);
  mpz_set_ui(pk->g, 4344);
  mpz_t n;
  mpz_init(n);
  mpz_t res;
  mpz_init(res);
  mpz_t y;
  mpz_init(pk->y);

  mpz_set(n, pk->p);

  gmp_printf ("%s %Zd\n", "p =", pk->p);
  gmp_printf ("%s %Zd\n", "g =", pk->g);

///////Debut random
  unsigned long int seed = time(NULL);
  gmp_randstate_t state;
  gmp_randinit_default (state);
  gmp_randseed_ui(state, seed);
  mpz_urandomm (res, state, n);
  mpz_powm_sec(pk->y, pk->g, res, pk->p);
  gmp_printf ("%s %Zd\n", "res =", res);

  gmp_printf ("%s %Zd\n", "g^x mod p =", pk->y);



  //mpz_clear(seed);
  mpz_clear(n);
  mpz_clear(pk->g);
  gmp_randclear(state);
  mpz_clear(res);
  mpz_clear(pk->p);
  mpz_clear(pk->y);

}

int main()
{
  pk pk;
  key_gen(&pk);

  return 0;
}
