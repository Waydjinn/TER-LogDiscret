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

typedef struct
{
  mpz_t c1;
  mpz_t c2;
} C;

void key_gen(pk *pk, sk *sk, mpz_t p, mpz_t g, mpz_t y, mpz_t x)
{
/////////initialisation////////


  mpz_set_ui(pk->p, 8999);
  mpz_set_ui(pk->g, 6426);


  gmp_printf ("%s %Zd\n", "p =", pk->p);
  gmp_printf ("%s %Zd\n", "g =", pk->g);

///////Debut random
  unsigned long int seed = time(NULL);
  gmp_randstate_t state;
  gmp_randinit_default (state);
  gmp_randseed_ui(state, seed);
  mpz_urandomm (sk->x, state, pk->p);
  mpz_powm_sec(pk->y, pk->g, sk->x, pk->p);

  gmp_printf ("%s %Zd\n", "x random", sk->x);
  gmp_printf ("%s %Zd\n", "g^x =", pk->y);




  gmp_randclear(state);


}

void encryption(pk *pk,sk *sk, C *C ,mpz_t p, mpz_t g, mpz_t y, mpz_t x, mpz_t m, mpz_t c1, mpz_t c2)
{
  gmp_printf ("%s %Zd\n", "x random=", sk->x);
  gmp_printf ("%s %Zd\n", "g^x =", pk->y);

  mpz_t n;
  mpz_init(n);
//  mpz_t s;
//  mpz_init(s);
  unsigned long int seed = time(NULL);
  gmp_randstate_t state;
  gmp_randinit_default (state);
  gmp_randseed_ui(state, seed);
  mpz_urandomm (n, state, pk->p);
  mpz_powm_sec(C->c1, pk->y, n, pk->p);
  gmp_printf ("%s %Zd\n", "g^xy", C->c1);


  mpz_mul(C->c2, m, C->c1);

  gmp_printf ("%s %Zd\n", "c2", C->c2);




}

int main()
{
  pk pk;
  sk sk;
  C C;
  mpz_t p,g,x,y,m,c1,c2;
  mpz_init (pk.p);
  mpz_init (pk.g);
  mpz_init(sk.x);
  mpz_init(pk.y);
  mpz_init(m);
  mpz_set_ui(m, 204);
  mpz_init(C.c1);
  mpz_init(C.c2);


  key_gen(&pk, &sk, p, g, y, x);
  encryption(&pk, &sk,&C, p, g, y, x, m, c1, c2);

  gmp_printf ("%s %Zd\n", "c2 v2 ", C.c2);


  mpz_clear(pk.g);
  mpz_clear(sk.x);
  mpz_clear(pk.p);
  mpz_clear(pk.y);
  mpz_clear(m);
  mpz_clear(C.c1);
  mpz_clear(C.c2);

  return 0;
}
