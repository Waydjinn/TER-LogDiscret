#include <gmp.h>
#include <stdio.h>
#include <time.h>

typedef struct
{
  mpz_t p;	    /* nb premier */
  mpz_t g;	    /* generateur du groupe */
  mpz_t h;	    /* g^x mod p */
} pk;


typedef struct
{
  mpz_t p;	    /* nb premier */
  mpz_t g;	    /* generateur du groupe */
  mpz_t h;	    /* g^x mod p */
  mpz_t x;	    /* exposant secret */
} sk;

typedef struct
{
  mpz_t c1;
  mpz_t c2;
} C;

void key_gen(pk *pk, sk *sk, mpz_t p, mpz_t g, mpz_t h, mpz_t x)
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
  mpz_powm_sec(pk->h, pk->g, sk->x, pk->p);

  gmp_printf ("%s %Zd\n", "x random", sk->x);
  gmp_printf ("%s %Zd\n", "g^x =", pk->h);

  gmp_randclear(state);

}

void encryption(pk *pk, sk *sk, C *C ,mpz_t p, mpz_t g, mpz_t h, mpz_t x, mpz_t m, mpz_t c1, mpz_t c2)
{


  mpz_t y;
  mpz_init(y);
  mpz_t s;
  mpz_init(s);
  mpz_t tmp;
  mpz_init(tmp);
  unsigned long int seed = time(NULL);
  gmp_randstate_t state;
  gmp_randinit_default (state);
  gmp_randseed_ui(state, seed);
  mpz_urandomm (y, state, pk->p);

  mpz_mul(tmp, sk->x, y); //tmp prend la valeur x*y
  mpz_powm_sec(C->c1, pk->g, y, pk->p); //c1 = g^y
  mpz_powm_sec(s, pk->g, tmp, pk->p); //s = g^x*y

  gmp_printf ("%s %Zd\n", "g^xy", s);
  gmp_printf ("%s %Zd\n", "g^xy", tmp);

  gmp_printf ("%s %Zd\n", "g^xy", y);

  mpz_mul(C->c2, m, s);
  gmp_printf ("%s %Zd\n", "c1", C->c1);
  gmp_printf ("%s %Zd\n", "c2", C->c2);

  mpz_clear(y);

}

void decryption(sk *sk, pk *pk, C *C, mpz_t p, mpz_t g, mpz_t h, mpz_t x, mpz_t m, mpz_t c1, mpz_t c2)
{
  // Square And Multiply

  // Euclide Etendu



}

int main()
{
  pk pk;
  sk sk;
  C C;
  mpz_t m;
  gmp_printf ("%s %Zd\n", "juste p", pk.p);

  mpz_init (pk.p);
  gmp_printf ("%s %Zd\n", "juste p", pk.p);

  mpz_init (pk.g);
  mpz_init(sk.x);
  mpz_init(pk.h);
  mpz_init(m);
  mpz_set_ui(m, 204);
  mpz_init(C.c1);
  mpz_init(C.c2);

  key_gen(&pk, &sk, pk.p, pk.g, pk.h, sk.x);
  encryption(&pk, &sk,&C, pk.p, pk.g, pk.h, sk.x, m, C.c1, C.c2);
  decryption(&sk, &pk,&C, pk.p, pk.g, pk.h, sk.x, m, C.c1, C.c2);
  gmp_printf ("%s %Zd\n", "c2 v2 ", C.c2);


  mpz_clear(pk.g);
  mpz_clear(sk.x);
  mpz_clear(pk.p);
  mpz_clear(pk.h);
  mpz_clear(m);
  mpz_clear(C.c1);
  mpz_clear(C.c2);

  return 0;
}
