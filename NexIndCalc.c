#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct PrimeFact
{
    mpz_t prime;
    mpz_t pow;
    mpz_t solution;
}PrimeFact;

void initTab(PrimeFact* P)
{
    int i;
    for (i=0; i<1000; i++)
    {
        mpz_init(P[i].prime);
        mpz_init(P[i].pow);
        mpz_init(P[i].solution);
        mpz_set_ui(P[i].prime, 1);
        mpz_set_ui(P[i].pow, 0);
    }
}

void affiche(PrimeFact* P)
{
    int i=0;
    while(mpz_cmp_ui(P[i].prime, 1) != 0)
    {
        gmp_printf("Facteur : %Zd\nPuissance : %Zd\n\n", P[i].prime, P[i].pow);
        i++;
    }
}

void primeFactDecomp(PrimeFact* P, mpz_t n)      // décomposition de n en facteurs premiers
{
    mpz_t copyn;
    mpz_init_set(copyn, n);

    mpz_t rootn;
    mpz_init(rootn);
    mpz_sqrt(rootn, copyn);
    mpz_add_ui(rootn, rootn, 1);

    mpz_t i;
    mpz_init_set_ui(i, 2);

    mpz_t puis;
    mpz_init_set_ui(puis, 0);

    int indice = 0;

    while (mpz_cmp(i, rootn) <= 0)          // Tant que les facteurs testés sont inférieurs ou égaux à la racine de n
    {
        if(mpz_probab_prime_p(copyn, 40) > 0)
        {                   // Si n est premier ou supposé premier après 40 tests on le renvoie
            mpz_init_set(P[indice].prime, copyn);
            mpz_init_set_ui(P[indice].pow, 1);
            indice++;

            break;
        }

        mpz_t test;
        mpz_init(test);
        do
        {
            mpz_mod(test, copyn, i);
            if(mpz_cmp_ui(test, 0) == 0)        // Si i divise n
            {
                mpz_div(copyn, copyn, i);
                mpz_add_ui(puis, puis, 1);
            }
        }while (mpz_cmp_ui(test, 0) == 0);

        if (mpz_cmp_ui(puis, 0) != 0)       // On n'affiche pas les puissance à 0.
        {
            mpz_set(P[indice].prime, i);
            mpz_set(P[indice].pow, puis);
            indice++;
        }
        mpz_nextprime(i, i);
        mpz_set_ui(puis, 0);

        mpz_clear(test);
    }
}

void putInMatrix(PrimeFact* P, mpz_t* FactorBase, mpz_t tailleFBase)
{

}

int testInFactorBase(PrimeFact* P, mpz_t* FactorBase, mpz_t tailleFBase, mpz_t primeVoulu)
{
    /* AMELIORABLE ************************ */

    int dernPos=0;

    //Parcours de tous les premiers de la décomposition
    int i=0;
    int j=0;

    int primeVouluOK = -1;

    while(mpz_cmp_ui(P[i].prime, 1) != 0)
    {
        //gmp_printf("Facteur : %Zd\nPuissance : %Zd\n\n", P[i].prime, P[i].pow);

        if(mpz_cmp(P[i].prime, primeVoulu) == 0)
        {
            primeVouluOK = 0;
        }
        
        //On parcours la Factor Base (à partir de la dernière position repérée)
        for (j = dernPos; j <= mpz_get_ui(tailleFBase); j++)
        {
            //Si après avoir tout parcouru on le l'a pas trouvé
            if( mpz_get_ui(tailleFBase) == j)
            {
                return -1;//on prévient que ce n'est pas ok
            }

            //Si on le trouve dans la Factor Base
            //if( mpz_cmp(mpz_get_ui(P[i].prime), mpz_get_ui(FactorBase[j]) ) == 0)
            if( mpz_cmp(P[i].prime, FactorBase[j] ) == 0)
            {
                //S'il n'y a qu'un seul prime dans la décomposition, on prévient avec un 2.
                if(mpz_cmp_ui(P[1].prime, 1) == 0)
                {
                   return 2; 
                }

                dernPos=j; break; //On retient la position
            }   
        }

        i++;
    }

    //Si on arrive ici, tous les nombres sont bien dans la FBase
    if(primeVouluOK == 0)
    {
        return 11;
    }
    else if(primeVouluOK == -1)
    {
        return 22;
    }

     
}


void index_calculus(mpz_t ordre, mpz_t generateur, mpz_t elt, mpz_t tailleFBase)
{

    /*mpz_t res_puissance;
    mpz_init(res_puissance);

    int i;
    for(i=1;i<100;i++){
    mpz_powm_ui(res_puissance,generateur,i,ordre);
    gmp_printf("%Zd ^(%d) = %Zd \n",generateur,i,res_puissance);
    }*/

   /*  Création matrice ******************************* */
   /*mpz_t* FactorBase = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));

    int x;
    for(x=0; x<mpz_get_ui(tailleFBase); x++) //Init
    {
        mpz_init(FactorBase[x]);
    }

    mpz_t prime;
    mpz_init(prime);
    mpz_set_ui(prime,2);

    for(x=0; x<mpz_get_ui(tailleFBase); x++) //
    {
        mpz_set(FactorBase[x], prime);
        mpz_nextprime (prime, prime);
    }

    gmp_printf("Factor base : ", FactorBase[x]);
    for(x=0; x<mpz_get_ui(tailleFBase); x++) //
    {
        gmp_printf("%Zd, ", FactorBase[x]);
    }
    printf("\n");*/


    /* 1. Factor base Init ******************************* 
     On veut avoir tous les nombres premiers de 2 à celui correspondant à 
    la taille du FBase, pour comparer par la suite plus simplement. */

    mpz_t* FactorBase = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));

    int x;
    for(x=0; x<mpz_get_ui(tailleFBase); x++) //Init
    {
        mpz_init(FactorBase[x]);
    }

    mpz_t prime;
    mpz_init(prime);
    mpz_set_ui(prime,2);

    for(x=0; x<mpz_get_ui(tailleFBase); x++) //
    {
        mpz_set(FactorBase[x], prime);
        mpz_nextprime (prime, prime);
    }

    gmp_printf("Factor base : ", FactorBase[x]);
    for(x=0; x<mpz_get_ui(tailleFBase); x++) //
    {
        gmp_printf("%Zd, ", FactorBase[x]);
    }
    printf("\n");

    /* 2. Calcul des équations  ************************************* */
    /* On calcule maintenant plein de g^k mod p. 
    On compare à chaque fois si la décomposition du nombre obtenu est faite de nombres
    du Factor base.
    Une fois qu'on a autant d'équations que la taille du FBase, on arrête.
    Aussi, on ne compte pas tous les g^k < p du début.

    On les place aussi dans la matrice.*/

    PrimeFact P[1000];
    initTab(P);

    int cnt=0;
    int test=1;
    int k=1;
    mpz_t resPowm;
    mpz_init(resPowm);

    //PROBLEME : JE CROIS QU'ON VEUT GARDER 3^1 = 3, A RESOUDRE !!! (facile)

    //Tant que g^k < p (sans mod), on continue
    while( mpz_cmp(  resPowm, ordre )  < 0   )
    {
        mpz_pow_ui(resPowm, generateur, k);
        k++;
    }

    mpz_t primeVoulu;
    mpz_init(primeVoulu);
    mpz_set_ui(primeVoulu,2);

    //while(cnt < mpz_get_ui(tailleFBase))
    while(cnt < 40)
    {
        mpz_powm_ui (resPowm, generateur, k, ordre);
        
        primeFactDecomp(P, resPowm);
        

        test = testInFactorBase(P, FactorBase, tailleFBase, primeVoulu);
        if (test == 11)
        {
            //printf("-> Dans la FBase !!\n");
            printf("-> Sélectionnée !\n");
             cnt++;
             gmp_printf("\nN°: %d - k : %d Nombre : %Zd\n", cnt, k, resPowm);
             affiche(P);

             mpz_nextprime (primeVoulu, primeVoulu);

             gmp_printf("test %Zd" , FactorBase[mpz_get_ui(tailleFBase)-1]);

             if(mpz_cmp(primeVoulu, FactorBase[mpz_get_ui(tailleFBase)-1]) > 0)
             {
                cnt = 40;
             }
        }
        else if (test == 22)
        {
            //printf("Nope\n");
        }
        /*else if (test == 2)
        {
            printf("Un seul facteur !");
            cnt++;
            gmp_printf("\nN°: %d - k : %d Nombre : %Zd\n", cnt, k, resPowm);
            affiche(P);
        }*/
        else if (test == -1)
        {
            //printf("Nope\n");
        }

        k++;
       
    }





    free(FactorBase);

}

void decomposition(mpz_t nombre)
{
    mpz_t nextPrime;
    mpz_init(nextPrime);
    mpz_set_ui(nextPrime,2);

    mpz_t temp;
    mpz_init(temp);

    int * cntTab;

    while(nombre > 0)
    {   while(mpz_get_ui(nombre)%mpz_get_ui(nextPrime)==0)
        {
            mpz_cdiv_r (nombre, temp, nextPrime);
            //ctn++;
        }
        mpz_nextprime(nextPrime,nextPrime);
    }

    //gmp_printf("La décomposition est : %");
}

/*void recherche_exhaustive_log(mpf_t ordre){

int i=0;

mpf_t res, expo;

mpf_init(res); mpf_init(expo);
double e = exp(1);
//double e = 2.71828182845904523536;
mpf_set_d(expo, e);
gmp_printf("exp(1): %lf gmp :%Zd", e, expo);
mpf_set_ui(res,0);


while(mpf_cmp(res,ordre)>0){
mpf_pow_ui(res,expo,i);
i++;
}
gmp_printf("Le résultat est %d. Pour ln(%Zlf)=%Zlf.\n",i,ordre,res);
}*/

void recherche_exhaustive_log_INT(mpz_t ordre){

    int i=0;

    double res;
    res = mpz_get_d(ordre);

    log(res);
    printf("Log est %lf\n", res);

    

    //log(5675675);

}



void elimination_Gauss(int** tab,mpz_t taille){
int lenght = mpz_get_ui(taille);
int r=0;
int j;
for(j=0;j<lenght;j++){

}
}

void calcul_tailleFactorBase(mpz_t res, mpz_t ordre){
mpz_t calcul;
mpz_init(calcul);



}

/*int main(){
mpf_t n;
mpf_init(n);
mpf_set_d(n,83);
recherche_exhaustive_log(n);
}*/


/*int main()
{
int i,j,k,n;
float A[20][20],c,x[10];
printf("\nEnter the size of matrix: ");
scanf("%d",&n);
printf("\nEnter the elements of augmented matrix row-wise:\n");
for(i=1; i<=n; i++)
{
    for(j=1; j<=(n+1); j++)
    {
        printf(" A[%d][%d]:", i,j);
        scanf("%f",&A[i][j]);
    }
}

for(j=1; j<=n; j++)
{
    for(i=1; i<=n; i++)
    {
        if(i!=j)
        {
            c=A[i][j]/A[j][j];
            for(k=1; k<=n+1; k++)
            {
                A[i][k]=A[i][k]-c*A[j][k];
            }
        }
    }
}
printf("\nThe solution is:\n");
for(i=1; i<=n; i++)
{
    x[i]=A[i][n+1]/A[i][i];
    printf("\n x%d=%f\n",i,x[i]);
}
return(0);
}

*/
int main(int argc, char **argv)
{
    /* INITS ********************************************** */
    mpz_t generateur,ordre,elt,tailleFBase;
    mpz_init(generateur);mpz_init(ordre);mpz_init(elt);mpz_init(tailleFBase);

    //Example "rouge"
   /* mpz_set_ui(generateur,3);
    mpz_set_ui(ordre,1217);
    mpz_set_ui(elt,37);
    mpz_set_ui(tailleFBase, 8);//HARDCODé ICI POUR L'INSTANT*/

    //Exemple "bleu"
    mpz_set_ui(generateur,2);
    mpz_set_ui(ordre,83);
    mpz_set_ui(elt,31);
    mpz_set_ui(tailleFBase, 4);//HARDCODé ICI POUR L'INSTANT

    

    index_calculus(ordre, generateur, elt, tailleFBase);

   // primeFactDecomp(P, ordre);
    //affiche(P);

    //mpz_mul_ui(elt, 5795825882885517379, 2);
    //gmp_printf("res = %Zd", elt);

//recherche_exhaustive_log_INT(ordre);

return 0;
}
