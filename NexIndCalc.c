#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void index_calculus(/*mpz_t resultat,*/ mpz_t ordre, mpz_t generateur /*mpz_t elt*//*, int* factorBase*/){

mpz_t res_puissance;
mpz_init(res_puissance);

int i;
for(i=1;i<100;i++){
mpz_powm_ui(res_puissance,generateur,i,ordre);
gmp_printf("%Zd ^(%d) = %Zd \n",generateur,i,res_puissance);
}
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
mpz_t generateur,ordre,elt;
mpz_init(generateur);mpz_init(ordre);mpz_init(elt);

mpz_set_ui(generateur,3);
mpz_set_ui(ordre,1217);
mpz_set_ui(elt,37);

//index_calculus(ordre,generateur);

//decomposition(ordre);

recherche_exhaustive_log_INT(ordre);

return 0;
}
