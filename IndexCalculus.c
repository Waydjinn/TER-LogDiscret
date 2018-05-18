#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void index_calculus(/*mpz_t resultat,*/ mpz_t ordre, mpz_t generateur, mpz_t elt/*, int* factorBase*/){
	
	mpz_t res_puissance;
	mpz_init(res_puissance);
	
	int i;
	for(i=1;i<100;i++){
		mpz_powm_ui(res_puissance,generateur,i,ordre);
		gmp_printf("%Zd ^(%d) = %Zd \n",generateur,i,res_puissance);
	}
}

void recherche_exhaustive_log(mpz_t ordre){
	
	int i=0;
	
	mpz_t res, expo;
	mpz_init(res); mpz_init(expo);
	mpz_set_ui(expo,exp(1));
	mpz_set_ui(res,0);
	
	
	while(mpz_cmp(res,ordre)>0){
		mpz_pow_ui(res,expo,i);
		i++;
	}
	gmp_printf("Le résultat est %d. Pour ln(%Zd)=%Zd.\n",i,ordre,res);
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

int main(){
	mpz_t n;
	mpz_init(n);
	mpz_set_ui(n,83);
	recherche_exhaustive_log(n);
}


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


//int main(int argc, char **argv)
//{
	//mpz_t generateur,ordre,elt;
	//mpz_init(generateur);mpz_init(ordre);mpz_init(elt);
	
	//mpz_set_ui(generateur,3);
	//mpz_set_ui(ordre,1217);
	//mpz_set_ui(elt,37);
	
	//index_calculus(ordre,generateur,elt);
	
	//return 0;
//}
*/
