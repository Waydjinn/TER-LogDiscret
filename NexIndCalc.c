#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*GAUSS */
//Affichage du système
void affiche_systeme(mpz_t** A ,mpz_t* b ,mpz_t n)
{
    int i , j ;
    printf(" ===>Affichage du systeme : \n\n\n");
    
    for(i = 0 ; i < mpz_get_ui(n) ; i++)
    {
        printf("  (");
        for(j = 0 ; j < mpz_get_ui(n) ; j++)
        {
            gmp_printf("  %Zd  ",A[i][j]);
        }
        printf(" )    (X%d)   =",i+1);
        gmp_printf("\t%Zd",b[i]);
        printf("\n\n");
    }
}

/* Affichage de la solution du système */
void affiche_solution(mpz_t *x, mpz_t n)
{
    int i ;
    printf(" ===>Affichage de la solution : \n\n\n");
    
    for(i = 0 ; i < mpz_get_ui(n) ; i++)
    {
        printf("(X%d)   =",i+1);
        gmp_printf("\t%Zd",x[i]);
        printf("\n\n");
    }
}

/* Méthode d'élimination de Gauss */
void gauss(mpz_t** A, mpz_t* b, mpz_t* x, mpz_t n)//++
{
     affiche_systeme(A,b,n);//++
     
     int i, j, k ;
     int imin ;
     //float p ;
    /* float sum valmin, tump1, tump2 */
     
     mpz_t valmin;//++
     mpz_init(valmin);//++
     mpz_t tump1;//++
     mpz_init(tump1);//++
     mpz_t tump2;//++
     mpz_init(tump2);//++
     mpz_t p;//++
     mpz_init(p);//++
     mpz_t sum;//++
     mpz_init(sum);//++
     
     for(k = 0 ; k < mpz_get_ui(n)-1 ; k++)//++
     {
        /* Dans un premier temps, on cherche l'élément minimum (non */
        /* nul) en valeur absolue dans la colonne k et d'indice i   */
        /* supérieur ou égal à k.                                   */
        
        /*valmin = A[k][k] ;*/ imin = k ;//++
        mpz_set(valmin,A[k][k]); //++
        for(i = k+1 ; i < mpz_get_ui(n) ; i++)//++
        {
           if (mpz_get_ui(valmin) != 0)//++
           {
              mpz_t test1;//++
              mpz_init(test1);//++
              mpz_t test2;//++
              mpz_init(test2);//++
              
              mpz_abs(test1,A[i][k]);//++
              mpz_abs(test2,valmin);//++
              //if(mpz_abs(A[i][k],A[i][k]) < mpz_abs(valmin,valmin))
              if (mpz_cmp(test1,test2) < 0 && mpz_get_ui(A[i][k]) != 0)//++
              {
                 /*valmin = A[i][k] ;*/ //++
                 mpz_set(valmin,A[i][k]);//++
                 imin = i ;
              }
              mpz_clear(test1);
              mpz_clear(test2);
           }
           else 
           {
                 /*valmin = A[i][k] ;*/
                 mpz_set(valmin,A[i][k]);//++
                 imin = i ;
           }     
        }
        
        /* Si l'élément minimum est nul, on peut en déduire */
        /* que la matrice est singulière. Le pogramme est   */
        /* alors interrompu.                                */
        
        if (mpz_get_ui(valmin) == 0) //++
        {
           printf("\n\n\nAttention! Matrice singuliere!\n\n\n") ;
           exit( EXIT_FAILURE ) ;
        }
        
        /* Si la matrice n'est pas singulière, on inverse    */
        /* les éléments de la ligne imax avec les éléments   */
        /* de la ligne k. On fait de même avec le vecteur b. */
        
        for(j = 0 ; j < mpz_get_ui(n) ; j++)//++
        {
           /*tump1 = A[imin][j] ;
           A[imin][j] = A[k][j] ;
           A[k][j] = tump1 ;*/
           mpz_set(tump1,A[imin][j]);//++
           mpz_set(A[imin][j],A[k][j]);//++
           mpz_set(A[k][j],tump1);//++
        }
        
        /*tump2 = b[imin] ;
        b[imin] = b[k] ;
        b[k] = tump2 ;*/
        mpz_set(tump2,b[imin]);//++
        mpz_set(b[imin],b[k]);//++
        mpz_set(b[k],tump2);//++
        
        /* On procède à la réduction de la matrice par la */
        /* méthode d'élimination de Gauss. */
        
        for(i = k+1 ; i < mpz_get_ui(n) ; i++)//++
        {
           //p = A[i][k]/A[k][k] ;
           mpz_cdiv_q(p,A[i][k],A[k][k]); //++
           
           for(j = 0 ; j < mpz_get_ui(n) ; j++)//++
           {
              //A[i][j] = A[i][j] - p*A[k][j] ;
              mpz_submul(A[i][j],p,A[k][j]);//++
           }
           
           //b[i] = b[i] - p*b[k] ; 
           mpz_submul(b[i],p,b[k]);//++
        }
     }   
     
     /* On vérifie que la matrice n'est toujours pas singulière. */
     /* Si c'est le cas, on interrompt le programme. */
     
     if (A[mpz_get_ui(n)-1][mpz_get_ui(n)-1] == 0)//++
     {
        printf("\n\n\nAttention! Matrice singuliere!\n\n\n") ;
        exit( EXIT_FAILURE ) ; 
     }
     
     /* Une fois le système réduit, on obtient une matrice triangulaire */
     /* supérieure et la résolution du système se fait très simplement. */
     
     //x[n-1] = b[n-1]/A[n-1][n-1] ;//++
     
     mpz_cdiv_q(x[mpz_get_ui(n)-1],b[mpz_get_ui(n)-1],A[mpz_get_ui(n)-1][mpz_get_ui(n)-1]);//++
     
     for(i = mpz_get_ui(n)-2 ; i > -1 ; i--)//++
     {
           //sum = 0 ;
           mpz_set_ui(sum,0); //++
           
           for(j = mpz_get_ui(n)-1 ; j > i ; j--)//++
           {
              //sum = sum + A[i][j]*x[j] ;
              mpz_t inter1;//++
              mpz_init(inter1);//++
              mpz_mul(inter1,A[i][j],x[j]);//++
              mpz_add(sum,sum,inter1);//++
              mpz_clear(inter1);//++
           }
           //x[i] = (b[i] - sum)/A[i][i] ;
           mpz_t inter2;//++
           mpz_init(inter2);//++
           mpz_sub(inter2,b[i],sum);//++
           mpz_cdiv_q(x[i],inter2,A[i][i]); //++
           mpz_clear(inter2);//++
           
     }
     
     affiche_solution(x,n);//++
     mpz_clear(valmin);//++
     mpz_clear(tump1);//++
     mpz_clear(tump2);//++
     mpz_clear(p);//++
     mpz_clear(sum);//++
}


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
    //On parcours tous les éléments de la décomposition.\
    Pour chaque prime, ...
    /*while(mpz_cmp_ui(P[i].prime, 1) != 0)
    {


        
        i++;
    }*/
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

    //Si on arrive ici, tous les nombres sont bien dans la FBase...

    //et le prime voulu est bien dans la FBase
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

    /*  Création matrice carrée ******************************* */
    //On alloue la mémoire
    mpz_t** matriceCarree = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));
    int y;
    int h;

    //Là aussi
    for (y = 0 ; y < mpz_get_ui(tailleFBase) ; y++)
    {
        matriceCarree[y] = malloc(mpz_get_ui(tailleFBase) * sizeof(mpz_t));
    }

    //On l'initialise
    for(y=0; y<mpz_get_ui(tailleFBase); y++) //Init
    {
       for (h = 0 ; h < mpz_get_ui(tailleFBase) ; h++)
        {
            mpz_init(matriceCarree[y][h]);
        }
    }

    gmp_printf("Matrice carrée : \n");

    for(y=0; y<mpz_get_ui(tailleFBase); y++) //Init
    {
       for (h = 0 ; h < mpz_get_ui(tailleFBase) ; h++)
        {
            gmp_printf("%Zd ", matriceCarree[y][h]);
        }
        gmp_printf("\n");
    }

    /*  Création tableau résultat ******************************* */
    int x;

    //On alloue la mémoire
    mpz_t* tableauResultat = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));

    //On l'initialise
    for(x=0; x<mpz_get_ui(tailleFBase); x++)
    { 
        mpz_init(tableauResultat[x]);
    }

    gmp_printf("Tableau résultat : \n");
    for(x=0; x<mpz_get_ui(tailleFBase); x++) 
    {
        gmp_printf("%Zd ", tableauResultat[x]);
    }
    printf("\n");

    /*  Création tableau résultat ******************************* */

    //On alloue la mémoire
    mpz_t* solutions = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));

    //On l'initialise
    for(x=0; x<mpz_get_ui(tailleFBase); x++)
    { 
        mpz_init(solutions[x]);
    }

    gmp_printf("Tableau résultat : \n");
    for(x=0; x<mpz_get_ui(tailleFBase); x++) 
    {
        gmp_printf("%Zd ", solutions[x]);
    }
    printf("\n");

    /* 1. Factor base Init ******************************* 
     On veut avoir tous les nombres premiers de 2 à celui correspondant à 
    la taille du FBase, pour comparer par la suite plus simplement. */

    mpz_t* FactorBase = malloc(mpz_get_ui(tailleFBase)*sizeof(mpz_t));

    
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

     mpz_t k2;
             mpz_init(k2);

    //Tant qu'on a pas atteint le nombre d'équations testées voulu, on boucle
    //while(cnt < mpz_get_ui(tailleFBase))
    while(cnt < 40)
    {
        //On calcule generateur^k
        mpz_powm_ui (resPowm, generateur, k, ordre);

        //On le décompose en facteurs premiers
        primeFactDecomp(P, resPowm);

        //Et on donne cette décomposition à une fonction qui va nous dire si ces éléments \
        de cette décomposition sont bien tous dans le factor base. \
        Elle nous dira aussi si le prime voulu est dans cette équation (sinon on la drop).
        test = testInFactorBase(P, FactorBase, tailleFBase, primeVoulu);

        //Si tous les éléments sont bien dans la FactorBase et que notre prime voulu y\
        est aussi, alors :
        if (test == 11)
        {
            //printf("-> Dans la FBase !!\n");
            //printf("-> Sélectionnée !\n");
             
             gmp_printf("\nN°: %d - k : %d Nombre : %Zd\n", cnt, k, resPowm);
             affiche(P);

             int curseurPrime = 0;
             int curseurMCetFB = 0;
             printf("test\n");
             //On va la mettre dans la matrice carrée. On la parcours (ainsi que FBase)
             for (curseurMCetFB = 0; curseurMCetFB < mpz_get_ui(tailleFBase); curseurMCetFB++)
             {
                //gmp_printf("test2 cnt=%d\n", cnt);
                gmp_printf("curseurPrime = %d, curseurMCetFB = %d, FactorBase[] = %Zd, P[] =  %Zd, matriceCarree[] = %Zd\n", \
                    curseurPrime, curseurMCetFB, FactorBase[curseurMCetFB], P[curseurPrime].prime, matriceCarree[cnt][curseurMCetFB]);

                //printf("test3\n");
                //si la colonne correspond bien au prime de la décomposition current
                if(mpz_cmp(FactorBase[curseurMCetFB], P[curseurPrime].prime) == 0)
                {
                    mpz_set(matriceCarree[cnt][curseurMCetFB], P[curseurPrime].pow);  
                    curseurPrime++;
                }
                else
                {
                    mpz_set_ui(matriceCarree[cnt][curseurMCetFB], 0);
                }

             }

             mpz_set_ui(k2, k);
             
             mpz_set(tableauResultat[cnt], k2);


             mpz_nextprime (primeVoulu, primeVoulu);

             //gmp_printf("test %Zd" , FactorBase[mpz_get_ui(tailleFBase)-1]);
             cnt++;

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



gmp_printf("Matrice carrée : \n");

    for(y=0; y<mpz_get_ui(tailleFBase); y++) //Init
    {
       for (h = 0 ; h < mpz_get_ui(tailleFBase) ; h++)
        {
            gmp_printf("%Zd ", matriceCarree[y][h]);
        }
        gmp_printf("\n");
    }

    gmp_printf("\n");

    for(y=0; y<mpz_get_ui(tailleFBase); y++) //Init
    {
       gmp_printf("%Zd ", tableauResultat[y]);
    }


    gauss(matriceCarree, tableauResultat, solutions, tailleFBase);//++






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
