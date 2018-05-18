#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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


main(){
	
	
	return 1;	
}
