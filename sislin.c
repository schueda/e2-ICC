#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"

// Alocaçao de matriz em memória. 
SistLinear_t* alocaSisLin (unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if ( SL ) {
    
    SL->n = n;
    SL->A = (real_t **) malloc(n * sizeof(real_t *));
    SL->b = (real_t *) malloc(n * sizeof(real_t));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SL->A[0] = (real_t *) malloc(n * n * sizeof(real_t));
    if (!(SL->A[0])) {
      liberaSisLin(SL);
      return NULL;
    }
    
    for (int i=1; i < n; ++i) {
      SL->A[i] = SL->A[i-1]+n;
    }
  }
  
  return (SL);
}

// Liberacao de memória
void liberaSisLin (SistLinear_t *SL)
{
  if (SL) {
    if (SL->A) {
      if (SL->A[0]) free (SL->A[0]);
    free(SL->A);
    }
    
    if (SL->b) free(SL->b);

    free(SL);
  }
}


/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: 
     comSolucao, eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void iniSisLin (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max)
{
  unsigned int n = SL->n;
  // para gerar valores no intervalo [0,coef_max]
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // inicializa vetor b
  for (unsigned int i=0; i<n; ++i) {
    SL->b[i] = (real_t)rand() * invRandMax;
  }
    
  if (tipo == hilbert) {
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
	SL->A[i][j] = 1.0 / (real_t)(i+j+1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
	SL->A[i][j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == diagDominante) {
      // aumenta o valor dos termos da diagonal principal
      for (unsigned int i=0; i<n; ++i) {
	real_t soma = 0.0;
	for (unsigned int j=0; j < i; ++j) soma += SL->A[i][j];
	for (unsigned int j=i+1; j < n; ++j) soma += SL->A[i][j];
        SL->A[i][i] += soma;
      }
    }
  }
}


SistLinear_t *lerSisLin ()
{
  unsigned int n;
  SistLinear_t *SL;
  
  scanf("%d",&n);

  SL = alocaSisLin (n);
  
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      scanf ("%lg", &SL->A[i][j]);

  for(int i=0; i < n; ++i)
    scanf ("%lg", &SL->b[i]);
  
  return SL;
}


void prnSisLin (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n  ");
    for(int j=0; j < n; ++j)
      printf ("%10g", SL->A[i][j]);
    printf ("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10g ", v[i]);
  printf ("\n\n");

}

void trocaLinhas(SistLinear_t *SL, int i, int j)
{
  real_t *aux;
  real_t aux2;
  
  aux = SL->A[i];
  SL->A[i] = SL->A[j];
  SL->A[j] = aux;

  aux2 = SL->b[i];
  SL->b[i] = SL->b[j];
  SL->b[j] = aux2;
}

void resolveSisTri(SistLinear_t *SL, real_t *x)
{
  for (int i=SL->n-1; i >= 0; i--) {
    x[i] = SL->b[i];
    for (int j=i+1; j < SL->n; j++) {
      x[i] -= SL->A[i][j] * x[j];
    }
    x[i] /= SL->A[i][i];
  }
}