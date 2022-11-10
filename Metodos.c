#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
  for (int i=0; i<SL->n; i++) {

    // Pivoteamento parcial
    int p = i;
    for (int j=i+1; j<SL->n; j++) {
      if (fabs(SL->A[j][i]) > fabs(SL->A[p][i])) {
        p = j;
      }
    }
    if (p != i) {
      trocaLinhas(SL, i, p);
    }

    for (int k=i+1; k<SL->n; k++) {
      double m = SL->A[k][i] / SL->A[i][i];
      SL->A[k][i] = 0;
      for (int j=i+1; j<SL->n; j++) {
        SL->A[k][j] = SL->A[k][j] - m * SL->A[i][j];
      }
      SL->b[k] = SL->b[k] - m * SL->b[i];
    }
  }

  resolveSisTri(SL, x);
}


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  int *R = (int *) malloc(SL->n * sizeof(int));

  for (int i=0; i<SL->n; i++) {
    R[i] = 0;
    for (int j=0; j<SL->n; j++) {
      R[i] += SL->A[i][j] * x[j];
    }
    R[i] -= SL->b[i];
  }

  real_t norma = 0;
  for (int i=0; i<SL->n; i++) {
    norma += R[i] * R[i];
  }
  norma = sqrt(norma);

  free(R);

  return norma;
}


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  int iteracoes = 0;
  while (normaL2Residuo(SL, x) > erro) {
    for (int i=0; i<SL->n; i++) {
      x[i] = SL->b[i];
      for (int j=0; j<i; j++) {
        x[i] -= SL->A[i][j] * x[j];
      }
      for (int j=i+1; j<SL->n; j++) {
        x[i] -= SL->A[i][j] * x[j];
      }
      x[i] /= SL->A[i][i];
    }
    iteracoes++;
  }
  return iteracoes;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{

}


