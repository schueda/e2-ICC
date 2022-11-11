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
  double t0 = timestamp();

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

  *tTotal = timestamp() - t0;

  return 0;
}


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  real_t *R = (real_t *) malloc(SL->n * sizeof(real_t));
  calcResiduo(SL, x, R);

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
  double t0 = timestamp();

  while ((normaL2Residuo(SL, x) >= erro) && (iteracoes < MAXIT)) {
    for (int i=0; i<SL->n; i++) {
      real_t soma = 0;
      for (int j=0; j<i; j++) {
        soma += SL->A[i][j] * x[j];
      }
      for (int j=i+1; j<SL->n; j++) {
        soma += SL->A[i][j] * x[j];
      }

      x[i] = (SL->b[i] - soma)/SL->A[i][i];
    }
    iteracoes++;
  }

  *tTotal = timestamp() - t0;

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
  double t0 = timestamp();
  int iteracoes = 0;
  eliminacaoGauss(SL, x, tTotal);

  real_t *res = (real_t *) malloc(SL->n * sizeof(real_t));

  while ((normaL2Residuo(SL, x) > erro) && (iteracoes < MAXIT)) {
    calcResiduo(SL, x, res);
    
    for (int i=0; i<SL->n; i++) {
      SL->b[i] = res[i];
    }

    resolveSisTri(SL, x);

    iteracoes++;
  }

  *tTotal = timestamp() - t0;

  free(res);
  return iteracoes;
}


