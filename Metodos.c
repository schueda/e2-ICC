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
    for (int j=i+1; j<SL->n; j++) { // Esse loop procura o maior valor na linha
      if (fabs(SL->A[j][i]) > fabs(SL->A[p][i])) {
        p = j;
      }
    }
    if (p != i) {
      trocaLinhas(SL, i, p); // Troca as linhas se forem diferentes
    }

    for (int k=i+1; k<SL->n; k++) {
      double m = SL->A[k][i] / SL->A[i][i]; // Calcula o multiplicador
      SL->A[k][i] = 0; // Zera o valor que está abaixo da diagonal principal
      for (int j=i+1; j<SL->n; j++) {
        SL->A[k][j] = SL->A[k][j] - m * SL->A[i][j]; // Aplica o multiplicador para todas as linhas abaixo da linha i
      }
      SL->b[k] = SL->b[k] - m * SL->b[i]; // Aplica o multiplicador nos itens de b
    }
  }

  resolveSisTri(SL, x); // Resolve o sistema triangular

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
  calcResiduo(SL, x, R); // Popula a matriz de resíduo

  real_t norma = 0; 
  for (int i=0; i<SL->n; i++) { // Soma os quadrados em um acumulador
    norma += R[i] * R[i];
  }
  norma = sqrt(norma); // Tira a raíz quadrada

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

  while ((normaL2Residuo(SL, x) >= erro) && (iteracoes < MAXIT)) { // Enquanto o resíduo for maior que o erro
    for (int i=0; i<SL->n; i++) { 
      real_t soma = 0;
      for (int j=0; j<i; j++) {
        soma += SL->A[i][j] * x[j];
      }
      for (int j=i+1; j<SL->n; j++) {
        soma += SL->A[i][j] * x[j];
      }

      //Coloca todos os valores em um acumulador e depois subtrai do valor de b e divide pelo valor da diagonal

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
  eliminacaoGauss(SL, x, tTotal); // Resolve o sistema pela primeira vez utilizando Eliminação de Gauss

  real_t *res = (real_t *) malloc(SL->n * sizeof(real_t));

  while ((normaL2Residuo(SL, x) > erro) && (iteracoes < MAXIT)) { // Se o resíduo for maior que o erro máximo, continua o loop
    calcResiduo(SL, x, res);
    
    for (int i=0; i<SL->n; i++) { //Troca o vetor B pelo vetor resíduo
      SL->b[i] = res[i];
    }

    resolveSisTri(SL, x); // Resolve o sistema triangular

    iteracoes++;
  }

  *tTotal = timestamp() - t0;

  free(res);
  return iteracoes;
}


