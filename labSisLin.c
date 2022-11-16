#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(202202);
  
  SistLinear_t *SL;
  unsigned int tam[9] = { 10, 30, 50, 128, 256, 512, 1000, 2000, 3000 };

  printf("Diagonal dominante\n");

  printf("    n|            t_egp| normaResiduo_egp|             t_gs|  it_gs|  normaResiduo_gs|            t_ref| it_ref| normaResiduo_ref\n");
  for (int i=0; i<9; i++) {
    SL = alocaSisLin(tam[i]);
    if (!SL) {
      perror("Erro na alocação do sistema linear\n");
    } else {
      real_t *x = (real_t *) calloc(tam[i], tam[i] * sizeof(real_t));

      double tempo;
      int iteracoes;

      printf("%5d", tam[i]);

      iniSisLin(SL, diagDominante, COEF_MAX);
      eliminacaoGauss(SL, x, &tempo);
      printf("|%17e", tempo);
      printf("|%17e", normaL2Residuo(SL, x));
      
      iniSisLin(SL, diagDominante, COEF_MAX);
      iteracoes = gaussSeidel(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e", normaL2Residuo(SL, x));


      iniSisLin(SL, diagDominante, COEF_MAX);
      iteracoes = refinamento(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e\n", normaL2Residuo(SL, x));
      
      liberaSisLin(SL);
      free(x);
    }
  }

  printf("\n");
  printf("Genérico\n");

  printf("    n|            t_egp| normaResiduo_egp|             t_gs|  it_gs|  normaResiduo_gs|            t_ref| it_ref| normaResiduo_ref\n");
  for (int i=0; i<9; i++) {
    SL = alocaSisLin(tam[i]);
    if (!SL) {
      perror("Erro na alocação do sistema linear\n");
    } else {
      real_t *x = (real_t *) calloc(tam[i], tam[i] * sizeof(real_t));

      double tempo;
      int iteracoes;

      printf("%5d", tam[i]);

      iniSisLin(SL, generico, COEF_MAX);
      eliminacaoGauss(SL, x, &tempo);
      printf("|%17e", tempo);
      printf("|%17e", normaL2Residuo(SL, x));
      
      iniSisLin(SL, generico, COEF_MAX);
      iteracoes = gaussSeidel(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e", normaL2Residuo(SL, x));


      iniSisLin(SL, generico, COEF_MAX);
      iteracoes = refinamento(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e\n", normaL2Residuo(SL, x));
      
      liberaSisLin(SL);
      free(x);
    }
  }

  printf("\n");
  printf("Hilbert\n");

  printf("    n|            t_egp| normaResiduo_egp|             t_gs|  it_gs|  normaResiduo_gs|            t_ref| it_ref| normaResiduo_ref\n");
  for (int i=0; i<9; i++) {
    SL = alocaSisLin(tam[i]);
    if (!SL) {
      perror("Erro na alocação do sistema linear\n");
    } else {
      real_t *x = (real_t *) calloc(tam[i], tam[i] * sizeof(real_t));

      double tempo;
      int iteracoes;

      printf("%5d", tam[i]);

      iniSisLin(SL, hilbert, COEF_MAX);
      eliminacaoGauss(SL, x, &tempo);
      printf("|%17e", tempo);
      printf("|%17e", normaL2Residuo(SL, x));
      
      iniSisLin(SL, hilbert, COEF_MAX);
      iteracoes = gaussSeidel(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e", normaL2Residuo(SL, x));


      iniSisLin(SL, hilbert, COEF_MAX);
      iteracoes = refinamento(SL, x, ERRO, &tempo);
      printf("|%17e", tempo);
      printf("|%7d", iteracoes);
      printf("|%17e\n", normaL2Residuo(SL, x));

      liberaSisLin(SL);
      free(x);
    }
  }


  return 0;
}

