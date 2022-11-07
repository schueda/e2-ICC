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

  for (int i=0; i<1; i++) {
    SL = alocaSisLin(tam[i]);
    if (!SL) {
      perror("Erro na alocação do sistema linear\n");
    } else {
      real_t *x = (real_t *) malloc(tam[i] * sizeof(real_t));
      iniSisLin(SL, diagDominante, 10);
      eliminacaoGauss(SL, x, NULL);
      liberaSisLin(SL);
    }
  }


  return 0;
}

