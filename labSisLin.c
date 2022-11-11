#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(time(NULL));
  
  SistLinear_t *SL;
  unsigned int tam[9] = { 10, 30, 50, 128, 256, 512, 1000, 2000, 3000 };

  for (int i=0; i<9; i++) {
    SL = alocaSisLin(tam[i]);
    if (!SL) {
      perror("Erro na alocação do sistema linear\n");
    } else {
      real_t *x = (real_t *) calloc(tam[i], tam[i] * sizeof(real_t));
      iniSisLin(SL, generico, COEF_MAX);

      real_t *res = (real_t *) calloc(tam[i], tam[i] * sizeof(real_t));
      double tempo;

      printf("%f\n", normaL2Residuo(SL, x));

      printf("iterações: %d\n", gaussSeidel(SL, x, ERRO, &tempo));

      printf("%f\n", normaL2Residuo(SL, x));
      
      
      liberaSisLin(SL);
      free(x);
      free(res);
    }
  }


  return 0;
}

