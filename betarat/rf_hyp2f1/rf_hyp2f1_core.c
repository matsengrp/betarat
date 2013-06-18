#include <gmp.h>
#include <stdio.h>
#include "rf_hyp2f1_core.h"

/* This is a special definition of hyp2f1 which only works when b is a negative integer. In particular this is
 * the _recursively factored_ form mentioned in section 5.5 of Concrete Mathematics by Graham, Knuth and
 * Patashnik. This recursively factored form terminates finitely when b is an integer. However, we do need to
 * use high precision arithmetic here, as there are significant parts of the domain where rounding errors lead
 * to numerically unstable solutions. */
double rf_hyp2f1_core(double a, int b, double c, double w, int mpf_prec)
{
  double d_result;
  mpf_t ma, mb, mc, mw;
  mpf_t result;

  /* i = actual iterator; start_iter - the start iter. We go backwards! */  
  int i;
  int start_iter = -(b+1);

  /* This one will stay the same throughout */
  mpf_init_set_d(mw, w);

  /* Init value for the fold is 1 */
  mp_bitcnt_t prec = mpf_prec;
  mpf_init2(result, prec);
  mpf_set_si(result, 1);
  /*int p = mpf_get_prec(result);*/
  /*printf("Current is %u\n", p);*/

  /* Going to be using and decrementing these in place... */
  mpf_init_set_d(ma, a + start_iter);
  mpf_init_set_si(mb, b + start_iter);
  mpf_init_set_d(mc, c + start_iter);

  /* Let the games begin */ 
  for (i = -(b+1); i >= 0; i--)
  {
    /* Update the result with the mofidications of this iteration */
    mpf_mul(result, result, mw);
    mpf_mul(result, result, ma);
    mpf_mul(result, result, mb);
    mpf_div_ui(result, result, i + 1);
    mpf_div(result, result, mc);
    mpf_add_ui(result, result, 1);

    /* Decrement these counters */
    mpf_sub_ui(ma, ma, 1);
    mpf_sub_ui(mb, mb, 1);
    mpf_sub_ui(mc, mc, 1);
  }
  /* extract as double */
  d_result = mpf_get_d(result);
  
  /* Cleaning up shop */
  mpf_clear(result);
  mpf_clear(ma);
  mpf_clear(mb);
  mpf_clear(mc);
  mpf_clear(mw);

  return d_result;
}

