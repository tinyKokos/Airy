#include <string.h>
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "flint/profiler.h"

// int main(int argc, char *argv[])
int main(int argc, char *argv[])
{
    slong prec, i;
    double dxa, dxb, dya, dyb, xd;
    arb_t xa, xb, ya, yb, w, x, y, z;
    acb_t a, b, c, d, e;

    double scale = 1.0;
    double multiplier = 3.0 * sqrt(3) / 8.0 / scale;

    dxa = 0.01;
    dxb = 2.0 + dxa;

    arb_init(w);
    arb_init(x);
    arb_init(y);
    arb_init(z);

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    acb_init(e);

prec = 30 * 3.33;
    for (xd = dxa; xd <= dxb; xd += dxa)
    {
//	printf("x = %16.15f\n", xd);
    	arb_set_d(x, xd);
//	flint_printf("\n x = ");
	arb_printd(x, 50);
	flint_printf(" ");

// do all of the math in arb (after copying in the value of x)
  	arb_set_d(z, 1.0 / cbrt(3.0 * multiplier * xd));
//	arb_mul(y, 9, x, prec);
//	arb_sqrt(w, 3, prec);
//	arb_mul(z, w, y, prec);
//	arb_div(x, z, 8, prec);
	arb_hypgeom_airy(w, NULL, NULL, NULL, z, prec);
//	flint_printf("\n arb airy only = ");
//	arb_printn(w, 50, 0);
	flint_printf(" ");
	arb_printd(w, 50);
	flint_printf(" ");
//	arb_print(w);
//	flint_printf(" ");

/*
    	acb_set_d(a, 1.0 / cbrt(3.0 * multiplier * xd));
	acb_hypgeom_airy(b, NULL, NULL, NULL, a, prec);
	flint_printf("\n acb airy only = ");
	acb_printn(b, 50, 0);
*/

    	arb_set_d(y, 1.0 / cbrt(3.0 * multiplier * xd) / xd);
	arb_mul(z, w, y, prec);
//	flint_printf("\n z = ");
	arb_printd(z, 50);
	flint_printf("\n");
    }
	flint_printf("\n");

    arb_clear(w);
    arb_clear(x);
    arb_clear(y);
    arb_clear(z);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    acb_clear(e);

    flint_cleanup();
    return 0;
}
