#include <string.h>
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "flint/profiler.h"

int main()
{
	slong prec;
	double dxa, dxb, xd;
	arb_t xa, xb, res, nu, x, z;

	double scale = 1.0;
	double multiplier = 3.0 * sqrt(3) / 8.0 / scale;

	dxa = 0.01;
	dxb = 2.0 + dxa;

	arb_init(res);
     	arb_init(x);
    	arb_init(nu);
	arb_init(z);

	prec = 30 * 3.33;
	// adjust this command to adjust nu in Bessel code
	arb_set_d(nu, 1/3);

	// I think that it might be possible to remove the for loop if you can set
	//the arb_t ball of x  to all the points that you need.

	for(xd = dxa; xd <= dxb; xd += dxa)
	{
		arb_set_d(x, xd);
		arb_printd(x, 50);
		flint_printf(" ");
		arb_set_d(z, 1.0 / cbrt(3.0 * multiplier * xd));

		arb_hypgeom_bessel_j(res, nu, z, prec);
		flint_printf(" ");
		arb_printd(res, 50);
		flint_printf("\n");
	}

	flint_printf("\n");

	arb_clear(res);
       	arb_clear(x);
    	arb_clear(nu);
	arb_clear(z);
	
	flint_cleanup();
	return 0;
}
