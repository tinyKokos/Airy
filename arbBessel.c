#include <string.h>
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "flint/profiler.h"

int main()
{
	slong prec, i;
	double dxa, dxb, dya, dyb, xd;
	arb_t xa, xb, ya, yb, w, x, y, z;
	acb_t res, nu, z;

	double scale = 1.0;
	double multiplier = 3.0 * sqrt(3) / 8.0 / scale;

	dxa = 0.01;
	dxb = 2.0 + dxa;

	arb_init(w);
     	arb_init(x);
    	arb_init(y);
    	arb_init(z);

    	arb_init(res);    
	arb_init(nu);    
	arb_init(z); 

	prec = 30 * 3.33;

	for(xd = dxa; xd <= dxb; xd += dxa)
	{
		arb_set_d(x, xd);
		arb_printd(x, 50);
		flint_printf(" ");
		arb_set_d(z, 1.0 / cbrt(3.0 * multiplier * xd));

		arb_hypgeom_bessel_j(res, nu, z, prec);
		flint_printf(" ");
		arb_printd(res, 50);
		flint_printf(" ");

		arb_set_d(y, 1.0 / cbrt(3.0 * multiplier * xd) / xd);
		arb_mul(z, res, y, prec);
		arb_printd(z, 50);
		flint_printf("\n");
	}

	flint_printf("\n");

	arb_clear(w);
       	arb_clear(x);
    	arb_clear(y);
    	arb_clear(z);

    	arb_clear(res);    
	arb_clear(nu);    
	arb_clear(z);	
	
	flint_cleanup();
	return 0;
}
