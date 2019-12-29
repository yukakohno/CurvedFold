#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Bezier.h"
#include "util.h"

int Bezier::set( double _x0, double _y0, double _dx0, double _dy0, double len0,
				double _x1, double _y1, double _dx1, double _dy1, double len1,
				char *fname=NULL )
{
	int i;
	double dx0, dy0, dx1, dy1, len;

	x0 = _x0;	y0 = _y0;
	x3 = _x1;	y3 = _y1;

	dx0 = _dx0;	dy0 = _dy0;
	len = sqrt(dx0*dx0 + dy0*dy0);
	x1 = x0 + dx0*len0/len;
	y1 = y0 + dy0*len0/len;

	dx1 = _dx1;	dy1 = _dy1;
	len = sqrt(dx1*dx1 + dy1*dy1);
	x2 = x3 + dx1*len1/len;
	y2 = y3 + dy1*len1/len;

	for( i=0; i<=BTCNT; i++ ){
		double t = (double)i / (double)BTCNT;
		double _t = 1.0-t, t2=t*t, t3=t2*t, _t2=_t*_t, _t3=_t2*_t;
		tt[i] = t;
		tx[i] = _t3*x0 + 3*_t2*t*x1 + 3*_t*t2*x2 + t3*x3;
		ty[i] = _t3*y0 + 3*_t2*t*y1 + 3*_t*t2*y2 + t3*y3;
	}

	if( fname != NULL ){
		FILE *fp = fopen(fname, "w");
		if( fp ){
			fprintf( fp, "i,tt,tx,ty\n" );
			for( i=0; i<=BTCNT; i++ ){
				fprintf( fp, "%d,%f,%f,%f\n", i,tt[i], tx[i], ty[i] );
			}
			fclose( fp );
		}
	}

	return 0;
}

double Bezier::gety( double x )
{
	int i;
	double t;

	// x‚Í’P’²‘‰Á
	for( i=1; i<BTCNT; i++ ){
		if( x < tx[i] ){
			i--;
			break;
		}
	}

	t = tt[i] + (tt[i+1]-tt[i]) * (x-tx[i]) / (tx[i+1]-tx[i]);
	double _t = 1.0-t, t2=t*t, t3=t2*t, _t2=_t*_t, _t3=_t2*_t;
	double xx = _t3*x0 + 3*_t2*t*x1 + 3*_t*t2*x2 + t3*x3;
	double yy = _t3*y0 + 3*_t2*t*y1 + 3*_t*t2*y2 + t3*y3;

	return yy;
}
