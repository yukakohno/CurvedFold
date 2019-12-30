#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#define ZERO 0.000001	

#include "papermodel.h"
#include "Bspline.h"
#include "util.h"

int papermodel::setRulLen_curve( double *cvx, double *cvy, double *rx, double *ry, int cvcnt, 
								double *cvx1, double *cvy1, int cvcnt1, double *_Cx, double *_Cy,
								double *rlen, ruling_end_type *rflg )
{
	int ret=0, Scnt = 10000, nv2[] = {0,0,0,0, 1, 2, 3,3,3,3};
	double *Sx=NULL, *Sy=NULL, *L=NULL, Cx[CCNT], Cy[CCNT], t0[MAX_SPCNT], *t1=NULL, maxlen=pw+ph;
	Bspline bs;
	t1 = ht1;
	Sx = hSx;
	Sy = hSy;
	L = hL;
	ret = bs.calcT( cvcnt1, 3, 1, t0, cvx1, cvy1, NULL );
	ret = bs.calcT( Scnt, 3, 0, t1);
	ret = bs.calcCP( cvx1, cvy1, NULL, cvcnt1, nv2, 10, 4/*degree*/, t0, CCNT, Cx, Cy, NULL );
	ret = bs.calcBSpline( Cx, Cy, NULL, CCNT, nv2, 10, 4/*degree*/, Scnt, t1, Sx, Sy, NULL );
	// Sx, Sy Ç∆ruling ÇÃåç∑ÇîªíËÅAruling í∑Ç≥(rllen, rrlen)ÇçXêV
	for( int i=0; i<cvcnt+1; i++ ){
		double x0, y0, x1, y1, ix,iy,l0,l1;
		x0 = cvx[i];
		y0 = cvy[i];
		x1 = cvx[i]+rx[i]*maxlen;
		y1 = cvy[i]+ry[i]*maxlen;
		rlen[i] = maxlen;
		for( int j=0; j<Scnt-1; j++ ){
			int ret = intersectionOfLine( x0, y0, x1, y1, Sx[j], Sy[j], Sx[j+1], Sy[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			if( rlen[i] > l0 ){
				rflg[i] = RETYPE_CURVE;
				rlen[i] = l0;
			}
		}
	}
	if( _Cx ){ memcpy( _Cx, Cx, sizeof(double)*CCNT ); }
	if( _Cy ){ memcpy( _Cy, Cy, sizeof(double)*CCNT ); }
	t1=NULL;
	Sx=NULL;
	Sy=NULL;
	L=NULL;
	return ret;
}

int papermodel::setRulLen_SplineCP( double *cvx, double *cvy, double *rx, double *ry, int cvcnt, 
								   double *Cx, double *Cy,
								   double *rlen, ruling_end_type *rflg )
{
	int ret=0, Scnt = 10000, nv2[] = {0,0,0,0, 1, 2, 3,3,3,3};
	double *Sx=NULL, *Sy=NULL, *t1=NULL, maxlen=pw+ph;
	Bspline bs;
	t1 = ht1;
	Sx = hSx;
	Sy = hSy;
	//ret = bs.calcT( Scnt-1, 3, 0, t1); t1[Scnt-1]=3.0;
	ret = bs.calcBSpline( Cx, Cy, NULL, CCNT, nv2, 10, 4/*degree*/, Scnt, t1, Sx, Sy, NULL );
	for( int i=0; i<cvcnt; i++ ){
		double x0, y0, x1, y1, ix,iy,l0,l1;
		x0 = cvx[i];
		y0 = cvy[i];
		x1 = cvx[i]+rx[i]*maxlen;
		y1 = cvy[i]+ry[i]*maxlen;
		rlen[i] = maxlen;
		for( int j=0; j<Scnt-1; j++ ){
			int ret = intersectionOfLine( x0, y0, x1, y1, Sx[j], Sy[j], Sx[j+1], Sy[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			if( rlen[i] > l0 ){
				rflg[i] = RETYPE_CURVE;
				rlen[i] = l0;
			}
		}
	}
	t1=NULL;
	Sx=NULL;
	Sy=NULL;
	return ret;
}

int papermodel::setRulLen_SplineCP( double *cvx, double *cvy, double *rx, double *ry, int si, int ei, 
								   double *Cx, double *Cy,
								   double *rlen, ruling_end_type *rflg )
{
	int ret = setRulLen_SplineCP( &(cvx[si]), &(cvy[si]), &(rx[si]), &(ry[si]), ei-si+1, Cx, Cy, &(rlen[si]), &(rflg[si]) );
	return ret;
}
