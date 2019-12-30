#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "papermodel.h"
#include "spline.h"
#include "util.h"

int papermodel::setSpherePieceBoundary( int divnum, double crv, curve_type ctype )
{
	int ret=0;

	dccnt = 0;
	crcnt = lcrcnt = rcrcnt = 1;

	if( divnum>4 ){
		double angle, cosa, sina;

		crcnt = rcrcnt = 2;
		crease *c0 = rcrs[0] = &(crs[0]);
		crease *c1 = rcrs[1] = &(crs[1]);
		c1->Pcnt = c0->Pcnt;
		c1->Xcnt = c0->Xcnt;

		// lrft
		angle = M_PI * ( 0.25 + 1.0/(double)divnum );
		cosa = cos(angle);
		sina = sin(angle);
		c1->m2[0] = cosa;	c1->m2[1] = sina;	c1->m2[2] = 0.0;
		c1->m2[3] = -sina;	c1->m2[4] = cosa;	c1->m2[5] = 0.0;
		c1->m2[6] = 0.0;	c1->m2[7] = 0.0;	c1->m2[8] = 1.0;
		//c1->m2[6] = -1.0;	c1->m2[7] = 1.0;	c1->m2[8] = 1.0; // è≠ÇµéÜÇ©ÇÁÇÕÇ›èoÇ∑

		for( int i=0; i<c1->Pcnt; i++ ){
			c1->Px2d[i] = crv;
		}
		crease::interpolate_spline( c1->Px2d, c1->Pcnt, c1->k2d, c1->Xcnt );
		c1->flg_xspc_def=1;
		c1->calcXTN2d( c1->m2 ); // k2d -> X2d -> calcTN2d();
		c1->flg_xspc_def=0;
#if 0
		memcpy( dcurve[dccnt].cvx, c1->Xx2d, sizeof(double)*c1->Xcnt ); // MIN( Xcnt, MAXDCX=1000 )
		memcpy( dcurve[dccnt].cvy, c1->Xy2d, sizeof(double)*c1->Xcnt );
		dcurve[dccnt].cvx[c1->Xcnt] = c1->Xx2d[c1->Xcnt-2] + 15.0*(c1->Xx2d[c1->Xcnt-1]-c1->Xx2d[c1->Xcnt-2]);
		dcurve[dccnt].cvy[c1->Xcnt] = c1->Xy2d[c1->Xcnt-2] + 15.0*(c1->Xy2d[c1->Xcnt-1]-c1->Xy2d[c1->Xcnt-2]);
#else
		dcurve[dccnt].cvx[0] = c1->Xx2d[c1->Xcnt-2] + 10.0*(c1->Xx2d[c1->Xcnt-1]-c1->Xx2d[c1->Xcnt-2]);
		dcurve[dccnt].cvy[0] = c1->Xy2d[c1->Xcnt-2] + 10.0*(c1->Xy2d[c1->Xcnt-1]-c1->Xy2d[c1->Xcnt-2]);	
		for( int i=1; i<=c1->Xcnt; i++ ){
			dcurve[dccnt].cvx[i] = c1->Xx2d[ c1->Xcnt - i ];
			dcurve[dccnt].cvy[i] = c1->Xy2d[ c1->Xcnt - i ];
		}
#endif
		dcurve[dccnt].ctype = ctype;
		dcurve[dccnt].cvcnt = c1->Xcnt+1;
		dccnt++;

		// right
		angle = -2.0*M_PI / (double)divnum;
		cosa = cos(angle);
		sina = sin(angle);
		for( int i=0; i<c1->Xcnt; i++ ){
			double tmpx = cosa*c1->Xx2d[i] - sina*c1->Xy2d[i];
			double tmpy = sina*c1->Xx2d[i] + cosa*c1->Xy2d[i];
			c1->Xx2d[i] = tmpx;
			c1->Xy2d[i] = tmpy;
//			c1->Xx2d[i] = tmpx + 2.0; // è≠ÇµéÜÇ©ÇÁÇÕÇ›èoÇ∑
//			c1->Xy2d[i] = tmpy - 2.0;
		}
		memcpy( dcurve[dccnt].cvx, c1->Xx2d, sizeof(double)*c1->Xcnt ); // MIN( Xcnt, MAXDCX=1000 )
		memcpy( dcurve[dccnt].cvy, c1->Xy2d, sizeof(double)*c1->Xcnt );
		dcurve[dccnt].cvx[c1->Xcnt] = c1->Xx2d[c1->Xcnt-2] + 10.0*(c1->Xx2d[c1->Xcnt-1]-c1->Xx2d[c1->Xcnt-2]);
		dcurve[dccnt].cvy[c1->Xcnt] = c1->Xy2d[c1->Xcnt-2] + 10.0*(c1->Xy2d[c1->Xcnt-1]-c1->Xy2d[c1->Xcnt-2]);
		dcurve[dccnt].ctype = ctype;
		dcurve[dccnt].cvcnt = c1->Xcnt+1;
		dccnt++;

		// reset
		crcnt = rcrcnt = lcrcnt = 1;
		crs[1].init();
		crs[2].init();
#if 1
		if( ctype==CTYPE_FOLD ){
			// left
			int dcidx = 0;
			c1 = lcrs[lcrcnt] = &(crs[crcnt]);	lcrcnt++;	crcnt++;
			c1->rl = -1;
			c1->org_idx = dcidx;
			c1->org_cnt = dcurve[dcidx].cvcnt;
			c1->org_x = dcurve[dcidx].cvx;
			c1->org_y = dcurve[dcidx].cvy;
			c1->org2CP();
			c1->flg_org = 1;

			// right
			dcidx = 1;
			c1 = rcrs[rcrcnt] = &(crs[crcnt]);	rcrcnt++;	crcnt++;
			c1->rl = 1;
			c1->org_idx = dcidx;
			c1->org_cnt = dcurve[dcidx].cvcnt;
			c1->org_x = dcurve[dcidx].cvx;
			c1->org_y = dcurve[dcidx].cvy;
			c1->org2CP();
			c1->flg_org = 1;

			//dccnt = 0;
		}
#endif
	}

	return ret;
}
