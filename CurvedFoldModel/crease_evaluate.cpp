#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "crease.h"
#include "util.h"

int crease::check180( double *errdata, int row, int col )
{
	int ret=0;

	for( int i=Xsidx; i<=Xeidx; i++ ){
		if( rllen[i]==0 || rrlen[i]==0.0 ){
			//fprintf( fp, "%d,0,0,0,0,0,0,0,0,0\n", i);
			memset( &(errdata[i*col]), 0, sizeof(double)*col );
			continue;
		}
		double ip0r, ip1r, ip0l, ip1l, a0r, a1r, a0l, a1l;
		ip0r = -Tx[i-1]*rrx[i] -Ty[i-1]*rry[i] -Tz[i-1]*rrz[i];
		ip1r =  Tx[i]  *rrx[i] +Ty[i]  *rry[i] +Tz[i]  *rrz[i];
		ip0l = -Tx[i-1]*rlx[i] -Ty[i-1]*rly[i] -Tz[i-1]*rlz[i];
		ip1l =  Tx[i]  *rlx[i] +Ty[i]  *rly[i] +Tz[i]  *rlz[i];
		a0r = acos( ip0r );
		a1r = acos( ip1r );
		a0l = acos( ip0l );
		a1l = acos( ip1l );
		//fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	i, ip0r, ip1r, ip0l, ip1l, a0r, a1r, a0l, a1l, a0r+a1r+a0l+a1l-2*M_PI);
		errdata[i*col+0] = ip0r;
		errdata[i*col+1] = ip1r;
		errdata[i*col+2] = ip0l;
		errdata[i*col+3] = ip1l;
		errdata[i*col+4] = a0r;
		errdata[i*col+5] = a1r;
		errdata[i*col+6] = a0l;
		errdata[i*col+7] = a1l;
		errdata[i*col+8] = a0r+a1r+a0l+a1l-2*M_PI;
		if( i-Xsidx > row ){
			printf( "row: overflow\n" );
		}
	}
end:
	return ret;
}

int crease::checkquadplane( double *errdata, int row, int col )
{
	int ret=0;

	for( int i=Xsidx; i<Xeidx; i++ ){
		double dx0,dy0,dz0, dx1,dy1,dz1, dxL,dyL,dzL, opxL,opyL,opzL, ipL, dxR,dyR,dzR, opxR,opyR,opzR, ipR;
		dx0 = dx1 = Xx[i+1] - Xx[i];
		dy0 = dy1 = Xy[i+1] - Xy[i];
		dz0 = dz1 = Xz[i+1] - Xz[i];
		normalize_v3( &dx0, &dy0, &dz0 );
		if( rllen[i]==0.0 || rllen[i+1]==0.0 ){
			opxL = opyL = opzL = dxL = dyL = dzL = ipL = 0.0;
		} else {
			opxL = rly[i]*dz0 - rlz[i]*dy0;
			opyL = rlz[i]*dx0 - rlx[i]*dz0;
			opzL = rlx[i]*dy0 - rly[i]*dx0;
			normalize_v3( &opxL, &opyL, &opzL ); // 面の法線方向
			dxL = rlx[i+1]*rllen[i+1] + dx1; // X[i]を原点としたあと1点の座標
			dyL = rly[i+1]*rllen[i+1] + dy1;
			dzL = rlz[i+1]*rllen[i+1] + dz1;
			ipL = opxL*dxL + opyL*dyL + opzL*dzL;
		}
		if( rrlen[i]==0.0 || rrlen[i+1]==0.0 ){
			opxR = opyR = opzR = dxR = dyR = dzR = ipR = 0.0;
		} else {
			opxR = rry[i]*dz0 - rrz[i]*dy0;
			opyR = rrz[i]*dx0 - rrx[i]*dz0;
			opzR = rrx[i]*dy0 - rry[i]*dx0;
			normalize_v3( &opxR, &opyR, &opzR ); // 面の法線方向
			dxR = rrx[i+1]*rrlen[i+1] + dx1;
			dyR = rry[i+1]*rrlen[i+1] + dy1;
			dzR = rrz[i+1]*rrlen[i+1] + dz1;
			ipR = opxR*dxR + opyR*dyR + opzR*dzR;
		}
		//fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
		//	i, dx0,dy0,dz0, opxL,opyL,opzL, dxL,dyL,dzL, ipL, opxR,opyR,opzR, dxR,dyR,dzR, ipR );
		errdata[i*col+0] = dx0;
		errdata[i*col+1] = dy0;
		errdata[i*col+2] = dz0;
		errdata[i*col+3] = opxL;
		errdata[i*col+4] = opyL;
		errdata[i*col+5] = opzL;
		errdata[i*col+6] = dxL;
		errdata[i*col+7] = dyL;
		errdata[i*col+8] = dzL;
		errdata[i*col+9] = ipL;
		errdata[i*col+10] = opxR;
		errdata[i*col+11] = opyR;
		errdata[i*col+12] = opzR;
		errdata[i*col+13] = dxR;
		errdata[i*col+14] = dyR;
		errdata[i*col+15] = dzR;
		errdata[i*col+16] = ipR;
	}
	return ret;
}

int crease::checkRulingCross( double *errdata, int row, int col )
{
	int ret=0, li=0;
	double lex[MAX_SPCNT], ley[MAX_SPCNT], rex[MAX_SPCNT], rey[MAX_SPCNT];
	memset( errdata, 0, sizeof(double)*row*col );

	for( int i=Xsidx;  i<=Xeidx; i++ ){
		lex[i] = Xx2d[i]+rllen[i]*rlx_cp[i];
		ley[i] = Xy2d[i]+rllen[i]*rly_cp[i];
		rex[i] = Xx2d[i]+rrlen[i]*rrx_cp[i];
		rey[i] = Xy2d[i]+rrlen[i]*rry_cp[i];
	}

	for( int li=Xsidx, i=Xsidx+1; i<=Xeidx; i++ ){
		double ix,iy, l0,l1, l0_,l1_, a,b,c,d, area;
		if( rllen[i]==0.0 ){
			continue;
		}
		int ret = intersectionOfLine( Xx2d[li],Xy2d[li],lex[li],ley[li], Xx2d[i],Xy2d[i],lex[i],ley[i], &ix, &iy, &l0, &l1 );
		if( ret<0 ){
			li=i;
			continue;
		}
		// サラスの公式 O(0,0),A(a,b),B(c,d) -> area of OAB = abs(a*d-b*c)*0.5
		l0_ = rllen[li]-l0;
		l1_ = rllen[i]-l1;
		a = l0_*rlx_cp[li];
		b = l0_*rly_cp[li];
		c = l1_*rlx_cp[i];
		d = l1_*rly_cp[i];
		errdata[li*col ] = fabs(a*d-b*c)*0.5;
		li=i;
	}

	for( int li=Xsidx, i=Xsidx+1; i<=Xeidx; i++ ){
		double ix,iy, l0,l1, l0_,l1_, a,b,c,d, area;
		if( rrlen[i]==0.0 ){
			continue;
		}
		int ret = intersectionOfLine( Xx2d[li],Xy2d[li],rex[li],rey[li], Xx2d[i],Xy2d[i],rex[i],rey[i], &ix, &iy, &l0, &l1 );
		if( ret<0 ){
			li=i;
			continue;
		}
		// サラスの公式 O(0,0),A(a,b),B(c,d) -> area of OAB = abs(a*d-b*c)*0.5
		l0_ = rrlen[li]-l0;
		l1_ = rrlen[i]-l1;
		a = l0_*rrx_cp[li];
		b = l0_*rry_cp[li];
		c = l1_*rrx_cp[i];
		d = l1_*rry_cp[i];
		errdata[li*col+1] = fabs(a*d-b*c)*0.5;
		li=i;
	}

end:
	return ret;
}
