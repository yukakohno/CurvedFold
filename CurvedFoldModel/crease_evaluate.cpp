#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

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

double crease::calcRulingCross(double* xx, double* xy, double* rx, double* ry, double* rlen, int cvcnt,
	double* errdata=NULL, int row=0, int col=0)
{
	double ex[MAX_SPCNT], ey[MAX_SPCNT], diff = 0.0;

	for (int i = 0; i < cvcnt; i++) {
		ex[i] = xx[i] + rlen[i] * rx[i];
		ey[i] = xy[i] + rlen[i] * ry[i];
	}
#if 0
	{
		FILE* fp = fopen("output/calcRulCross.csv", "w");
		if (fp) {
			for (int i = 4; i < cvcnt - 4; i++) {
				fprintf(fp, "%f,%f\n%f,%f\n\n", xx[i], xy[i], ex[i], ey[i]);
			}
			fclose(fp);
		}
	}
#endif
	for (int li = 0, i = 1; i < cvcnt; i++) {
		double ix, iy, l0, l1, l0_, l1_, a, b, c, d, area;
		if (rlen[i] == 0.0) {
			continue;
		}
		int ret = intersectionOfLine(xx[li], xy[li], ex[li], ey[li], xx[i], xy[i], ex[i], ey[i], &ix, &iy, &l0, &l1);
		if (ret < 0) {
			li = i;
			continue;
		}
		// サラスの公式 O(0,0),A(a,b),B(c,d) -> area of OAB = abs(a*d-b*c)*0.5
		l0_ = rlen[li] - l0;
		l1_ = rlen[i] - l1;
		a = l0_ * rx[li];
		b = l0_ * ry[li];
		c = l1_ * rx[i];
		d = l1_ * ry[i];
		area = sqrt(fabs(a * d - b * c)) + 0.01;
		diff += area;
		if (errdata) {
			errdata[li * col] = fabs(a * d - b * c) * 0.5;
		}
		li = i;
	}
end:
	//return diff*0.1 / (double)cvcnt;
	return diff / (double)cvcnt;
}

int crease::checkRulingCross()
{
	double area = 0.0;
	return crease::checkRulingCross(area);
}

int crease::checkRulingCross(double& area)
{
	double area_l = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), &(rllen[Xsidx]), Xeidx - Xsidx + 1);
	double area_r = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), &(rrlen[Xsidx]), Xeidx - Xsidx + 1);
	area = area_l + area_r;
	if (area > 0.0) {
		return 1;
	}
	else {
		return 0;
	}
}

int crease::checkRulingCross(int rl, double& area)
{
	if (rl < 0) {
		area = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), &(rllen[Xsidx]), Xeidx - Xsidx + 1);
	}
	else {
		area = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), &(rrlen[Xsidx]), Xeidx - Xsidx + 1);
	}
	if (area > 0.0) {
		return 1;
	}
	else {
		return 0;
	}
}

// return index of rulings crossing, -1 if no crossing
int crease::checkRulingCross( int rl ) // -1:left, 1:right
{
	double area = 0.0;
	if (rl < 0) {
		area = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), &(rllen[Xsidx]), Xeidx - Xsidx + 1);
	}
	else {
		area = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), &(rrlen[Xsidx]), Xeidx - Xsidx + 1);
	}
	if (area > 0.0) {
		return 1;
	}
	else {
		return 0;
	}
}

int crease::checkRulingCross( double *errdata, int row, int col )
{
	double area_l = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), &(rllen[Xsidx]), Xeidx - Xsidx + 1,
		errdata, row, col);
	double area_r = calcRulingCross(&(Xx2d[Xsidx]), &(Xy2d[Xsidx]), &(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), &(rrlen[Xsidx]), Xeidx - Xsidx + 1,
		&(errdata[1]), row, col);
	if (area_l + area_r > 0.0) {
		return 1;
	}
	else {
		return 0;
	}
}

int crease::checkRulingAngle()
{
	int ret = 0;

	for (int i = Xsidx; i <= Xeidx; i++) {
		double dx0 = 0, dy0 = 0, dx1 = 0, dy1 = 0;
		if (Xsidx < i) {
			dx0 = Xx2d[i] - Xx2d[i - 1];
			dy0 = Xy2d[i] - Xy2d[i - 1];
		}
		if (i < Xeidx) {
			dx1 = Xx2d[i + 1] - Xx2d[i];
			dy1 = Xy2d[i + 1] - Xy2d[i];
		}
		double ip = rlx_cp[i] * Tx2d[i] + rly_cp[i] * Ty2d[i];
		double rl = 0;
		if (ip > 0.0 && i < Xeidx) {
			rl = rlx_cp[i] * dy1 - rly_cp[i] * dx1;
		}
		else if (ip < 0.0 && Xsidx < i) {
			rl = rlx_cp[i] * dy0 - rly_cp[i] * dx0;
		}
		else {
			rl = rlx_cp[i] * Ty2d[i] - rly_cp[i] * Tx2d[i];
		}
		if (rl >= 0.0) {
			ret = -1;
			break;
		}

		ip = rrx_cp[i] * Tx2d[i] + rry_cp[i] * Ty2d[i];
		if (ip > 0.0 && i < Xeidx) {
			rl = rrx_cp[i] * dy1 - rry_cp[i] * dx1;
		}
		else if (ip < 0.0 && Xsidx < i) {
			rl = rrx_cp[i] * dy0 - rry_cp[i] * dx0;
		}
		else {
			rl = rrx_cp[i] * Ty2d[i - 1] - rry_cp[i] * Tx2d[i - 1];
		}
		if (rl < 0.0) {
			ret = 1;
			break;
		}
	}

	return ret;
}

int crease::checkRulCreaseCross()
{
	int ret = 0;

	for (int i = Xsidx; i <= Xeidx; i++) {
		double lex, ley, rex, rey;
		lex = Xx2d[i] + rllen[i] * rlx_cp[i];
		ley = Xy2d[i] + rllen[i] * rly_cp[i];
		rex = Xx2d[i] + rrlen[i] * rrx_cp[i];
		rey = Xy2d[i] + rrlen[i] * rry_cp[i];
		for (int j = Xsidx; j < Xeidx; j++) {
			double ix, iy, l0, l1;
			int res = intersectionOfLine(Xx2d[i], Xy2d[i], lex, ley,
				Xx2d[j], Xy2d[j], Xx2d[j + 1], Xy2d[j + 1], &ix, &iy, &l0, &l1);
			if (res == 0 && l0 > 0.0 && l1 > 0.0) {
				ret = -1;
				break;
			}
		}
		if (ret < 0) {
			break;
		}
		for (int j = Xsidx; j < Xeidx; j++) {
			double ix, iy, l0, l1;
			int res = intersectionOfLine(Xx2d[i], Xy2d[i], rex, rey,
				Xx2d[j], Xy2d[j], Xx2d[j + 1], Xy2d[j + 1], &ix, &iy, &l0, &l1);
			if (res == 0 && l0 > 0.0 && l1 > 0.0) {
				ret = -1;
				break;
			}
		}
		if (ret < 0) {
			break;
		}
	}

	return ret;
}
