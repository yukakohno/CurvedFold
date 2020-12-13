#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "crease.h"
#include "util.h"


int crease::getSamplePoints2D( double space, double *x, double *y )
{
	int pcnt=0;
	double len_nextp=space, len_prevv=0.0, len_nextv;
	for( int i=0; i<Xeidx; i++ )
	{
		len_nextv = len_prevv + d2d[i];
		while( len_nextp < len_nextv ){
			double mod = (len_nextp-len_prevv) / d2d[i];
			x[pcnt] = Xx2d[i] + mod * (Xx2d[i+1]-Xx2d[i]);
			y[pcnt] = Xy2d[i] + mod * (Xy2d[i+1]-Xy2d[i]);
			pcnt++;
			len_nextp += space;
		}
		len_prevv = len_nextv;
	} 
	return pcnt;
}

int crease::getSamplePoints3D( double space, double *x, double *y, double *z, double *len_sp )
{
	int pcnt=0;
	double len_nextp=space, len_prevv=0.0, len_nextv;
	double dxs, dxe;

	{
		double dx,dy,dz;
		dx = Xx[Xsidx]-Xxs;	dy = Xy[Xsidx]-Xys;	dz = Xz[Xsidx]-Xzs;	dxs = sqrt(dx*dx+dy*dy+dz*dz);
		dx = Xxe-Xx[Xeidx];	dy = Xye-Xy[Xeidx];	dz = Xze-Xz[Xeidx];	dxe = sqrt(dx*dx+dy*dy+dz*dz);
	}

	len_nextv = len_prevv + dxs;	// 次の頂点
	while( len_nextp < len_nextv ){
		double mod = (len_nextp-len_prevv) / dxs;
		x[pcnt] = Xxs + mod * (Xx[Xsidx]-Xxs);
		y[pcnt] = Xys + mod * (Xy[Xsidx]-Xys);
		z[pcnt] = Xzs + mod * (Xz[Xsidx]-Xzs);
		if( len_sp ){ len_sp[pcnt] = len_nextp; }
		pcnt++;
		len_nextp += space;	// 次のサンプル点
	}
	len_prevv = len_nextv;

	for( int i=Xsidx; i<Xeidx; i++ )
	{
		len_nextv = len_prevv + dx[i];	// 次の頂点
		while( len_nextp < len_nextv ){
			double mod = (len_nextp-len_prevv) / dx[i];
			x[pcnt] = Xx[i] + mod * (Xx[i+1]-Xx[i]);
			y[pcnt] = Xy[i] + mod * (Xy[i+1]-Xy[i]);
			z[pcnt] = Xz[i] + mod * (Xz[i+1]-Xz[i]);
			if( len_sp ){ len_sp[pcnt] = len_nextp; }
			pcnt++;
			len_nextp += space;	// 次のサンプル点
		}
		len_prevv = len_nextv;
	}

	len_nextv = len_prevv + dxe;	// 次の頂点
	while( len_nextp < len_nextv ){
		double mod = (len_nextp-len_prevv) / dxe;
		x[pcnt] = Xx[Xeidx] + mod * (Xxe-Xx[Xeidx]);
		y[pcnt] = Xy[Xeidx] + mod * (Xye-Xy[Xeidx]);
		z[pcnt] = Xz[Xeidx] + mod * (Xze-Xz[Xeidx]);
		if( len_sp ){ len_sp[pcnt] = len_nextp; }
		pcnt++;
		len_nextp += space;	// 次のサンプル点
	}
	len_prevv = len_nextv;

	return pcnt;
}

int crease::getSamplePoints2D3D(double space, double* x, double* y, double* z, double* x2d, double* y2d, double* len_sp)
{
	int pcnt = 0;
	double len_nextp = space, len_prevv = 0.0, len_nextv;
	double dxs, dxe;

	{
		double dx, dy, dz;
		dx = Xx[Xsidx] - Xxs;	dy = Xy[Xsidx] - Xys;	dz = Xz[Xsidx] - Xzs;	dxs = sqrt(dx * dx + dy * dy + dz * dz);
		dx = Xxe - Xx[Xeidx];	dy = Xye - Xy[Xeidx];	dz = Xze - Xz[Xeidx];	dxe = sqrt(dx * dx + dy * dy + dz * dz);
	}

	len_nextv = len_prevv + dxs;	// 次の頂点
	while (len_nextp < len_nextv) {
		double mod = (len_nextp - len_prevv) / dxs;
		x[pcnt] = Xxs + mod * (Xx[Xsidx] - Xxs);
		y[pcnt] = Xys + mod * (Xy[Xsidx] - Xys);
		z[pcnt] = Xzs + mod * (Xz[Xsidx] - Xzs);
		x2d[pcnt] = Xxs2d + mod * (Xx2d[Xsidx] - Xxs2d);
		y2d[pcnt] = Xys2d + mod * (Xy2d[Xsidx] - Xys2d);
		if (len_sp) { len_sp[pcnt] = len_nextp; }
		pcnt++;
		len_nextp += space;	// 次のサンプル点
	}
	len_prevv = len_nextv;

	for (int i = Xsidx; i < Xeidx; i++)
	{
		len_nextv = len_prevv + dx[i];	// 次の頂点
		while (len_nextp < len_nextv) {
			double mod = (len_nextp - len_prevv) / dx[i];
			x[pcnt] = Xx[i] + mod * (Xx[i + 1] - Xx[i]);
			y[pcnt] = Xy[i] + mod * (Xy[i + 1] - Xy[i]);
			z[pcnt] = Xz[i] + mod * (Xz[i + 1] - Xz[i]);
			x2d[pcnt] = Xx2d[i] + mod * (Xx2d[i + 1] - Xx2d[i]);
			y2d[pcnt] = Xy2d[i] + mod * (Xy2d[i + 1] - Xy2d[i]);
			if (len_sp) { len_sp[pcnt] = len_nextp; }
			pcnt++;
			len_nextp += space;	// 次のサンプル点
		}
		len_prevv = len_nextv;
	}

	len_nextv = len_prevv + dxe;	// 次の頂点
	while (len_nextp < len_nextv) {
		double mod = (len_nextp - len_prevv) / dxe;
		x[pcnt] = Xx[Xeidx] + mod * (Xxe - Xx[Xeidx]);
		y[pcnt] = Xy[Xeidx] + mod * (Xye - Xy[Xeidx]);
		z[pcnt] = Xz[Xeidx] + mod * (Xze - Xz[Xeidx]);
		x2d[pcnt] = Xx2d[Xeidx] + mod * (Xxe2d - Xx2d[Xeidx]);
		y2d[pcnt] = Xy2d[Xeidx] + mod * (Xye2d - Xy2d[Xeidx]);
		if (len_sp) { len_sp[pcnt] = len_nextp; }
		pcnt++;
		len_nextp += space;	// 次のサンプル点
	}
	len_prevv = len_nextv;

	return pcnt;
}

int crease::getVertexPoints3D( double space, double *x, double *y, double *z, double *len_sp )
{
	int pcnt=0;
	double dxs, dxe;
	{
		double dx,dy,dz;
		dx = Xx[Xsidx]-Xxs;	dy = Xy[Xsidx]-Xys;	dz = Xz[Xsidx]-Xzs;	dxs = sqrt(dx*dx+dy*dy+dz*dz);
		dx = Xxe-Xx[Xeidx];	dy = Xye-Xy[Xeidx];	dz = Xze-Xz[Xeidx];	dxe = sqrt(dx*dx+dy*dy+dz*dz);
	}
	x[pcnt]=Xxs;		y[pcnt]=Xys;		z[pcnt]=Xzs;		len_sp[pcnt]=0.0;	pcnt++;
	x[pcnt]=Xx[Xsidx];	y[pcnt]=Xy[Xsidx];	z[pcnt]=Xz[Xsidx];	len_sp[pcnt]=dxs;	pcnt++;
	for( int i=Xsidx+1; i<=Xeidx; i++ ){
		x[pcnt]=Xx[i];	y[pcnt]=Xy[i];	z[pcnt]=Xz[i];
		len_sp[pcnt] = len_sp[pcnt-1] + dx[i];
		pcnt++;
	}
	x[pcnt]=Xxe;	y[pcnt]=Xye;	z[pcnt]=Xze;
	len_sp[pcnt]=len_sp[pcnt-1]+dxe;
	pcnt++;

	return pcnt;
}
