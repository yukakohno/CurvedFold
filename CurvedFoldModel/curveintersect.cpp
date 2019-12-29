#include <memory.h>
#include "curveintersect.h"

curveintersect::curveintersect()
{
	init();
}

curveintersect::~curveintersect()
{
	//init();
}

void curveintersect::init()
{
	ctype = CTYPE_UNDEF;
	cvcnt = 0;
	memset( cvx, 0, sizeof(double)*MAXFCX );
	memset( cvy, 0, sizeof(double)*MAXFCX );
	for( int i=0; i<MAXFCX; i++ ){
		cvtype[i] = CVTYPE_UNDEF;
	}
	memset( cvidx, 0, sizeof(int)*MAXFCX );
	memset( cvlen, 0, sizeof(double)*MAXFCX );
}

int curveintersect::create_drawcurve( curvedraw *dc, crease *c, int psx, int psy, int pex, int pey )
{
	int ret=0;

	int IScnt=0;
	double ISx[MAXFCX], ISy[MAXFCX];
	cvertex_type Xtype[MAXFCX];
	int Xidx[MAXFCX], Cidx[MAXFCX], ISodr[MAXFCX];
	double Xlen[MAXFCX];

	ret = curvedraw::IS_drawcurve( dc, c, psx, psy, pex, pey, ISx, ISy, MAXFCX, IScnt, Xtype, Xidx, Xlen, Cidx );
	ret = curvedraw::IS_check_cw( dc, c, psx, psy, pex, pey, ISx, ISy, MAXFCX, IScnt, Xtype, Xidx, Xlen, Cidx );

	cvcnt = IScnt;
	if( ret>0 ){ // clockwise
		for( int i=0; i<IScnt; i++ ){
			cvx[i] = ISx[i];
			cvy[i] = ISy[i];
			cvtype[i] = Xtype[i];
			cvidx[i] = Xidx[i];
			cvlen[i] = Xlen[i];
		}
	} else { // ccw
		for( int i=0; i<IScnt; i++ ){
			int j=IScnt-1-i;
			cvx[i] = ISx[j];
			cvy[i] = ISy[j];
			cvtype[i] = Xtype[j];
			cvidx[i] = Xidx[j];
			cvlen[i] = Xlen[j];
		}
	}

end:
	return ret;
}

// ‡”Ô”½“]
void curveintersect::inverse_c()
{
	double cvx_[MAXFCX], cvy_[MAXFCX];
	cvertex_type cvtype_[MAXFCX];
	int cvidx_[MAXFCX];
	double cvlen_[MAXFCX];

	memcpy( cvx_, cvx, sizeof(double)*cvcnt );
	memcpy( cvy_, cvy, sizeof(double)*cvcnt );
	memcpy( cvtype_, cvtype, sizeof(cvertex_type)*cvcnt );
	memcpy( cvidx_, cvidx, sizeof(int)*cvcnt );
	memcpy( cvlen_, cvlen, sizeof(double)*cvcnt );

	for( int i=0 ; i<cvcnt; i++ ){
		int j = cvcnt-1-i;
		cvx[i] = cvx_[j];
		cvy[i] = cvy_[j];
		cvtype[i] = cvtype[j];
		cvidx[i] = cvidx[j];
		cvlen[i] = cvlen[j];
	}
}

void curveintersect::closecv( double th_d /* 5.0 mm */ )
{
	double th_d2 = th_d*th_d; // mm
	double sx = cvx[0];
	double sy = cvy[0];
	double ex = cvx[cvcnt-1];
	double ey = cvy[cvcnt-1];
	double dx=ex-sx, dy=ey-sy, d2=dx*dx+dy*dy;
	if( d2<th_d2 ){
		if( cvcnt<MAXFCX-1 ){
			cvx[cvcnt] = sx;
			cvy[cvcnt] = sy;
			cvcnt++;
		} else {
			cvx[cvcnt-1] = sx;
			cvy[cvcnt-1] = sy;
		}
	}
}
