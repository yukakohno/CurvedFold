#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"
#include "util.h"

// creaseéZèo
//	double cosa[MAX_SPCNT], tana[MAX_SPCNT], da[MAX_SPCNT];
//	double cotbl[MAX_SPCNT], cotbr[MAX_SPCNT];
//	double cosbl[MAX_SPCNT], sinbl[MAX_SPCNT], cosbr[MAX_SPCNT], sinbr[MAX_SPCNT];
//	double rlx[MAX_SPCNT], rly[MAX_SPCNT], rlz[MAX_SPCNT];
//	double rrx[MAX_SPCNT], rry[MAX_SPCNT], rrz[MAX_SPCNT];
//	double rllen[MAX_SPCNT], rrlen[MAX_SPCNT]; // crease pattern ÇÇ‡Ç∆Ç…í∑Ç≥ãÅÇﬂÇÈ
int crease::calcRuling( int flg_rectifyR, double rectifyR_kvthres )
{
	int i, ret=0;
	memset( betal, 0, sizeof(double)*MAX_SPCNT );
	memset( betar, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbl, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbr, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbl, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbr, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbl, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbr, 0, sizeof(double)*MAX_SPCNT );
	memset(   rlx, 0, sizeof(double)*MAX_SPCNT );
	memset(   rly, 0, sizeof(double)*MAX_SPCNT );
	memset(   rlz, 0, sizeof(double)*MAX_SPCNT );
	memset(   rrx, 0, sizeof(double)*MAX_SPCNT );
	memset(   rry, 0, sizeof(double)*MAX_SPCNT );
	memset(   rrz, 0, sizeof(double)*MAX_SPCNT );
	memset( rllen, 0, sizeof(double)*MAX_SPCNT );
	memset( rrlen, 0, sizeof(double)*MAX_SPCNT );
	memset( rlx_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rly_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rrx_cp, 0, sizeof(double)*MAX_SPCNT );
	memset( rry_cp, 0, sizeof(double)*MAX_SPCNT );

	//double Bl[MAX_SPCNT],Br[MAX_SPCNT];
	//memset( Bl, 0, sizeof(double)*MAX_SPCNT );
	//memset( Br, 0, sizeof(double)*MAX_SPCNT );

	calcRuling2D( flg_rectifyR, rectifyR_kvthres );
	calcRuling3D();

	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		if( flg_rectifyR && fabs(k2d[i]) < rectifyR_kvthres ){
			rllen[i] = rrlen[i] = 0;
		} else{
			rllen[i] = rrlen[i] = maxrlen;
		}
	}
#if 0
	{
		FILE *fp=fopen("output/dump_ruling.csv","w");
		if( fp ){
			fprintf(fp, "dX2d,curv2d,dX,curv,tortion,A,sinA,cosA,tanA,dA,Bl,cotBl,sinBl,cosBl,Br,cotBr,sinBr,cosBr\n" );
			for( i=0; i<Xcnt; i++ ){

				fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					d2d[i],k2d[i], dx[i],kv[i],tr[i],
					alpha[i]*180/M_PI, sina[i],cosa[i],tana[i],da[i],
					betal[i]*180/M_PI, cotbl[i],sinbl[i],cosbl[i],
					betar[i]*180/M_PI, cotbr[i],sinbr[i],cosbr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	return ret;
}

int crease::calcRuling( int cxcnt0,
					   double *Tx0, double *Ty0, double *Tz0, double *Tx2d0, double *Ty2d0, 
					   double *Nx0, double *Ny0, double *Nz0,  double *Bx0, double *By0, double *Bz0,					   
					   double *da0, double *kv0, double *tr0,
					   double *betal0, double *betar0, double *sinbl0, double *sinbr0,
					   double *cosbl0, double *cosbr0, double *cotbl0, double *cotbr0, double *sina0, double *cosa0, 
					   double *rlx0, double *rly0, double *rlz0, double *rlx_cp0, double *rly_cp0, double *rllen0,
					   double *rrx0, double *rry0, double *rrz0, double *rrx_cp0, double *rry_cp0, double *rrlen0, int rl )
{
	int ret=0;
	double maxrlen = 400;
	memset( betal0, 0, sizeof(double)*MAX_SPCNT );
	memset( betar0, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbl0, 0, sizeof(double)*MAX_SPCNT );
	memset( cotbr0, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbl0, 0, sizeof(double)*MAX_SPCNT );
	memset( cosbr0, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbl0, 0, sizeof(double)*MAX_SPCNT );
	memset( sinbr0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rlx0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rly0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rlz0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rrx0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rry0, 0, sizeof(double)*MAX_SPCNT );
	memset(   rrz0, 0, sizeof(double)*MAX_SPCNT );
	//memset( rllen0, 0, sizeof(double)*MAX_SPCNT );
	//memset( rrlen0, 0, sizeof(double)*MAX_SPCNT );
	memset( rlx_cp0, 0, sizeof(double)*MAX_SPCNT );
	memset( rly_cp0, 0, sizeof(double)*MAX_SPCNT );
	memset( rrx_cp0, 0, sizeof(double)*MAX_SPCNT );
	memset( rry_cp0, 0, sizeof(double)*MAX_SPCNT );

	calcRuling2D(cxcnt0, da0, sina0, kv0, tr0,
		betal0, betar0, sinbl0, sinbr0, cosbl0, cosbr0, cotbl0, cotbr0, rl);
	calcRuling3D( cxcnt0, Tx0, Ty0, Tz0, Tx2d0, Ty2d0, Nx0, Ny0, Nz0, Bx0, By0, Bz0,
		sina0, cosa0, sinbl0, sinbr0, cosbl0, cosbr0,
		rlx0, rly0, rlz0, rrx0, rry0, rrz0, rlx_cp0, rly_cp0, rrx_cp0, rry_cp0, rl );
	for( int i=0; i<cxcnt0; i++ ){
		rllen0[i] = maxrlen;
		rrlen0[i] = maxrlen;
	}
	return ret;
}


int crease::calcRuling2D( int flg_rectifyR, double rectifyR_kvthres )
{
	int i, ret=0;

	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		if( flg_rectifyR && fabs(k2d[i]) < rectifyR_kvthres ){
			cotbl[i] = 0.0;
			cotbr[i] = 0.0;
			betal[i] = M_PI/2.0;
			betar[i] = M_PI/2.0;
			continue;
		}
		if( sina[i]*kv[i]==0.0 ){
			betal[i] = betar[i] = M_PI/2.0;
			cotbl[i] = cotbr[i] = 10000.0; // infinity, not used after this
			sinbl[i] = sinbr[i] = 1.0;
			cosbl[i] = cosbr[i] = 0.0;
			continue;
		}
			cotbl[i] = (da[i]+tr[i])/(sina[i]*kv[i]); // +:[0,pi/2], -:[pi/2,pi]
			cotbr[i] = (-da[i]+tr[i])/(sina[i]*kv[i]);
			betal[i] = atan(1.0/cotbl[i]);
			if( betal[i] < 0.0 ){
				betal[i] += M_PI;
			}
			betar[i] = atan(1.0/cotbr[i]);
			if( betar[i] < 0.0 ){
				betar[i] += M_PI;
			}
		sinbl[i] = 1.0/sqrt(1.0+cotbl[i]*cotbl[i]); // +:[0,pi]
		sinbr[i] = 1.0/sqrt(1.0+cotbr[i]*cotbr[i]);
		cosbl[i] = cotbl[i]/sqrt(1.0+cotbl[i]*cotbl[i]); // +:[0,pi/2], -:[pi/2,pi]
		cosbr[i] = cotbr[i]/sqrt(1.0+cotbr[i]*cotbr[i]);
	}
end:
	return ret;
}

int crease::calcRuling2D( int cxcnt0, double *da0, double *sina0, double *kv0, double *tr0,
						 double *betal0, double *betar0, double *sinbl0, double *sinbr0,
						 double *cosbl0, double *cosbr0, double *cotbl0, double *cotbr0, int rl )
{
	int ret=0;
	if( rl < 0){
		for( int i=0; i<cxcnt0; i++ ){
			if( sina0[i]==0.0 || kv0[i]==0.0 ){
#if 0
				cotbl0[i] = 0.0;
				cotbr0[i] = 0.0;
				betal0[i] = 0.0;
				betar0[i] = 0.0;
				sinbl0[i] = 0.0;
				sinbr0[i] = 0.0;
				cosbl0[i] = 0.0;
				cosbr0[i] = 0.0;
#endif
			} else {
				cotbl0[i] = (-da0[i]+tr0[i])/(sina0[i]*kv0[i]); // +:[0,pi/2], -:[pi/2,pi]
				cotbr0[i] = (da0[i]+tr0[i])/(sina0[i]*kv0[i]);
				betal0[i] = atan(1.0/cotbl0[i]);
				if( betal0[i] < 0.0 ){
					betal0[i] += M_PI;
				}
				betar0[i] = atan(1.0/cotbr0[i]);
				if( betar0[i] < 0.0 ){
					betar0[i] += M_PI;
				}
				sinbl0[i] = 1.0/sqrt(1.0+cotbl0[i]*cotbl0[i]); // +:[0,pi]
				sinbr0[i] = 1.0/sqrt(1.0+cotbr0[i]*cotbr0[i]);
				cosbl0[i] = cotbl0[i]/sqrt(1.0+cotbl0[i]*cotbl0[i]); // +:[0,pi/2], -:[pi/2,pi]
				cosbr0[i] = cotbr0[i]/sqrt(1.0+cotbr0[i]*cotbr0[i]);
			}
		}
	} else {
		for( int i=0; i<cxcnt0; i++ ){
			if( sina0[i]==0.0 || kv0[i]==0.0 ){
			} else {
				cotbl0[i] = (da0[i]+tr0[i])/(sina0[i]*kv0[i]); // +:[0,pi/2], -:[pi/2,pi]
				cotbr0[i] = (-da0[i]+tr0[i])/(sina0[i]*kv0[i]);
				betal0[i] = atan(1.0/cotbl0[i]);
				if( betal0[i] < 0.0 ){
					betal0[i] += M_PI;
				}
				betar0[i] = atan(1.0/cotbr0[i]);
				if( betar0[i] < 0.0 ){
					betar0[i] += M_PI;
				}
				sinbl0[i] = 1.0/sqrt(1.0+cotbl0[i]*cotbl0[i]); // +:[0,pi]
				sinbr0[i] = 1.0/sqrt(1.0+cotbr0[i]*cotbr0[i]);
				cosbl0[i] = cotbl0[i]/sqrt(1.0+cotbl0[i]*cotbl0[i]); // +:[0,pi/2], -:[pi/2,pi]
				cosbr0[i] = cotbr0[i]/sqrt(1.0+cotbr0[i]*cotbr0[i]);
			}
		}
	}
end:
	return ret;
}

int crease::calcRuling2DH( int flg_rectifyR, double rectifyR_kvthres )
{
	int i, ret=0;
	double rx,ry, dx,dy,len,ip,op;
	//rx=1.0;	ry=0.0;
	rx=cos(M_PI*0.25);	ry=-sin(M_PI*0.25);
	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		if( flg_rectifyR && fabs(k2d[i]) < rectifyR_kvthres ){
			cotbl[i] = cotbr[i] = 0.0;
			cosbl[i] = cosbr[i] = 0.0;
			sinbl[i] = sinbr[i] = 1.0;
			betal[i] = betar[i] = M_PI/2.0;
		} else {
			dx = Xx2d[i+1]-Xx2d[i];
			dy = Xy2d[i+1]-Xy2d[i];
			len = sqrt(dx*dx+dy*dy); dx/=len; dy/=len;
			ip = dx*rx+dy*ry;
			op = dx*ry-dy*rx;
			if( op<0 ){ 
				cosbl[i] = -ip;
				cosbr[i] = ip;
			} else {
				cosbl[i] = ip;
				cosbr[i] = -ip;
			}
			sinbl[i] = sinbr[i] = sqrt(1.0-ip*ip);
			cotbl[i] = cosbl[i]/sinbl[i];
			cotbr[i] = cosbr[i]/sinbr[i];

			betal[i] = atan(1.0/cotbl[i]);
			if( betal[i] < 0.0 ){
				betal[i] += M_PI;
			}
			betar[i] = atan(1.0/cotbr[i]);
			if( betar[i] < 0.0 ){
				betar[i] += M_PI;
			}
		}
	}
end:
	return ret;
}

int crease::calcRuling3D()
{
	int i, ret=0;
	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		rlx[i] = cosbl[i]*Tx[i] +sinbl[i]*cosa[i]*Nx[i] +sinbl[i]*sina[i]*Bx[i];
		rly[i] = cosbl[i]*Ty[i] +sinbl[i]*cosa[i]*Ny[i] +sinbl[i]*sina[i]*By[i];
		rlz[i] = cosbl[i]*Tz[i] +sinbl[i]*cosa[i]*Nz[i] +sinbl[i]*sina[i]*Bz[i];
		normalize_v3( &rlx[i], &rly[i], &rlz[i] );
		rrx[i] = cosbr[i]*Tx[i] -sinbr[i]*cosa[i]*Nx[i] +sinbr[i]*sina[i]*Bx[i];
		rry[i] = cosbr[i]*Ty[i] -sinbr[i]*cosa[i]*Ny[i] +sinbr[i]*sina[i]*By[i];
		rrz[i] = cosbr[i]*Tz[i] -sinbr[i]*cosa[i]*Nz[i] +sinbr[i]*sina[i]*Bz[i];
		normalize_v3( &rrx[i], &rry[i], &rrz[i] );
		rlx_cp[i] = cosbl[i]*Tx2d[i] -sinbl[i]*Ty2d[i];
		rly_cp[i] = sinbl[i]*Tx2d[i] +cosbl[i]*Ty2d[i];
		normalize_v2( &rlx_cp[i], &rly_cp[i] );
		rrx_cp[i] = cosbr[i]*Tx2d[i] +sinbr[i]*Ty2d[i];
		rry_cp[i] = -sinbr[i]*Tx2d[i] + cosbr[i]*Ty2d[i];
		normalize_v2( &rrx_cp[i], &rry_cp[i] );
	}
end:
	return ret;
}

int crease::calcRuling3D( int cvcnt0, double *Tx0, double *Ty0, double *Tz0, double *Tx2d0, double *Ty2d0, 
						 double *Nx0, double *Ny0, double *Nz0,  double *Bx0, double *By0, double *Bz0,
						 double *sina0, double *cosa0, double *sinbl0, double *sinbr0, double *cosbl0, double *cosbr0,
						 double *rlx0, double *rly0, double *rlz0, double *rrx0, double *rry0, double *rrz0,
						 double *rlx_cp0, double *rly_cp0, double *rrx_cp0, double *rry_cp0, int rl )
{
	int i, ret=0;
	if( rl < 0){
		for( i=0; i<cvcnt0; i++ ){
			rlx0[i] = cosbl0[i]*Tx0[i] -sinbl0[i]*cosa0[i]*Nx0[i] +sinbl0[i]*sina0[i]*Bx0[i];
			rly0[i] = cosbl0[i]*Ty0[i] -sinbl0[i]*cosa0[i]*Ny0[i] +sinbl0[i]*sina0[i]*By0[i];
			rlz0[i] = cosbl0[i]*Tz0[i] -sinbl0[i]*cosa0[i]*Nz0[i] +sinbl0[i]*sina0[i]*Bz0[i];
			normalize_v3( &rlx0[i], &rly0[i], &rlz0[i] );
			rrx0[i] = cosbr0[i]*Tx0[i] +sinbr0[i]*cosa0[i]*Nx0[i] +sinbr0[i]*sina0[i]*Bx0[i];
			rry0[i] = cosbr0[i]*Ty0[i] +sinbr0[i]*cosa0[i]*Ny0[i] +sinbr0[i]*sina0[i]*By0[i];
			rrz0[i] = cosbr0[i]*Tz0[i] +sinbr0[i]*cosa0[i]*Nz0[i] +sinbr0[i]*sina0[i]*Bz0[i];
			normalize_v3( &rrx0[i], &rry0[i], &rrz0[i] );
			rlx_cp0[i] = cosbl0[i]*Tx2d0[i] -sinbl0[i]*Ty2d0[i];
			rly_cp0[i] = sinbl0[i]*Tx2d0[i] +cosbl0[i]*Ty2d0[i];
			normalize_v2( &rlx_cp0[i], &rly_cp0[i] );
			rrx_cp0[i] = cosbr0[i]*Tx2d0[i] +sinbr0[i]*Ty2d0[i];
			rry_cp0[i] = -sinbr0[i]*Tx2d0[i] + cosbr0[i]*Ty2d0[i];
			normalize_v2( &rrx_cp0[i], &rry_cp0[i] );
		}
	} else {
		for( i=0; i<cvcnt0; i++ ){
			rlx0[i] = cosbl0[i]*Tx0[i] +sinbl0[i]*cosa0[i]*Nx0[i] +sinbl0[i]*sina0[i]*Bx0[i];
			rly0[i] = cosbl0[i]*Ty0[i] +sinbl0[i]*cosa0[i]*Ny0[i] +sinbl0[i]*sina0[i]*By0[i];
			rlz0[i] = cosbl0[i]*Tz0[i] +sinbl0[i]*cosa0[i]*Nz0[i] +sinbl0[i]*sina0[i]*Bz0[i];
			normalize_v3( &rlx0[i], &rly0[i], &rlz0[i] );
			rrx0[i] = cosbr0[i]*Tx0[i] -sinbr0[i]*cosa0[i]*Nx0[i] +sinbr0[i]*sina0[i]*Bx0[i];
			rry0[i] = cosbr0[i]*Ty0[i] -sinbr0[i]*cosa0[i]*Ny0[i] +sinbr0[i]*sina0[i]*By0[i];
			rrz0[i] = cosbr0[i]*Tz0[i] -sinbr0[i]*cosa0[i]*Nz0[i] +sinbr0[i]*sina0[i]*Bz0[i];
			normalize_v3( &rrx0[i], &rry0[i], &rrz0[i] );
			rlx_cp0[i] = cosbl0[i]*Tx2d0[i] -sinbl0[i]*Ty2d0[i];
			rly_cp0[i] = sinbl0[i]*Tx2d0[i] +cosbl0[i]*Ty2d0[i];
			normalize_v2( &rlx_cp0[i], &rly_cp0[i] );
			rrx_cp0[i] = cosbr0[i]*Tx2d0[i] +sinbr0[i]*Ty2d0[i];
			rry_cp0[i] = -sinbr0[i]*Tx2d0[i] + cosbr0[i]*Ty2d0[i];
			normalize_v2( &rrx_cp0[i], &rry_cp0[i] );
		}
	}
end:
	return ret;
}

