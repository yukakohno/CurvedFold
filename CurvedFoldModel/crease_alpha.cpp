#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"

int crease::setA(double a)
{
	int i,ret=0;
	if( a<0 ){ a=0.0; }
	if( a>M_PI/2.0 ){ a=M_PI/2.0; }
	for( i=0; i<Xcnt; i++ ){
		alpha[i] = a;
	}
end:
	return ret;
}

// TNB -> α算出
int crease::calcAlpha( int flg_rectifyA, rectify_param *rp )
{
	int ret=0;
	ret = calcAlpha( Xcnt, Xx, Xy, Xz, kv, k2d, cosa, sina, tana, alpha, da, flg_rectifyA, rp );
#if 0
	{
		FILE *fp=fopen("output/dump_alpha.csv","w");
		if( fp ){
			fprintf(fp, "k2d[i], kv[i], alpha[i], cosa[i], sina[i], da[i]\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f,%f\n", k2d[i], kv[i], alpha[i]*180/M_PI, cosa[i], sina[i], da[i]*180/M_PI);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;

}

int crease::calcAlpha( int Xcnt0, double *Xx0, double *Xy0, double *Xz0, double *kv0, double *k2d0,
					  double *cosa0, double *sina0, double *tana0, double *alpha0, double *da0,
					  int flg_rectifyA, rectify_param *rp )
{
	int i, ret=0;
	memset(  cosa0, 0, sizeof(double)*Xcnt0 );
	memset(  sina0, 0, sizeof(double)*Xcnt0 );
	memset(  tana0, 0, sizeof(double)*Xcnt0 );
	memset( alpha0, 0, sizeof(double)*Xcnt0 );
	memset(    da0, 0, sizeof(double)*Xcnt0 );

	for( i=0; i<Xcnt0; i++ ){
		if( kv0[i]==0 ){
			continue;
		}
		cosa0[i] = k2d0[i]/kv0[i];
		if( cosa0[i]>1.0 ){ cosa0[i]=1.0; }
		if( cosa0[i]<-1.0 ){ cosa0[i]=-1.0; }
		alpha0[i] = acos(cosa0[i]);
	}
	if( flg_rectifyA && rp ){
		rp->scnt = gets1e1( Xcnt0, CEMGN, Xcnt0-CEMGN-1, k2d0, rp->kvthres, rp->s1mgn, rp->s1, rp->e1, RECT_ASIZE );
		if( rp->scnt>0 ){
			rectifyAlphaBezier( Xcnt0, alpha0, rp );
		}
	}
	for( i=0; i<Xcnt0; i++ ){
		cosa0[i] = cos(alpha0[i]);
		sina0[i] = sin(alpha0[i]);
		tana0[i] = sina0[i]/cosa0[i];
	}
	for( i=1; i<Xcnt0-1; i++ ){
		double dx0 = Xx0[i+1] - Xx0[i];
		double dx1 = Xx0[i]   - Xx0[i-1];
		double dy0 = Xy0[i+1] - Xy0[i];
		double dy1 = Xy0[i]   - Xy0[i-1];
		double dz0 = Xz0[i+1] - Xz0[i];
		double dz1 = Xz0[i]   - Xz0[i-1];
		double d0 = sqrt( dx0*dx0 + dy0*dy0 + dz0*dz0 );
		double d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
		double da0_ = alpha0[i+1] - alpha0[i];
		double da1_ = alpha0[i]   - alpha0[i-1];
		da0[i] = (da0_/d0 + da1_/d1)*0.5;
	}
end:
	return ret;

}

// alpha, kv -> sina, cosa, tana, da, k2d
int crease::calcAK_K2D()
{
	int i, ret=0;

	memset( cosa, 0, sizeof(double)*MAX_SPCNT );
	memset( sina, 0, sizeof(double)*MAX_SPCNT );
	memset( tana, 0, sizeof(double)*MAX_SPCNT );
	memset(   da, 0, sizeof(double)*MAX_SPCNT );
	memset(  k2d, 0, sizeof(double)*MAX_SPCNT );

	for( i=0; i<Xcnt; i++ ){
		cosa[i] = cos(alpha[i]);
		sina[i] = sin(alpha[i]); // -:山折り
		tana[i] = tan(alpha[i]);
		k2d[i] = kv[i]*cosa[i];
#if 0
		//if( Nx[i]*Nx[0] + Ny[i]*Ny[0] + Nz[i]*Nz[0] < 0  ){
		if( Nx[i]*0 + Ny[i]*1 + Nz[i]*0 < 0  ){
			k2d[i] = -k2d[i];
		}
#endif
	}
#if 1 // 距離に応じた3点微分
	for( i=1; i<Xcnt-1; i++ ){
		double dx0 = Xx[i+1] - Xx[i];
		double dx1 = Xx[i]   - Xx[i-1];
		double dy0 = Xy[i+1] - Xy[i];
		double dy1 = Xy[i]   - Xy[i-1];
		double dz0 = Xz[i+1] - Xz[i];
		double dz1 = Xz[i]   - Xz[i-1];
		double d0 = sqrt( dx0*dx0 + dy0*dy0 + dz0*dz0 );
		double d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
		double da0 = alpha[i+1] - alpha[i];
		double da1 = alpha[i]   - alpha[i-1];
		da[i] = (da0/d0 + da1/d1)*0.5;
	}
#else // 5点微分
	for( i=6; i<Xcnt-6; i++ ){
		da[i] = alpha[i-2] -8*alpha[i-1] +8*alpha[i+1] -alpha[i+2];
	}
#endif
#if 0
	{
		FILE *fp=fopen("output/dump_alpha.csv","w");
		if( fp ){
			fprintf(fp, "k2d[i], kv[i], alpha[i], cosa[i], sina[i], da[i]\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f,%f\n", k2d[i], kv[i], alpha[i]*180/M_PI, cosa[i], sina[i], da[i]*180/M_PI);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}

// alpha, k2d -> (rectifyAlpha) -> kv
int crease::calcAK2D_K()
{
	int i, ret=0;
	memset( cosa, 0, sizeof(double)*MAX_SPCNT );
	memset( sina, 0, sizeof(double)*MAX_SPCNT );
	memset( tana, 0, sizeof(double)*MAX_SPCNT );
	memset(   kv, 0, sizeof(double)*MAX_SPCNT );

	for( i=0; i<Xcnt; i++ ){
			cosa[i] = cos(alpha[i]);
			sina[i] = sin(alpha[i]); // -:山折り
			tana[i] = tan(alpha[i]);
			kv[i] = k2d[i]/cosa[i];
	}
#if 0
	{
		FILE *fp=fopen("output/dump_alpha.csv","w");
		if( fp ){
			fprintf(fp, "k2d[i], kv[i], alpha[i], cosa[i], sina[i], da[i]\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f\n", k2d[i], kv[i], alpha[i]*180/M_PI, cosa[i], sina[i]);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}

// alpha, X -> da
int crease::calcDA()
{
	int i, ret=0;

	memset( da, 0, sizeof(double)*MAX_SPCNT );

#if 1 // 距離に応じた3点微分
	for( i=1; i<Xcnt-1; i++ ){
		double dx0 = Xx[i+1] - Xx[i];
		double dx1 = Xx[i]   - Xx[i-1];
		double dy0 = Xy[i+1] - Xy[i];
		double dy1 = Xy[i]   - Xy[i-1];
		double dz0 = Xz[i+1] - Xz[i];
		double dz1 = Xz[i]   - Xz[i-1];
		double d0 = sqrt( dx0*dx0 + dy0*dy0 + dz0*dz0 );
		double d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
		double da0 = alpha[i+1] - alpha[i];
		double da1 = alpha[i]   - alpha[i-1];
		if( d0<=0.0 && d1<=0.0 ){
			da[i] = 0.0;
			continue;
		}
		if( d0>0.0 && d1>0.0 ){
			da[i] = (da0/d0 + da1/d1)*0.5;
		} else if( d0>0.0 ){
			da[i] = da0;
		} else {
			da[i] = da1;
		}
	}
#else // 5点微分
	for( i=6; i<Xcnt-6; i++ ){
		da[i] = alpha[i-2] -8*alpha[i-1] +8*alpha[i+1] -alpha[i+2];
	}
#endif
end:
	return ret;
}

#define MAX_ALPHA_LOOP 100
#define DIF_ALPHA 2.0/180.0*M_PI
#define MIN_ALPHA -90.0/180.0*M_PI
#define MAX_ALPHA 90.0/180.0*M_PI
//#define MIN_ALPHA -85.0/180.0*M_PI
//#define MAX_ALPHA 85.0/180.0*M_PI

// crease->tr,alpha算出
int crease::calcRul2TA( int retry_mode,	int mini, int maxi )
// retry_mode
// 0: quit and return ret!=0 
// 1: use same alpha and torsion at the edge, return ret!=0
// 2: change alpha and retry
// return 0: OK, -1/-2: NG
{
	int i, m=Xcnt/2, flg=0, ret=0;
	double datana[MAX_SPCNT], ttana[MAX_SPCNT], alpha_[MAX_SPCNT], da_[MAX_SPCNT], tr_[MAX_SPCNT];
	double tan_a, minalpha, maxalpha;

	if( mini < 0 ){ mini = 0; }
	if( maxi < 0 ){ maxi = Xcnt; }
	m = (maxi+mini)/2;

	memset( datana, 0, sizeof(double)*MAX_SPCNT ); 
	memset( ttana, 0, sizeof(double)*MAX_SPCNT );
	for( i=0; i<Xcnt; i++ ){
		if( cotbl[i]==0.0 || cotbr[i]==0.0 ){
			datana[i] = 0.0;
			ttana[i] = 0.0;
		} else {
			datana[i] = 0.5 * k2d[i] * (cotbl[i]-cotbr[i]);
			ttana[i] = 0.5 * k2d[i] * (cotbl[i]+cotbr[i]);
		}
		alpha_[i] = da_[i] = tr_[i] = 0.0;
	}

	if( alpha[m] > 0){
		minalpha = MIN_ALPHA;
		maxalpha = MAX_ALPHA;
	} else {
		minalpha = -MAX_ALPHA;
		maxalpha = -MIN_ALPHA;
	}

	flg=0;
	for( int c=0; c<MAX_ALPHA_LOOP; c++ ){
		//printf( "c = %d\n", c );
		switch( flg ){
			case 0:
				alpha_[m] = alpha[m];	// 今の値をそのまま使う（1回目）
				//printf( " alpha_[m=%d] = %f\n", m, alpha_[m] );
				break;
			case -1:
				alpha_[m] -= DIF_ALPHA;
				if( alpha_[m]>=minalpha ){
					//printf( " alpha_[m=%d] = %f\n", m, alpha_[m] );
					flg = 0;
				} else {
					//printf( " alpha_[m=%d] = %f < minalpha()\n", m, alpha_[m], minalpha );
				}
				break;
			case -2:
				alpha_[m] += DIF_ALPHA;
				if( alpha_[m]<maxalpha ){
					//printf( " alpha_[m=%d] = %f\n", m, alpha_[m] );
					flg = 0;
				} else {
					//printf( " alpha_[m=%d] = %f > maxalpha()\n", m, alpha_[m], maxalpha );
				}
				break;
		}
		if( flg!=0 ){
			ret = flg;
			break; // for(c)
		}

		tan_a = tan(alpha_[m]);
		da_[m] = datana[m]*tan_a;
		tr_[m] = ttana[m]*tan_a;

		// set range [mini, m-1]
		for( i=m-1; i>=mini; i-- ){
			alpha_[i] = alpha_[i+1] - da_[i+1]*XspcOrg;
			if( alpha_[i] >= maxalpha ){
				//printf( " alpha_[%d] = %f >= maxalpha(%f)\n", i, alpha_[i], maxalpha );
				flg = -1;	break; // for(i)
			} else if( alpha_[i] < minalpha ){
				//printf( " alpha_[%d] = %f < minalpha(%f)\n", i, alpha_[i], minalpha );
				flg = -2;	break; // for(i)
			} else {
				tan_a = tan(alpha_[i]);
				da_[i] = datana[i]*tan_a;
				tr_[i] = ttana[i]*tan_a;
			}
		}
		if( flg<0 ){
			// 0: quit and return ret!=0
			if (retry_mode == 0) {
				ret = flg;
				break; // for(c)
			}
			// 1: use same alpha and torsion at the edge, return ret!=0
			if (retry_mode == 1) {
				int ii = i + 1;
				for (; i >= mini; i--) {
					alpha_[i] = alpha_[ii];
					da_[i] = 0;
					tr_[i] = tr_[ii];
				}
				ret = flg;
				flg = 0;
			}
			// 2: change alpha and retry
			if (retry_mode == 2) {
				continue; // for(c)
			}
		}

		// set range [m+1, maxi]
		for( i=m+1; i<maxi; i++ ){
			alpha_[i] = alpha_[i-1] + da_[i-1]*XspcOrg;
			if( alpha_[i] >= maxalpha ){
				//printf( " alpha_[%d] = %f >= maxalpha(%f)\n", i, alpha_[i], maxalpha );
				flg = -1;	break; // for(i)
			} else if( alpha_[i] < minalpha ){
				//printf( " alpha_[%d] = %f < minalpha(%f)\n", i, alpha_[i], minalpha );
				flg = -2;	break; // for(i)
			} else {
				tan_a = tan(alpha_[i]);
				da_[i] = datana[i]*tan_a;
				tr_[i] = ttana[i]*tan_a;
			}
		}
		if( flg<0 ){
			// 0: quit and return ret!=0
			if (retry_mode == 0) {
				ret = flg;
				break; // for(c)
			}
			// 1: use same alpha and torsion at the edge, return ret!=0
			if (retry_mode == 1) {
				int ii = i - 1;
				for (; i < maxi; i++) {
					alpha_[i] = alpha_[ii];
					da_[i] = 0;
					tr_[i] = tr_[ii];
				}
				ret = flg;
				flg = 0;
				break; // for(c)
			}
			// 2: change alpha and retry
			if (retry_mode == 2) {
				continue; // for(c)
			}
		}
	}

#if 0
	FILE *fp = fopen("output/Rul2TA0.csv","w");
	if(fp){
		fprintf(fp, "cotbl,cotbr,cotbl-cotbr,cotbl+cotbr,k2d,datana,ttana,alpha,tana,da,tr,alpha1,da1,tr1\n");
		for( int i=0; i<Xcnt; i++ ){
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				cotbl[i], cotbr[i], cotbl[i]-cotbr[i], cotbl[i]+cotbr[i], k2d[i],
				datana[i], ttana[i], alpha[i], tan(alpha[i]), da[i], tr[i],
				alpha_[i], da_[i], tr_[i]);
		}
		fclose(fp); fp=NULL;
	}
#endif
	if( flg==0 ){
		// set outside of the range [mini, maxi]
		// とりあえず同じ値で埋める
		for (i = mini; i >= 0; i--) {
			alpha_[i] = alpha_[mini];
			da_[i] = 0;
			tr_[i] = tr_[mini];
		}
		for (i = maxi; i < Xcnt; i++) {
			alpha_[i] = alpha_[maxi - 1];
			da_[i] = 0;
			tr_[i] = tr_[maxi - 1];
		}
		memcpy( alpha, alpha_, sizeof(double)*MAX_SPCNT );
		memcpy( da, da_, sizeof(double)*MAX_SPCNT );
		memcpy( tr, tr_, sizeof(double)*MAX_SPCNT );
	} else {
		//printf("crease::calcRul2TA(): using original parameters\n");
		ret = flg;
	}
#if 0
	FILE *fp = fopen("output/Rul2TA.csv","w");
	if(fp){
		fprintf(fp, "cotbl,cotbr,cotbl-cotbr,cotbl+cotbr,k2d,datana,ttana,alpha,tana,da,tr\n");
		for( int i=0; i<Xcnt; i++ ){
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				cotbl[i], cotbr[i], cotbl[i]-cotbr[i], cotbl[i]+cotbr[i], k2d[i],
				datana[i], ttana[i], alpha[i], tan(alpha[i]), da[i], tr[i]);
		}
		fclose(fp); fp=NULL;
	}
#endif
	return ret;
}

