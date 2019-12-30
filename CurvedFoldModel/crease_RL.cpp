#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "crease.h"
#include "util.h"

int crease::setL2R( crease *c )
{
	int ret=0;
	Xcnt = c->Xcnt;

	Xsidx = -1;
	for( int i=0; i<Xcnt; i++ ){
		if( c->rlflg[i] == RETYPE_CURVE ){
			if( Xsidx<0 ){ Xsidx=i; }
			Xeidx = i;
			Xx2d[i] = c->Xx2d[i] + c->rlx_cp[i] * c->rllen[i];
			Xy2d[i] = c->Xy2d[i] + c->rly_cp[i] * c->rllen[i];
		}
	}
	// 曲線と右側rulingをコピー
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		Xx[i] = c->Xx[i] + c->rlx[i] * c->rllen[i];
		Xy[i] = c->Xy[i] + c->rly[i] * c->rllen[i];
		Xz[i] = c->Xz[i] + c->rlz[i] * c->rllen[i];
		Xx2d[i] = c->Xx2d[i] + c->rlx_cp[i] * c->rllen[i];
		Xy2d[i] = c->Xy2d[i] + c->rly_cp[i] * c->rllen[i];

		rrlen[i] = c->rllen[i];
		rrflg[i] = RETYPE_CURVE;
		rrx[i] = -c->rlx[i];
		rry[i] = -c->rly[i];
		rrz[i] = -c->rlz[i];
		rrx_cp[i] = -c->rlx_cp[i];
		rry_cp[i] = -c->rly_cp[i];
	}
	return ret;
}

int crease::setR2L( crease *c )
{
	int ret=0;
	Xcnt = c->Xcnt;

	Xsidx = -1;
	for( int i=0; i<Xcnt; i++ ){
		if( c->rrflg[i]==RETYPE_CURVE ){
			if( Xsidx<0 ){ Xsidx=i; }
			Xeidx = i;
			Xx2d[i] = c->Xx2d[i] + c->rrx_cp[i] * c->rrlen[i];
			Xy2d[i] = c->Xy2d[i] + c->rry_cp[i] * c->rrlen[i];
		}
	}
	// 曲線と右側rulingをコピー
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		Xx[i] = c->Xx[i] + c->rrx[i] * c->rrlen[i];
		Xy[i] = c->Xy[i] + c->rry[i] * c->rrlen[i];
		Xz[i] = c->Xz[i] + c->rrz[i] * c->rrlen[i];
		Xx2d[i] = c->Xx2d[i] + c->rrx_cp[i] * c->rrlen[i];
		Xy2d[i] = c->Xy2d[i] + c->rry_cp[i] * c->rrlen[i];

		rllen[i] = c->rrlen[i];
		rlflg[i] = RETYPE_CURVE;
		rlx[i] = -c->rrx[i];
		rly[i] = -c->rry[i];
		rlz[i] = -c->rrz[i];
		rlx_cp[i] = -c->rrx_cp[i];
		rly_cp[i] = -c->rry_cp[i];
	}
	return ret;
}
#if 0
int crease::calcTNB()
{
	int size=Xeidx-Xsidx+1, ret=0;

	calcTN2d( size, &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
		&(Tx2d[Xsidx]), &(Ty2d[Xsidx]), &(d2d[Xsidx]),
		&(Nx2d[Xsidx]), &(Ny2d[Xsidx]), &(k2d[Xsidx]) );

	calcTNB( size, &(Xx[Xsidx]), &(Xy[Xsidx]), &(Xz[Xsidx]),
		&(Tx[Xsidx]), &(Ty[Xsidx]), &(Tz[Xsidx]), &(dx[Xsidx]),
		&(Nx[Xsidx]), &(Ny[Xsidx]), &(Nz[Xsidx]), &(kv[Xsidx]),
		&(Bx[Xsidx]), &(By[Xsidx]), &(Bz[Xsidx]), &(tr[Xsidx]), 0/*flg_rectifyT*/ );
end:
	return ret;
}
#endif
int crease::calcLeft0()
{
	int ret=0;
	maxrlen = 20; // 仮
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		alpha[i] = sina[i] = tana[MAX_SPCNT] = da[MAX_SPCNT] = 0.0;
		cosa[i] = 1.0; 
		rlx[i] = -rrx[i];	rly[i] = -rry[i];	rlz[i] = -rrz[i];
		rlx_cp[i] = -rrx_cp[i];	rly_cp[i] = -rry_cp[i];

		cosbr[i] = rrx_cp[i]*Tx2d[i] + rry_cp[i]*Ty2d[i];
		sinbr[i] = 1.0 - cosbr[i]*cosbr[i];
		cotbr[i] = cosbr[i]/sinbr[i];
		betar[i] = acos( cosbr[i] );

		cosbl[i] = rlx_cp[i]*Tx2d[i] + rly_cp[i]*Ty2d[i]; // inner product
		sinbl[i] = 1.0 - cosbl[i]*cosbl[i];
		cotbl[i] = cosbl[i]/sinbl[i];
		betal[i] = acos( cosbl[i] );

		rlflg[i] = RETYPE_UNDEF;
		rllen[i] = 0; 		
	}

	return ret;
}

int crease::calcRight0()
{
	int ret=0;
	maxrlen = 20; // 仮
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		alpha[i] = sina[i] = tana[MAX_SPCNT] = da[MAX_SPCNT] = 0.0;
		cosa[i] = 1.0; 
		rrx[i] = -rlx[i];	rry[i] = -rly[i];	rrz[i] = -rlz[i];
		rrx_cp[i] = -rlx_cp[i];	rry_cp[i] = -rly_cp[i];

		cosbr[i] = Tx2d[i]*rrx_cp[i] + Ty2d[i]*rry_cp[i]; // inner product
		sinbr[i] = 1.0 - cosbr[i]*cosbr[i];
		betar[i] = acos( cosbr[i] );
		cotbr[i] = cosbr[i]/sinbr[i];

		cosbl[i] = Tx2d[i]*rrx_cp[i] + Ty2d[i]*rry_cp[i]; // inner product
		sinbl[i] = 1.0 - cosbl[i]*cosbl[i];
		betal[i] = acos( cosbl[i] );
		cotbl[i] = cosbl[i]/sinbl[i];

		rrflg[i] = RETYPE_UNDEF;
		rrlen[i] = 20; 
	}

	return ret;
}

int crease::calcLeft1()
{
	int ret=0;
	maxrlen = 20; // 仮

	int Xcnt0 = Xeidx - Xsidx +1;
	calcTNB( Xcnt, Xx, Xy, Xz, Tx, Ty, Tz, dx, Nx, Ny, Nz, kv, Bx, By, Bz, tr, 0 ); // flg_rectifyT
	calcTN2d( Xcnt, Xx2d, Xy2d, Tx2d, Ty2d, d2d, Nx2d, Ny2d, k2d );
	calcAlpha( Xcnt, Xx, Xy, Xz, kv, k2d, cosa, sina, tana, alpha, da, 0 ); // flg_rectifyA
	calcRuling( Xcnt, Tx, Ty, Tz, Tx2d, Ty2d, Nx, Ny, Nz, Bx, By, Bz, da, kv, tr,
		betal, betar, sinbl, sinbr, cosbl, cosbr, cotbl, cotbr, sina, cosa, 
		rlx, rly, rlz, rlx_cp, rly_cp, rllen, rrx, rry, rrz, rrx_cp, rry_cp, rrlen, -1 );

	return ret;
}

int crease::calcLeft()
{
	int ret=0;
	maxrlen = 400; // 仮
	//	double alpha[MAX_SPCNT], cosa[MAX_SPCNT], sina[MAX_SPCNT], tana[MAX_SPCNT], da[MAX_SPCNT];
	//	double betal[MAX_SPCNT], cotbl[MAX_SPCNT], cosbl[MAX_SPCNT], sinbl[MAX_SPCNT];
	//	double rlx[MAX_SPCNT], rly[MAX_SPCNT], rlz[MAX_SPCNT], rlx_cp[MAX_SPCNT], rly_cp[MAX_SPCNT];
	//	int rlflg[MAX_SPCNT]; // -1:未設定, 0:曲線, 1:紙端上, 2:紙端右, 3:紙端下, 4:紙端左
	//	double rllen[MAX_SPCNT]; // crease pattern をもとに長さ求める
	//	double maxrlen;
	// (rrx,rry,rrz), TNB → α
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		if( Nx[i]==0.0 && Ny[i]==0.0 && Nz[i]==0.0 ){
			cosa[i] = sina[i] = tana[i] = alpha[i] = 0.0;
		} else {
			// BN面に射影
			double ip,len,bnx,bny,bnz;
			ip = rrx[i]*Tx[i] + rry[i]*Ty[i] + rrz[i]*Tz[i];
			bnx = rrx[i] - ip*Tx[i];
			bny = rry[i] - ip*Ty[i];
			bnz = rrz[i] - ip*Tz[i];
			len = sqrt( bnx*bnx + bny*bny + bnz*bnz );
			bnx /= len;	bny /= len;	bnz /= len;
			//cosa[i] = bnx*Nx[i] + bny*Ny[i] + bnz*Nz[i];
			//sina[i] = 1.0 - cosa[i]*cosa[i];
			cosa[i] = bnx*(-Nx[i]) + bny*(-Ny[i]) + bnz*(-Nz[i]); // r なので -N
			sina[i] = bnx*Bx[i] + bny*By[i] + bnz*Bz[i];
			tana[i] = sina[i]/cosa[i];
			alpha[i] = atan(tana[i]);
		}
		if( Nx2d[i]==0.0 && Ny2d[i]==0.0 ){
			cosbr[i] = sinbr[i] = cotbr[i] = betar[i] = 0.0;
		} else {
			cosbr[i] = rrx_cp[i]*Tx2d[i] + rry_cp[i]*Ty2d[i];
			//sinbr[i] = 1.0 - cosbr[i]*cosbr[i];
			sinbr[i] = rrx_cp[i]*(-Nx2d[i]) + rry_cp[i]*(-Ny2d[i]); // r なので -N
			cotbr[i] = cosbr[i]/sinbr[i];
			betar[i] = acos( cosbr[i] );
		}
	}

	// da: 距離に応じた3点微分
	memset( da, 0, sizeof(double)*MAX_SPCNT );
	for( int i=Xsidx+1; i<Xeidx; i++ ){
		if( cosa[i-1]==0.0 && sina[i-1]==0.0
			|| cosa[i  ]==0.0 && sina[i  ]==0.0
			|| cosa[i+1]==0.0 && sina[i+1]==0.0 )
		{
			da[i] = 0.0;
			continue;
		}
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

	double betar0[MAX_SPCNT], cotbr0[MAX_SPCNT], cosbr0[MAX_SPCNT], sinbr0[MAX_SPCNT];
	double rrx0[MAX_SPCNT], rry0[MAX_SPCNT], rrz0[MAX_SPCNT], rrx_cp0[MAX_SPCNT], rry_cp0[MAX_SPCNT];
	memcpy( betar0, betar, sizeof(double)*MAX_SPCNT );
	memcpy( cotbr0, cotbr, sizeof(double)*MAX_SPCNT );
	memcpy( cosbr0, cosbr, sizeof(double)*MAX_SPCNT );
	memcpy( sinbr0, sinbr, sizeof(double)*MAX_SPCNT );
	memcpy( rrx0, rrx, sizeof(double)*MAX_SPCNT );
	memcpy( rry0, rry, sizeof(double)*MAX_SPCNT );
	memcpy( rrz0, rrz, sizeof(double)*MAX_SPCNT );
	memcpy( rrx_cp0, rrx_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rry_cp0, rry_cp, sizeof(double)*MAX_SPCNT );

	// TNBα → cotbl, cotbr
	for( int i=Xsidx+1; i<Xeidx; i++ ){
		if( cosa[i-1]==0.0 && sina[i-1]==0.0
			|| cosa[i  ]==0.0 && sina[i  ]==0.0
			|| cosa[i+1]==0.0 && sina[i+1]==0.0 )
		{
			betal[i] = betar[i] = sinbl[i] = sinbr[i] = cosbl[i] = cosbr[i] = 0.0;
			continue;
		}
		//cotbl[i] = (da[i]-tr[i])/(sina[i]*kv[i]); // +:[0,pi/2], -:[pi/2,pi]
		//cotbr[i] = (-da[i]-tr[i])/(sina[i]*kv[i]);
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

		rlflg[i] = RETYPE_UNDEF;
		rllen[i] = maxrlen; 		
	}
#if 0
	{
		FILE *fp=NULL;
		fp = fopen("output/calcleft_betar.csv","w");
		if( fp ){
			//cotbl[i] = (da[i]+tr[i])/(sina[i]*kv[i]); // +:[0,pi/2], -:[pi/2,pi]
			//cotbr[i] = (-da[i]+tr[i])/(sina[i]*kv[i]);
			fprintf(fp, "alpha,alpha_deg,da,tr,da+tr,'-da+tr,sina,kv,cotbl,bl,bl_deg,cotbr,br,br_deg\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					alpha[i], alpha[i]*180.0/M_PI, da[i], tr[i], da[i]+tr[i], -da[i]+tr[i], sina[i], kv[i],
					cotbl[i], betal[i], betal[i]*180.0/M_PI, cotbr[i], betar[i], betar[i]*180.0/M_PI);
			}
			fclose( fp ); fp = NULL;
		}
		fp = fopen("output/calcleft_alpha.csv","w");
		if( fp ){
			fprintf(fp, "alpha,alpha_deg,cosa,sina,tana,da,k2d,kv,tr\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					alpha[i], alpha[i]*180.0/M_PI, cosa[i], sina[i], tana[i], da[i], k2d[i], kv[i], tr[i] );
			}
			fclose( fp ); fp = NULL;
		}
		// check cotbr
		fp = fopen("output/calcleft.csv","w");
		if( fp ){
			fprintf(fp, "betar,betar_deg,cotbr,cosbr,sinbr,betar0,betar_deg0,cotbr0,cosbr0,sinbr0,diff_betar,,rrx, rry, rrz, rrx0, rry0, rrz0, ip_rr*,, rrx_cp, rry_cp,rrx_cp0, rry_cp0, ip_rr*_cp\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f\n",
					betar[i], betar[i]*180.0/M_PI, cotbr[i], cosbr[i], sinbr[i],
					betar0[i], betar0[i]*180.0/M_PI, cotbr0[i], cosbr0[i], sinbr0[i], (betar[i]-betar0[i])*180.0/M_PI,
					rrx[i], rry[i], rrz[i], rrx0[i], rry0[i], rrz0[i], rrx[i]*rrx0[i] + rry[i]*rry0[i] + rrz[i]*rrz0[i],
					rrx_cp[i], rry_cp[i],rrx_cp0[i], rry_cp0[i], rrx_cp[i]*rrx_cp0[i] + rry_cp[i]*rry_cp0[i] );
			}
			fclose( fp ); fp = NULL;
		}
	}
#endif
	memcpy( betar, betar0, sizeof(double)*MAX_SPCNT );
	memcpy( cotbr, cotbr0, sizeof(double)*MAX_SPCNT );
	memcpy( cosbr, cosbr0, sizeof(double)*MAX_SPCNT );
	memcpy( sinbr, sinbr0, sizeof(double)*MAX_SPCNT );
	memcpy( rrx, rrx0, sizeof(double)*MAX_SPCNT );
	memcpy( rry, rry0, sizeof(double)*MAX_SPCNT );
	memcpy( rrz, rrz0, sizeof(double)*MAX_SPCNT );
	memcpy( rrx_cp, rrx_cp0, sizeof(double)*MAX_SPCNT );
	memcpy( rry_cp, rry_cp0, sizeof(double)*MAX_SPCNT );

	return ret;
}
int crease::calcRight()
{
	int ret=0;
	maxrlen = 400; // 仮

	// (rrx,rry,rrz), TNB → α
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		if( Nx[i]==0.0 && Ny[i]==0.0 && Nz[i]==0.0 ){
			cosa[i] = sina[i] = tana[i] = alpha[i] = 0.0;
		} else {
			// BN面に射影
			double ip,len,bnx,bny,bnz;
			ip = rlx[i]*Tx[i] + rly[i]*Ty[i] + rlz[i]*Tz[i];
			bnx = rlx[i] - ip*Tx[i];
			bny = rly[i] - ip*Ty[i];
			bnz = rlz[i] - ip*Tz[i];
			len = sqrt( bnx*bnx + bny*bny + bnz*bnz );
			bnx /= len;	bny /= len;	bnz /= len;
			cosa[i] = bnx*Nx[i] + bny*Ny[i] + bnz*Nz[i];
			sina[i] = bnx*Bx[i] + bny*By[i] + bnz*Bz[i];
			tana[i] = sina[i]/cosa[i];
			alpha[i] = atan(tana[i]);
		}
		if( Nx2d[i]==0.0 && Ny2d[i]==0.0 ){
			cosbr[i] = sinbr[i] = cotbr[i] = betar[i] = 0.0;
		} else {
			cosbl[i] = rlx_cp[i]*Tx2d[i] + rly_cp[i]*Ty2d[i];
			sinbl[i] = rlx_cp[i]*Nx2d[i] + rly_cp[i]*Ny2d[i];
			cotbl[i] = cosbl[i]/sinbl[i];
			betal[i] = acos( cosbl[i] );
		}
	}

	// da: 距離に応じた3点微分
	memset( da, 0, sizeof(double)*MAX_SPCNT );
	for( int i=Xsidx+1; i<Xeidx; i++ ){
		if( cosa[i-1]==0.0 && sina[i-1]==0.0
			|| cosa[i  ]==0.0 && sina[i  ]==0.0
			|| cosa[i+1]==0.0 && sina[i+1]==0.0 )
		{
			da[i] = 0.0;
			continue;
		}
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

	double betal0[MAX_SPCNT], cotbl0[MAX_SPCNT], cosbl0[MAX_SPCNT], sinbl0[MAX_SPCNT];
	double rlx0[MAX_SPCNT], rly0[MAX_SPCNT], rlz0[MAX_SPCNT], rlx_cp0[MAX_SPCNT], rly_cp0[MAX_SPCNT];
	memcpy( betal0, betal, sizeof(double)*MAX_SPCNT );
	memcpy( cotbl0, cotbl, sizeof(double)*MAX_SPCNT );
	memcpy( cosbl0, cosbl, sizeof(double)*MAX_SPCNT );
	memcpy( sinbl0, sinbl, sizeof(double)*MAX_SPCNT );
	memcpy( rlx0, rlx, sizeof(double)*MAX_SPCNT );
	memcpy( rly0, rly, sizeof(double)*MAX_SPCNT );
	memcpy( rlz0, rlz, sizeof(double)*MAX_SPCNT );
	memcpy( rlx_cp0, rlx_cp, sizeof(double)*MAX_SPCNT );
	memcpy( rly_cp0, rly_cp, sizeof(double)*MAX_SPCNT );

	// TNBα → cotbl, cotbr
	for( int i=Xsidx+1; i<Xeidx; i++ ){
		if( cosa[i-1]==0.0 && sina[i-1]==0.0
			|| cosa[i  ]==0.0 && sina[i  ]==0.0
			|| cosa[i+1]==0.0 && sina[i+1]==0.0 )
		{
			betal[i] = betar[i] = sinbl[i] = sinbr[i] = cosbl[i] = cosbr[i] = 0.0;
			continue;
		}
		//cotbl[i] = (da[i]-tr[i])/(sina[i]*kv[i]); // +:[0,pi/2], -:[pi/2,pi]
		//cotbr[i] = (-da[i]-tr[i])/(sina[i]*kv[i]);
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

		rrflg[i] = RETYPE_UNDEF;
		rrlen[i] = maxrlen; 		
	}
#if 0
	{
		FILE *fp=NULL;
		fp = fopen("output/calcright_betar.csv","w");
		if( fp ){
			//cotbl[i] = (da[i]+tr[i])/(sina[i]*kv[i]); // +:[0,pi/2], -:[pi/2,pi]
			//cotbr[i] = (-da[i]+tr[i])/(sina[i]*kv[i]);
			fprintf(fp, "alpha,alpha_deg,da,tr,da+tr,'-da+trsina,kv,cotbl,bl,bl_deg,cotbr,br,br_deg\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					alpha[i], alpha[i]*180.0/M_PI, da[i], tr[i], da[i]+tr[i], -da[i]+tr[i], sina[i], kv[i],
					cotbl[i], betal[i], betal[i]*180.0/M_PI, cotbr[i], betar[i], betar[i]*180.0/M_PI);
			}
			fclose( fp ); fp = NULL;
		}
		fp = fopen("output/calcright_alpha.csv","w");
		if( fp ){
			fprintf(fp, "alpha,alpha_deg,cosa,sina,tana,da,k2d,kv,tr\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					alpha[i], alpha[i]*180.0/M_PI, cosa[i], sina[i], tana[i], da[i], k2d[i], kv[i], tr[i] );
			}
			fclose( fp ); fp = NULL;
		}
		// check cotbr
		fp = fopen("output/calcright.csv","w");
		if( fp ){
			fprintf(fp, "betal,betal_deg,cotbl,cosbl,sinbl,betal0,betal_deg0,cotbl0,cosbl0,sinbl0,diff_betal,,rlx, rly, rlz, rlx0, rly0, rlz0, ip_rl*,, rlx_cp, rly_cp,rlx_cp0, rly_cp0, ip_rl*_cp\n");
			for( int i=Xsidx+1; i<Xeidx; i++ ){
				fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f\n",
					betal[i], betal[i]*180.0/M_PI, cotbl[i], cosbl[i], sinbl[i],
					betal0[i], betal0[i]*180.0/M_PI, cotbl0[i], cosbl0[i], sinbl0[i], (betal[i]-betal0[i])*180.0/M_PI,
					rlx[i], rly[i], rlz[i], rlx0[i], rly0[i], rlz0[i], rlx[i]*rlx0[i] + rly[i]*rly0[i] + rlz[i]*rlz0[i],
					rlx_cp[i], rly_cp[i],rlx_cp0[i], rly_cp0[i], rlx_cp[i]*rlx_cp0[i] + rly_cp[i]*rly_cp0[i] );
			}
			fclose( fp ); fp = NULL;
		}
	}
#endif
	memcpy( betal, betal0, sizeof(double)*MAX_SPCNT );
	memcpy( cotbl, cotbl0, sizeof(double)*MAX_SPCNT );
	memcpy( cosbl, cosbl0, sizeof(double)*MAX_SPCNT );
	memcpy( sinbl, sinbl0, sizeof(double)*MAX_SPCNT );
	memcpy( rlx, rlx0, sizeof(double)*MAX_SPCNT );
	memcpy( rly, rly0, sizeof(double)*MAX_SPCNT );
	memcpy( rlz, rlz0, sizeof(double)*MAX_SPCNT );
	memcpy( rlx_cp, rlx_cp0, sizeof(double)*MAX_SPCNT );
	memcpy( rly_cp, rly_cp0, sizeof(double)*MAX_SPCNT );

	return ret;
}

int crease::calcLeftK()
{
	int ret=0;
	double Xopx[MAX_SPCNT], Xopy[MAX_SPCNT], Xopz[MAX_SPCNT], Eopx[MAX_SPCNT], Eopy[MAX_SPCNT], Eopz[MAX_SPCNT];
	double Etpx0[MAX_SPCNT], Etpy0[MAX_SPCNT], Etpz0[MAX_SPCNT], Etpx1[MAX_SPCNT], Etpy1[MAX_SPCNT], Etpz1[MAX_SPCNT];
	maxrlen = 20; //仮

	//	double alpha[MAX_SPCNT], cosa[MAX_SPCNT], sina[MAX_SPCNT], tana[MAX_SPCNT], da[MAX_SPCNT];
	//	double betal[MAX_SPCNT], cotbl[MAX_SPCNT], cosbl[MAX_SPCNT], sinbl[MAX_SPCNT];
	//	double rlx[MAX_SPCNT], rly[MAX_SPCNT], rlz[MAX_SPCNT], rlx_cp[MAX_SPCNT], rly_cp[MAX_SPCNT];
	//	int rlflg[MAX_SPCNT]; // -1:未設定, 0:曲線, 1:紙端上, 2:紙端右, 3:紙端下, 4:紙端左
	//	double rllen[MAX_SPCNT]; // crease pattern をもとに長さ求める
	//	double maxrlen;

	// X[i-1], X[i], X[i+1] -> osculating planes on X[i]
	for( int i=1; i<Xcnt-1; i++ ){
		double x0,y0,z0,x1,y1,z1;
		x0 = Xx[i-1]-Xx[i];		y0 = Xy[i-1]-Xy[i];		z0 = Xz[i-1]-Xz[i];
		x1 = Xx[i+1]-Xx[i];		y1 = Xy[i+1]-Xy[i];		z1 = Xz[i+1]-Xz[i];
		Xopx[i] = y0*z1-y1*z0;	Xopy[i] = z0*x1-z1*x0;	Xopz[i] = x0*y1-x1*y0;
		normalize_v3( &(Xopx[i]), &(Xopy[i]), &(Xopz[i]) );
	}

	// osculating planes on edge (X[i], X[i+1])
	Eopx[0] = Xopx[1];
	Eopy[0] = Xopy[1];
	Eopz[0] = Xopz[1];
	for( int i=1; i<Xcnt-2; i++ ){
		//double x0,y0,z0,x1,y1,z1;
		Eopx[i] = Xopx[i]+Xopx[i+1];
		Eopy[i] = Xopy[i]+Xopy[i+1];
		Eopz[i] = Xopz[i]+Xopz[i+1];
		normalize_v3( &(Eopx[i]), &(Eopy[i]), &(Eopz[i]) );
	}
	Eopx[Xcnt-2] = Xopx[Xcnt-2];
	Eopy[Xcnt-2] = Xopy[Xcnt-2];
	Eopz[Xcnt-2] = Xopz[Xcnt-2];

	// tangent planes on edge (X[i], X[i+1]) -> other side
	for( int i=0; i<Xcnt-1; i++ ){
		Etpx0[i] = Ty[i]*rrz[i]-rry[i]*Tz[i];
		Etpy0[i] = Tz[i]*rrx[i]-rrz[i]*Tx[i];
		Etpz0[i] = Tx[i]*rry[i]-rrx[i]*Ty[i];
		normalize_v3( &(Etpx0[i]), &(Etpy0[i]), &(Etpz0[i]) );
		double ip = Etpx0[i]*Eopx[i] + Etpy0[i]*Eopy[i] + Etpz0[i]*Eopz[i];
		Etpx1[i] = Etpx0[i] - 2.0*(Etpx0[i]-ip*Eopx[i]);
		Etpy1[i] = Etpy0[i] - 2.0*(Etpy0[i]-ip*Eopy[i]);
		Etpz1[i] = Etpz0[i] - 2.0*(Etpz0[i]-ip*Eopz[i]);
		normalize_v3( &(Etpx1[i]), &(Etpy1[i]), &(Etpz1[i]) );
	}

	// intersecting lines of two tangent planes = outer product of two normal vectors
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		rlx[i] = Etpy1[i+1]*Etpz1[i]-Etpy1[i]*Etpz1[i+1];
		rly[i] = Etpz1[i+1]*Etpx1[i]-Etpz1[i]*Etpx1[i+1];
		rlz[i] = Etpx1[i+1]*Etpy1[i]-Etpx1[i]*Etpy1[i+1];
		normalize_v3( &(rlx[i]), &(rly[i]), &(rlz[i]) );
	}

	// 合計180°になるように cone に射影


	// (rrx,rry,rrz), TNB → alpha, beta*, rl*_cp
	FILE *fp=fopen("output/calcleftK.csv","w");
	if(fp){
		fprintf(fp, "alpha[i],cosa[i],sina[i],tana[i],alpha1,cosa1,sina1,tana1,cosbr[i],cosbr1,cosbl[i],sinbl[i],cosbl1,sinbl1\n");
	}
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		double ip, bnx0,bny0,bnz0, bnx1,bny1,bnz1;
		double cosa1,sina1,tana1,alpha1, cosbr1, cosbl1,sinbl1;
		// BN面に射影
		ip = rrx[i]*Tx[i] + rry[i]*Ty[i] + rrz[i]*Tz[i];
		bnx0 = rrx[i] - ip*Tx[i];
		bny0 = rry[i] - ip*Ty[i];
		bnz0 = rrz[i] - ip*Tz[i];
		normalize_v3( &bnx0, &bny0, &bnz0 );
		cosa[i] = bnx0*(-Nx[i]) + bny0*(-Ny[i]) + bnz0*(-Nz[i]); // r なので -N
		sina[i] = bnx0*Bx[i] + bny0*By[i] + bnz0*Bz[i];
		tana[i] = sina[i]/cosa[i];
		alpha[i] = atan(tana[i]);

		ip = rlx[i]*Tx[i] + rly[i]*Ty[i] + rlz[i]*Tz[i];
		bnx1 = rlx[i] - ip*Tx[i];
		bny1 = rly[i] - ip*Ty[i];
		bnz1 = rlz[i] - ip*Tz[i];
		normalize_v3( &bnx1, &bny1, &bnz1 );
		cosa1 = bnx1*Nx[i] + bny1*Ny[i] + bnz1*Nz[i];
		sina1 = bnx1*Bx[i] + bny1*By[i] + bnz1*Bz[i];
		tana1 = sina1/cosa1;
		alpha1 = atan(tana1);

		cosbr[i] = rrx_cp[i]*Tx2d[i] + rry_cp[i]*Ty2d[i];
		sinbr[i] = rrx_cp[i]*(-Nx2d[i]) + rry_cp[i]*(-Ny2d[i]); // r なので -N
		cotbr[i] = cosbr[i]/sinbr[i];
		betar[i] = acos( cosbr[i] );
		cosbr1 = rrx[i]*Tx[i] + rry[i]*Ty[i];

		cosbl1 = rlx[i]*Tx[i] + rly[i]*Ty[i];	// <- 違うかも
		sinbl1 = sqrt(1.0-cosbl1*cosbl1);
		rlx_cp[i] = cosbl1*Tx2d[i] -sinbl1*Ty2d[i];
		rly_cp[i] = sinbl1*Tx2d[i] +cosbl1*Ty2d[i];

		cosbl[i] = rlx_cp[i]*Tx2d[i] + rly_cp[i]*Ty2d[i];
		sinbl[i] = rlx_cp[i]*Nx2d[i] + rly_cp[i]*Ny2d[i];
		cotbl[i] = cosbl[i]/sinbl[i];
		betal[i] = acos( cosbl[i] );

		rlflg[i] = RETYPE_UNDEF;
		rllen[i] = maxrlen;

		if(fp){
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				alpha[i], cosa[i], sina[i], tana[i], alpha1, cosa1, sina1, tana1,
				cosbr[i], cosbr1, cosbl[i], sinbl[i], cosbl1, sinbl1);
		}
	}
	if(fp){ fclose(fp); fp=NULL; }

	return ret;
}

int crease::calcRightK()
{
	int ret=0;
	double Xopx[MAX_SPCNT], Xopy[MAX_SPCNT], Xopz[MAX_SPCNT], Eopx[MAX_SPCNT], Eopy[MAX_SPCNT], Eopz[MAX_SPCNT];
	double Etpx0[MAX_SPCNT], Etpy0[MAX_SPCNT], Etpz0[MAX_SPCNT], Etpx1[MAX_SPCNT], Etpy1[MAX_SPCNT], Etpz1[MAX_SPCNT];
	maxrlen = 20; //仮

	//	double alpha[MAX_SPCNT], cosa[MAX_SPCNT], sina[MAX_SPCNT], tana[MAX_SPCNT], da[MAX_SPCNT];
	//	double betal[MAX_SPCNT], cotbl[MAX_SPCNT], cosbl[MAX_SPCNT], sinbl[MAX_SPCNT];
	//	double rlx[MAX_SPCNT], rly[MAX_SPCNT], rlz[MAX_SPCNT], rlx_cp[MAX_SPCNT], rly_cp[MAX_SPCNT];
	//	int rlflg[MAX_SPCNT]; // -1:未設定, 0:曲線, 1:紙端上, 2:紙端右, 3:紙端下, 4:紙端左
	//	double rllen[MAX_SPCNT]; // crease pattern をもとに長さ求める
	//	double maxrlen;

	// X[i-1], X[i], X[i+1] -> osculating planes on X[i]
	for( int i=1; i<Xcnt-1; i++ ){
		double x0,y0,z0,x1,y1,z1;
		x0 = Xx[i-1]-Xx[i];		y0 = Xy[i-1]-Xy[i];		z0 = Xz[i-1]-Xz[i];
		x1 = Xx[i+1]-Xx[i];		y1 = Xy[i+1]-Xy[i];		z1 = Xz[i+1]-Xz[i];
		Xopx[i] = y0*z1-y1*z0;	Xopy[i] = z0*x1-z1*x0;	Xopz[i] = x0*y1-x1*y0;
		normalize_v3( &(Xopx[i]), &(Xopy[i]), &(Xopz[i]) );
	}

	// osculating planes on edge (X[i], X[i+1])
	Eopx[0] = Xopx[1];
	Eopy[0] = Xopy[1];
	Eopz[0] = Xopz[1];
	for( int i=1; i<Xcnt-2; i++ ){
		Eopx[i] = Xopx[i]+Xopx[i+1];
		Eopy[i] = Xopy[i]+Xopy[i+1];
		Eopz[i] = Xopz[i]+Xopz[i+1];
		normalize_v3( &(Eopx[i]), &(Eopy[i]), &(Eopz[i]) );
	}
	Eopx[Xcnt-2] = Xopx[Xcnt-2];
	Eopy[Xcnt-2] = Xopy[Xcnt-2];
	Eopz[Xcnt-2] = Xopz[Xcnt-2];

	// tangent planes on edge (X[i], X[i+1]) -> other side
	for( int i=0; i<Xcnt-1; i++ ){
		Etpx0[i] = Ty[i]*rlz[i]-rly[i]*Tz[i];
		Etpy0[i] = Tz[i]*rlx[i]-rlz[i]*Tx[i];
		Etpz0[i] = Tx[i]*rly[i]-rlx[i]*Ty[i];
		normalize_v3( &(Etpx0[i]), &(Etpy0[i]), &(Etpz0[i]) );
		double ip = Etpx0[i]*Eopx[i] + Etpy0[i]*Eopy[i] + Etpz0[i]*Eopz[i];
		Etpx1[i] = Etpx0[i] - 2.0*(Etpx0[i]-ip*Eopx[i]);
		Etpy1[i] = Etpy0[i] - 2.0*(Etpy0[i]-ip*Eopy[i]);
		Etpz1[i] = Etpz0[i] - 2.0*(Etpz0[i]-ip*Eopz[i]);
		normalize_v3( &(Etpx1[i]), &(Etpy1[i]), &(Etpz1[i]) );
	}

	// intersecting lines of two tangent planes = outer product of two normal vectors
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		rrx[i] = Etpy1[i+1]*Etpz1[i]-Etpy1[i]*Etpz1[i+1];
		rry[i] = Etpz1[i+1]*Etpx1[i]-Etpz1[i]*Etpx1[i+1];
		rrz[i] = Etpx1[i+1]*Etpy1[i]-Etpx1[i]*Etpy1[i+1];
		normalize_v3( &(rrx[i]), &(rry[i]), &(rrz[i]) );
	}

	// 合計180°になるように cone に射影


	// (rrx,rry,rrz), TNB → alpha, beta*, rl*_cp
	FILE *fp=fopen("output/calcrightK.csv","w");
	if(fp){
		fprintf(fp, "alpha[i],cosa[i],sina[i],tana[i],alpha1,cosa1,sina1,tana1,cosbl[i],cosbl1,cosbr[i],sinbr[i],cosbr1,sinbr1\n");
	}
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		double ip, bnx0,bny0,bnz0, bnx1,bny1,bnz1;
		double cosa1,sina1,tana1,alpha1, cosbr1,sinbr1, cosbl1;
		// BN面に射影
		ip = rlx[i]*Tx[i] + rly[i]*Ty[i] + rlz[i]*Tz[i];
		bnx0 = rlx[i] - ip*Tx[i];
		bny0 = rly[i] - ip*Ty[i];
		bnz0 = rlz[i] - ip*Tz[i];
		normalize_v3( &bnx0, &bny0, &bnz0 );
		cosa[i] = bnx0*Nx[i] + bny0*Ny[i] + bnz0*Nz[i];
		sina[i] = bnx0*Bx[i] + bny0*By[i] + bnz0*Bz[i];
		tana[i] = sina[i]/cosa[i];
		alpha[i] = atan(tana[i]);

		ip = rrx[i]*Tx[i] + rry[i]*Ty[i] + rrz[i]*Tz[i];
		bnx1 = rrx[i] - ip*Tx[i];
		bny1 = rry[i] - ip*Ty[i];
		bnz1 = rrz[i] - ip*Tz[i];
		normalize_v3( &bnx1, &bny1, &bnz1 );
		cosa1 = bnx1*(-Nx[i]) + bny1*(-Ny[i]) + bnz1*(-Nz[i]); // r なので -N
		sina1 = bnx1*Bx[i] + bny1*By[i] + bnz1*Bz[i];
		tana1 = sina1/cosa1;
		alpha1 = atan(tana1);

		cosbl[i] = rlx_cp[i]*Tx2d[i] + rly_cp[i]*Ty2d[i];
		sinbl[i] = rlx_cp[i]*Nx2d[i] + rly_cp[i]*Ny2d[i];
		cotbl[i] = cosbl[i]/sinbl[i];
		betal[i] = acos( cosbl[i] );
		cosbl1 = rlx[i]*Tx[i] + rly[i]*Ty[i];

		cosbr1 = rrx[i]*Tx[i] + rry[i]*Ty[i];	// <- 違うかも
		sinbr1 = sqrt(1.0-cosbr1*cosbr1);
		rrx_cp[i] = cosbr1*Tx2d[i] +sinbr1*Ty2d[i];
		rry_cp[i] = -sinbr1*Tx2d[i] + cosbr1*Ty2d[i];

		cosbr[i] = rrx_cp[i]*Tx2d[i] + rry_cp[i]*Ty2d[i];
		sinbr[i] = rrx_cp[i]*(-Nx2d[i]) + rry_cp[i]*(-Ny2d[i]);
		cotbr[i] = cosbr[i]/sinbr[i];
		betar[i] = acos( cosbr[i] );

		rrflg[i] = RETYPE_UNDEF;
		rllen[i] = maxrlen;

		if(fp){
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				alpha[i], cosa[i], sina[i], tana[i], alpha1, cosa1, sina1, tana1,
				cosbl[i], cosbl1, cosbr[i], sinbr[i], cosbr1, sinbr1);
		}
	}
	if(fp){ fclose(fp); fp=NULL; }

	return ret;
}

