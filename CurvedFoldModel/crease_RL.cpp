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
			//Xx2d[i] = c->Xx2d[i] + c->rlx_cp[i] * c->rllen[i];
			//Xy2d[i] = c->Xy2d[i] + c->rly_cp[i] * c->rllen[i];
		}
	}
	// 曲線と右側rulingをコピー
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		if( c->rlflg[i] == RETYPE_CURVE ){
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
		} else {
			Xx[i] = Xx[i-1];	Xy[i] = Xy[i-1];	Xz[i] = Xz[i-1];
			Xx2d[i] = Xx2d[i-1];	Xy2d[i] = Xy2d[i-1];
			rrlen[i] = 0.0;
			rrflg[i] = RETYPE_UNDEF;
			rrx[i] = rry[i] = rrz[i] = rrx_cp[i] = rry_cp[i] = 0.0;
		}
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
			//Xx2d[i] = c->Xx2d[i] + c->rrx_cp[i] * c->rrlen[i];
			//Xy2d[i] = c->Xy2d[i] + c->rry_cp[i] * c->rrlen[i];
		}
	}
	// 曲線と右側rulingをコピー
	for( int i=Xsidx; i<Xeidx+1; i++ ){
		if( c->rrflg[i] == RETYPE_CURVE ){
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
		} else {
			Xx[i] = Xx[i-1];	Xy[i] = Xy[i-1];	Xz[i] = Xz[i-1];
			Xx2d[i] = Xx2d[i-1];	Xy2d[i] = Xy2d[i-1];
			rllen[i] = 0.0;
			rlflg[i] = RETYPE_UNDEF;
			rlx[i] = rly[i] = rlz[i] = rlx_cp[i] = rly_cp[i] = 0.0;
		}
	}
	return ret;
}

int crease::calcLeft()
{
	int ret=0;
	//maxrlen = 400; // 仮
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
	//maxrlen = 400; // 仮

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
