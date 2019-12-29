#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "crease.h"
#include "util.h"

// P*をファイル読み込み
int crease::load( char *fname )
{
	int i, ret=0;
	char buf[512], *pbuf=NULL, retbuf[256];
	FILE *fp = fopen( fname ,"r" );
	if( fp==NULL ){ ret=-1; goto end; }
	fgets(buf,1024,fp);
	if( buf[0] == '1' && buf[1] != '0' ){
		fgets(buf,1024,fp);	pbuf = buf;	csvread(&pbuf, retbuf);	Pcnt = atoi(retbuf);
		fgets(buf,1024,fp);
		for( i=0; i<Pcnt; i++ ){
			fgets(buf,1024,fp);
			pbuf = buf;
			csvread(&pbuf, retbuf);	Px2d[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Py2d[i] = atof(retbuf);
		}
		fgets(buf,1024,fp);
		for( i=0; i<Pcnt; i++ ){
			fgets(buf,1024,fp);
			pbuf = buf;
			csvread(&pbuf, retbuf);	Px[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Py[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Pz[i] = atof(retbuf);
		}
		ret = 1;
	} else if( buf[0] == '2' ){
		// int Xcnt
		// double kv[MAX_SPCNT], tr[MAX_SPCNT], alpha[MAX_SPCNT];
		fgets(buf,1024,fp);	pbuf = buf;	csvread(&pbuf, retbuf);	Pcnt = atoi(retbuf);
		for( i=0; i<Pcnt; i++ ){
			double tmpk, tmpt, tmpa, tmpk2;
			fgets(buf,1024,fp);
			sscanf( buf, "%lf %lf %lf %lf", &tmpk, &tmpt, &tmpa, &tmpk2 );
			Px[i] = tmpk;
			Py[i] = tmpt;
			Pa[i] = tmpa*M_PI/180.0;
			Px2d[i] = tmpk2;
		}
		ret = 2;
	} else if( buf[0] == '3' ){
		fgets(buf,1024,fp);	pbuf = buf;	csvread(&pbuf, retbuf);	Pcnt = atoi(retbuf);
		fgets(buf,1024,fp);
		for( i=0; i<Pcnt; i++ ){
			fgets(buf,1024,fp);
			pbuf = buf;
			csvread(&pbuf, retbuf);	Px2d[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Py2d[i] = atof(retbuf);
		}
		ret = 3;
	} else if( buf[0] == '4' ){
		fgets(buf,1024,fp);	pbuf = buf;	csvread(&pbuf, retbuf);	Pcnt = atoi(retbuf);
		fgets(buf,1024,fp);
		for( i=0; i<Pcnt; i++ ){
			fgets(buf,1024,fp);
			pbuf = buf;
			csvread(&pbuf, retbuf);	Px[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Py[i] = atof(retbuf);
			csvread(&pbuf, retbuf);	Pz[i] = atof(retbuf);
		}
		ret = 4;
	} else if( buf[0] == '1' && buf[1] == '0' ){ // 10
		int mode;
		sscanf( buf, "%d %d", &ret, &mode );
		fgets(buf,1024,fp);	pbuf = buf;	csvread(&pbuf, retbuf);	Pcnt = atoi(retbuf);
		for( i=0; i<Pcnt; i++ ){
			double tmpk, tmpt, tmpa, tmpk2;
			fgets(buf,1024,fp);
			sscanf( buf, "%lf %lf %lf %lf", &tmpk, &tmpt, &tmpa, &tmpk2 );
			Px[i] = tmpk;
			Py[i] = tmpt;
			Pa[i] = tmpa;
			Px2d[i] = tmpk2;
		}
		ret += mode;
	}
end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}

int crease::dumpP( char *fname, int mode )
{
	int i, ret=0;

	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "10 %d # file type, mode\n", mode );
		fprintf( fp, "%d # plot count\n", Pcnt );
		for( i=0; i<Pcnt; i++ ){
			fprintf( fp, "%f\t%f\t%f\t%f\n", Px[i], Py[i], Pa[i], Px2d[i] );
		}
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

int crease::dumpP_( char *fname, int mode )
{
	int i, ret=0, idxk[MAX_CPCNT], idxt[MAX_CPCNT], idxa[MAX_CPCNT];

	for( i=0; i<Pcnt; i++ ){
		idxk[i] = idxt[i] = idxa[i] = (int)((double)i*(double)Xcnt/(double)(Pcnt-1) + 0.5);
	}
	if( idxk[0] < 2 ){ idxk[0] = 2; }
	if( idxt[0] < 3 ){ idxt[0] = 3; }
	if( idxa[0] < 0 ){ idxa[0] = 0; }
	if( idxk[Pcnt-1] > Xcnt-3 ){ idxk[Pcnt-1] = Xcnt-3; }
	if( idxt[Pcnt-1] > Xcnt-4 ){ idxt[Pcnt-1] = Xcnt-4; }
	if( idxa[Pcnt-1] > Xcnt-1 ){ idxa[Pcnt-1] = Xcnt-1; }

	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "2	%d # file type, mode\n", mode );
		fprintf( fp, "%d # plot count\n", Pcnt );
		for( i=0; i<Pcnt; i++ ){
			fprintf( fp, "%f\t%f\t%d\t%f\n",
				kv[idxk[i]], tr[idxt[i]], (int)(alpha[idxa[i]]*180.0/M_PI+0.5), kv[idxk[i]] );
		}
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

int crease::dumpX( char *fname )
{
	int i, ret=0;
	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "2	# file type, mode\n" );
		fprintf( fp, "%d # plot count\n", Xcnt );
		fprintf( fp, "Xx2d,Xy2d,d2d,k2d,Xx,Xy,Xz,dx,kv,tr,alpha\n" );
		for( i=0; i<Xcnt; i++ ){
			fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
				Xx2d[i], Xy2d[i], d2d[i], k2d[i], Xx[i], Xy[i], Xz[i], dx[i], kv[i], tr[i], alpha[i] );
		}
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

// 頂点周りの角度
int crease::check180( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");

	if( !fp ){
		ret = -1; goto end;
	}

	ret = check180( fp );

end:
	if( fp ){ fclose(fp); }
	return ret;
}

int crease::check180( FILE *fp )
{
	int ret=0;

	if( !fp ){
		ret = -1; goto end;
	}
	fprintf( fp, "i, ip0r, ip1r, ip0l, ip1l, a0r, a1r, a0l, a1l, atotal-2pi\n");
	for( int i=Xsidx; i<=Xeidx; i++ ){
		if( rllen[i]==0 || rrlen[i]==0.0 ){
			fprintf( fp, "%d,0,0,0,0,0,0,0,0,0\n", i);
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
		fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			i, ip0r, ip1r, ip0l, ip1l, a0r, a1r, a0l, a1l, a0r+a1r+a0l+a1l-2*M_PI);
	}
end:
	return ret;
}

// quad の平面度
int crease::checkquadplane( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");

	if( !fp ){
		ret = -1; goto end;
	}

	ret = checkquadplane( fp );

end:
	if( fp ){ fclose(fp); }
	return ret;
}

int crease::checkquadplane( FILE *fp )
{
	int ret=0;

	if( !fp ){
		ret = -1; goto end;
	}
	fprintf( fp, "i, dx0,dy0,dz0, opxL,opyL,opzL, dxL,dyL,dzL, ipL, opxR,opyR,opzR, dxR,dyR,dzR, ipR\n" );
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
		fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			i, dx0,dy0,dz0, opxL,opyL,opzL, dxL,dyL,dzL, ipL, opxR,opyR,opzR, dxR,dyR,dzR, ipR );
	}
end:
	return ret;
}

int crease::loadm2m3( char *fname )
{
	int ret=0;
	char buf[1024];

	// position & rotation (3D, CP)
	//double m3[16], m2[9];

	FILE *fp=fopen( fname, "r" );
	if( !fp ){
		ret = -1;
	} else {
		fgets( buf, 1024, fp );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf", &(m2[0]), &(m2[1]), &(m2[2]) );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf", &(m2[3]), &(m2[4]), &(m2[5]) );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf", &(m2[6]), &(m2[7]), &(m2[8]) );
		fgets( buf, 1024, fp );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf %lf", &(m3[0]), &(m3[1]), &(m3[2]), &(m3[3]) );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf %lf", &(m3[4]), &(m3[5]), &(m3[6]), &(m3[7]) );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf %lf", &(m3[8]), &(m3[9]), &(m3[10]), &(m3[11]) );
		fgets( buf, 1024, fp );	sscanf( buf, "%lf %lf %lf %lf", &(m3[12]), &(m3[13]), &(m3[14]), &(m3[15]) );
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

int crease::dumpm2m3( char *fname )
{
	int ret=0;

	FILE *fp=fopen( fname, "w" );
	if( !fp ){
		ret = -1;
	} else {		fprintf( fp, "m2\n" );
	fprintf( fp, "%f\t%f\t%f\n", m2[0], m2[1], m2[2] );
	fprintf( fp, "%f\t%f\t%f\n", m2[3], m2[4], m2[5] );
	fprintf( fp, "%f\t%f\t%f\n", m2[6], m2[7], m2[8] );
	fprintf( fp, "m3\n" );
	fprintf( fp, "%f\t%f\t%f\t%f\n", m3[0], m3[1], m3[2], m3[3] );
	fprintf( fp, "%f\t%f\t%f\t%f\n", m3[4], m3[5], m3[6], m3[7] );
	fprintf( fp, "%f\t%f\t%f\t%f\n", m3[8], m3[9], m3[10], m3[11] );
	fprintf( fp, "%f\t%f\t%f\t%f\n", m3[12], m3[13], m3[14], m3[15] );
	fclose(fp); fp=NULL;
	}

end:
	return ret;
}

int crease::dumpAll( char *fname )
{
	int i,j,ret=0;
	FILE *fp = fopen( fname, "w" );
	if( fp==NULL ){
		ret = -1;
		goto end;
	}

	//	int fileioflg[30];
	//	enum{ PX2, XX2, TX2, NX2, D2, K2, PX3, XX3, TX3, NX3, BX3, D3, K3, T3, PA, ALPHA, ALPHA0, DA, PB, BETA, BETA0, R3, R2 };

	for( i=0; i<30; i++ ){
		if( fileioflg[i]==0 ){ 
			continue;
		}
		switch( i ){
			case PX2: fprintf( fp, "Px2d,Py2d,"); break;
			case XX2: fprintf( fp, "Xx2d,Xy2d,"); break;
			case TX2: fprintf( fp, "Tx2d,Ty2d,"); break;
			case NX2: fprintf( fp, "Nx2d,Ny2d,"); break;
			case D2: fprintf( fp, "d2d,"); break;
			case K2: fprintf( fp, "k2d,"); break;
			case PX3: fprintf( fp, "Px,Py,Pz,"); break;
			case XX3: fprintf( fp, "Xx,Xy,Xz,"); break;
			case TX3: fprintf( fp, "Tx,Ty,Tz,"); break;
			case NX3: fprintf( fp, "Nx,Ny,Nz,"); break;
			case BX3: fprintf( fp, "Bx,By,Bz,"); break;
			case D3: fprintf( fp, "dx,"); break;
			case K3: fprintf( fp, "kv,"); break;
			case T3: fprintf( fp, "tr,"); break;
			case PA: fprintf( fp, "Pz,"); break;
			case ALPHA: fprintf( fp, "alpha,"); break;
			case ALPHA0: fprintf( fp, "cosa,sina,tana,da,"); break;
			case DA: fprintf( fp, "da,"); break;
			case PB: fprintf( fp, "Pbl,Pbr,"); break;
			case BETA: fprintf( fp, "betal, betar,"); break;
			case BETA0: fprintf( fp, "cotbl,cotbr,cosbl,sinbl,cosbr,sinbr,"); break;
			case R3: fprintf( fp, "rlx,rly,rlz,rrx,rry,rrz,"); break;
			case R2: fprintf( fp, "rlx_cp,rly_cp,rrx_cp,rry_cp,"); break;
			default: break;
		}
	}
	fprintf( fp, "\n");
#if 1
	for( j=0; j<Xcnt; j++ ){
		for( i=0; i<30; i++ ){
			if( fileioflg[i]==0 ){ 
				continue;
			}
			switch( i ){
				case PX2:
					if( j<Pcnt ){
						fprintf( fp, "%f,%f,", Px2d[j], Py2d[j] );
					} else {
						fprintf( fp, ",,");
					}
					break;
				case XX2: fprintf( fp, "%f,%f,", Xx2d[j], Xy2d[j] ); break;
				case TX2: fprintf( fp, "%f,%f,", Tx2d[j], Ty2d[j] ); break;
				case NX2: fprintf( fp, "%f,%f,", Nx2d[j], Ny2d[j] ); break;
				case D2: fprintf( fp, "%f,", d2d[j] ); break;
				case K2: fprintf( fp, "%f,", k2d[j] ); break;
				case PX3:
					if( j<Pcnt ){
						fprintf( fp, "%f,%f,%f,", Px[j], Py[j], Px[j] );
					} else {
						fprintf( fp, ",,,");
					}
					break;
				case XX3: fprintf( fp, "%f,%f,%f,", Xx[j], Xy[j], Xz[j] ); break;
				case TX3: fprintf( fp, "%f,%f,%f,", Tx[j], Ty[j], Tz[j] ); break;
				case NX3: fprintf( fp, "%f,%f,%f,", Nx[j], Ny[j], Nz[j] ); break;
				case BX3: fprintf( fp, "%f,%f,%f,", Bx[j], By[j], Bz[j] ); break;
				case D3: fprintf( fp, "%f,", dx[j] ); break;
				case K3: fprintf( fp, "%f,", kv[j] ); break;
				case T3: fprintf( fp, "%f,", tr[j] ); break;
				case PA:
					if( j<Pcnt ){
						fprintf( fp, "%f,", Pa[j] );
					} else {
						fprintf( fp, ",");
					}
					break;
				case ALPHA: fprintf( fp, "%f,", alpha[j] ); break;
				case ALPHA0: fprintf( fp, "%f,%f,%f,%f,", cosa[j], sina[j], tana[j], da[j] ); break;
				case DA: fprintf( fp, "%f,", da[j] ); break;
				case PB: break;
				case BETA: fprintf( fp, "%f,%f,", betal[j], betar[j] ); break;
				case BETA0: fprintf( fp, "%f,%f,%f,%f,%f,%f,", cotbl[j], cotbr[j], cosbl[j], sinbl[j], cosbr[j], sinbr[j] ); break;
				case R3: fprintf( fp, "%f,%f,%f,%f,%f,%f,", rlx[j], rly[j], rlz[j], rrx[j], rry[j], rrz[j] ); break;
				case R2: fprintf( fp, "%f,%f,%f,%f,", rlx_cp[j], rly_cp[j], rrx_cp[j], rry_cp[j] ); break;
				default: break;
			}
		}
		fprintf( fp, "\n");
	}
#endif
end:
	if( fp ){ fclose( fp ); fp=NULL; }
	return ret;
}

#if 1
int crease::checkCP3D( char *fname )
{
	int i,j, ret=0;
	// 角度誤差
	double aerr_l0[MAX_SPCNT], aerr_l1[MAX_SPCNT], aerr_r0[MAX_SPCNT], aerr_r1[MAX_SPCNT], aerr_t[MAX_SPCNT];
	// 誤差平均
	double ddiff;
	double averr_l0, averr_l1, averr_r0, averr_r1, averr_t;
	double pkdiff, ptdiff, padiff, pk2diff;
	FILE *fp=NULL;

	if(fname){
		fp = fopen( fname, "w" );
	}

	//d誤差, 角度誤差×4, 合計角度, 代表点誤差
	int acnt=0;
	double dx,dy,dz, tx0,ty0,tz0, tx1,ty1,tz1, tx2d0,ty2d0, tx2d1,ty2d1, ip;

	averr_l0=0.0; averr_l1=0.0; averr_r0=0.0; averr_r1=0.0; averr_t=0.0;
	ddiff=0.0; pkdiff=0.0; ptdiff=0.0; padiff=0.0; pk2diff=0.0;
	memset( aerr_l0, 0, sizeof(double)*MAX_SPCNT );
	memset( aerr_l1, 0, sizeof(double)*MAX_SPCNT );
	memset( aerr_r0, 0, sizeof(double)*MAX_SPCNT );
	memset( aerr_r1, 0, sizeof(double)*MAX_SPCNT );
	memset( aerr_t, 0, sizeof(double)*MAX_SPCNT );

	if(fp){
		fprintf( fp, "i,,d2d,d3d,,a2d_tl0,a2d_tr0,a2d_tl1,a2d_tr1,a2d_ttl,,a3d_tl0,a3d_tr0,a3d_tl1,a3d_tr1,a3d_ttl,,df_tl0,df_tr0,df_tl1,df_tr1,df_ttl\n" );
	}
	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){

		if( i==CEMGN ){
			dx = Xx2d[i]-Xx2d[i-1];
			dy = Xy2d[i]-Xy2d[i-1];
			double d2d = sqrt( dx*dx + dy*dy );
			tx2d0 = dx/d2d;
			ty2d0 = dy/d2d;
			dx = Xx[i]-Xx[i-1];
			dy = Xy[i]-Xy[i-1];
			dz = Xz[i]-Xz[i-1];
			double d3d = sqrt( dx*dx + dy*dy + dz*dz );
			tx0 = dx/d3d;
			ty0 = dy/d3d;
			tz0 = dz/d3d;
		} else {
			tx2d0 = tx2d1;
			ty2d0 = ty2d1;
			tx0 = tx1;
			ty0 = ty1;
			tz0 = tz1;
		}

		// check length
		dx = Xx2d[i+1]-Xx2d[i];
		dy = Xy2d[i+1]-Xy2d[i];
		double d2d = sqrt( dx*dx + dy*dy );
		tx2d1 = dx/d2d;
		ty2d1 = dy/d2d;
		dx = Xx[i+1]-Xx[i];
		dy = Xy[i+1]-Xy[i];
		dz = Xz[i+1]-Xz[i];
		double d3d = sqrt( dx*dx + dy*dy + dz*dz );
		tx1 = dx/d3d;
		ty1 = dy/d3d;
		tz1 = dz/d3d;

		// check ruling angle
		ip = -tx2d0*rlx_cp[i] -ty2d0*rly_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a2d_tl0 = acos( ip );
		ip = -tx2d0*rrx_cp[i] -ty2d0*rry_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a2d_tr0 = acos( ip );
		ip = tx2d1*rlx_cp[i] +ty2d1*rly_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a2d_tl1 = acos( ip );
		ip = tx2d1*rrx_cp[i] +ty2d1*rry_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a2d_tr1 = acos( ip );
		double a2d_ttl = a2d_tl0 + a2d_tr0 + a2d_tl1 + a2d_tr1;

		ip = -tx0*rlx[i] -ty0*rly[i] -tz0*rlz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a3d_tl0 = acos( ip );
		ip = -tx0*rrx[i] -ty0*rry[i] -tz0*rrz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a3d_tr0 = acos( ip );
		ip = tx1*rlx[i] +ty1*rly[i] +tz1*rlz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a3d_tl1 = acos( ip );
		ip = tx1*rrx[i] +ty1*rry[i] +tz1*rrz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double a3d_tr1 = acos( ip );
		double a3d_ttl = a3d_tl0 + a3d_tr0 + a3d_tl1 + a3d_tr1;

		aerr_l0[i] = a2d_tl0-a3d_tl0;
		aerr_l1[i] = a2d_tl1-a3d_tl1;
		aerr_r0[i] = a2d_tr0-a3d_tr0;
		aerr_r1[i] = a2d_tr1-a3d_tr1;
		aerr_t[i] = a2d_ttl-a3d_ttl;

		if(fp){
			fprintf( fp, "%d,,%f,%f,,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f,,%f,%f,%f,%f,%f\n",
				i, d2d, d3d, a2d_tl0, a2d_tr0, a2d_tl1, a2d_tr1, a2d_ttl,
				a3d_tl0, a3d_tr0, a3d_tl1, a3d_tr1, a3d_ttl,
				aerr_l0[i], aerr_l1[i], aerr_r0[i], aerr_r1[i], aerr_t[i]);
		}
		// 平均値
		ddiff += fabs( d2d - d3d );
		averr_l0 += fabs( aerr_l0[i] );
		averr_l1 += fabs( aerr_l1[i] );
		averr_r0 += fabs( aerr_r0[i] );
		averr_r1 += fabs( aerr_r1[i] );
		averr_t += fabs( aerr_t[i] );
		acnt++;
	}
	ddiff  /= (double)acnt;
	averr_l0 /= (double)acnt;
	averr_l1 /= (double)acnt;
	averr_r0 /= (double)acnt;
	averr_r1 /= (double)acnt;
	averr_t /= (double)acnt;

#if 0 // ポリゴンごとに角度比較
	if( fp ){ fprintf( fp, "\ni,,ll2,rr2,,ll3,rr3,,dll,drr\n"); }

	for( i=CEMGN+1; i<Xcnt-CEMGN; i++ ){

		ip = rlx_cp[i-1]*rlx_cp[i] + rly_cp[i-1]*rly_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double ll2 = acos( ip );
		ip = rrx_cp[i-1]*rrx_cp[i] + rry_cp[i-1]*rry_cp[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double rr2 = acos( ip );

		ip = rlx[i-1]*rlx[i] + rly[i-1]*rly[i] + rlz[i-1]*rlz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double ll3 = acos( ip );
		ip = rrx[i-1]*rrx[i] + rry[i-1]*rry[i] + rrz[i-1]*rrz[i];
		ip = ip < 1.0 ? (ip > -1.0 ? ip : -1.0) : 1.0;
		double rr3 = acos( ip );

		double dll = ll2-ll3;
		double drr = rr2-rr3;

		if(fp){
			fprintf( fp, "%d,,%f,%f,,%f,%f,,%f,%f\n", i, ll2,rr2,ll3,rr3,dll,drr);
		}
	}
#endif

	double tmpPx2d[MAX_CPCNT], tmpPy2d[MAX_CPCNT], tmpPx[MAX_CPCNT], tmpPy[MAX_CPCNT], tmpPz[MAX_CPCNT], tmpPa[MAX_CPCNT];
	double dffPx2d[MAX_CPCNT], dffPy2d[MAX_CPCNT], dffPx[MAX_CPCNT], dffPy[MAX_CPCNT], dffPz[MAX_CPCNT], dffPa[MAX_CPCNT];
	memcpy( tmpPx2d, Px2d, sizeof(double)*MAX_CPCNT);
	memcpy( tmpPy2d, Py2d, sizeof(double)*MAX_CPCNT);
	memcpy( tmpPx, Px, sizeof(double)*MAX_CPCNT);
	memcpy( tmpPy, Py, sizeof(double)*MAX_CPCNT);
	memcpy( tmpPz, Pz, sizeof(double)*MAX_CPCNT);
	memcpy( tmpPa, Pa, sizeof(double)*MAX_CPCNT);

	setP_k(-1); // kv -> Px[_Pcnt]
	setP_t(-1); // tr -> Py[_Pcnt]
	setP_a(-1);  // alpha -> Pa[_Pcnt]
	setP_k2(-1); // kt -> Px2d[_Pcnt]

	if(fp){
		fprintf( fp, "\ni,dffPx2d,dffPy2d,dffPx,dffPy,dffPz,dffPa\n" );
	}

	for( i=0; i<Pcnt; i++){
		dffPx2d[i] = tmpPx2d[i] - Px2d[i];
		dffPy2d[i] = tmpPy2d[i] - Py2d[i];
		dffPx[i] = tmpPx[i] - Px[i];
		dffPy[i] = tmpPy[i] - Py[i];
		dffPz[i] = tmpPz[i] - Pz[i];
		dffPa[i] = tmpPa[i] - Pa[i];
		if(fp){
			fprintf( fp, "%d,%f,%f,%f,%f,%f,%f\n", i, dffPx2d[i], dffPy2d[i], dffPx[i], dffPy[i], dffPz[i], dffPa[i] );
		}
		// 平均値
		pkdiff += fabs( dffPx[i] );
		ptdiff += fabs( dffPy[i] );
		padiff += fabs( dffPa[i] );
		pk2diff+= fabs( dffPx2d[i] );
	}
	pkdiff /= (double)Pcnt;
	ptdiff /= (double)Pcnt;
	padiff /= (double)Pcnt;
	pk2diff/= (double)Pcnt;

	// 戻す
	memcpy( Px2d, tmpPx2d, sizeof(double)*MAX_CPCNT);
	memcpy( Py2d, tmpPy2d, sizeof(double)*MAX_CPCNT);
	memcpy( Px, tmpPx, sizeof(double)*MAX_CPCNT);
	memcpy( Py, tmpPy, sizeof(double)*MAX_CPCNT);
	memcpy( Pz, tmpPz, sizeof(double)*MAX_CPCNT);
	memcpy( Pa, tmpPa, sizeof(double)*MAX_CPCNT);

end:
	if(fp){ fclose(fp); fp=NULL; }
	return ret;
}
#endif