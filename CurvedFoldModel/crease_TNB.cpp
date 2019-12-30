#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"
#include "util.h"

// k2déZèo
//	double Xx2d[MAX_SPCNT], Xy2d[MAX_SPCNT];
//	double Tx2d[MAX_SPCNT], Ty2d[MAX_SPCNT];
//	double Nx2d[MAX_SPCNT], Ny2d[MAX_SPCNT], k2d[MAX_SPCNT];
int crease::calcTN2d()
{
	int ret=0;

	ret = calcTN2d( Xcnt, Xx2d, Xy2d, Tx2d, Ty2d, d2d, Nx2d, Ny2d, k2d );
#if 0
	{
		FILE *fp=fopen("output/dump_K2D.csv","w");
		if(fp){
			fprintf(fp, "Xx2d,Xy2d,,Tx2d,Ty2d,,Nx2d,Ny2d,,d2d,k2d\n");
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,,%f,%f,,%f,%f,,%f,%f\n", Xx2d[i], Xy2d[i], Tx2d[i], Ty2d[i], Nx2d[i], Ny2d[i], d2d[i], k2d[i]);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	return ret;
}

int crease::calcTN2d( int Xsidx0, int Xeidx0 )
{
	int ret=0;
	int size=Xeidx-Xsidx+1;

	memset( Tx2d, 0, sizeof(double)*Xcnt );
	memset( Ty2d, 0, sizeof(double)*Xcnt );
	memset( Nx2d, 0, sizeof(double)*Xcnt );
	memset( Ny2d, 0, sizeof(double)*Xcnt );
	memset(  d2d, 0, sizeof(double)*Xcnt );
	memset(  k2d, 0, sizeof(double)*Xcnt );

	ret = calcTN2d( size, &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
		&(Tx2d[Xsidx]), &(Ty2d[Xsidx]), &(d2d[Xsidx]),
		&(Nx2d[Xsidx]), &(Ny2d[Xsidx]), &(k2d[Xsidx]) );
	//Xsidx = Xsidx0;
	//Xeidx = Xeidx0;
#if 0
	{
		FILE *fp=fopen("output/dump_K2D.csv","w");
		if(fp){
			fprintf(fp, "Xx2d,Xy2d,,Tx2d,Ty2d,,Nx2d,Ny2d,,d2d,k2d\n");
			for( i=0; i<Xcnt; i++ ){
				fprintf(fp, "%f,%f,,%f,%f,,%f,%f,,%f,%f\n", Xx2d[i], Xy2d[i], Tx2d[i], Ty2d[i], Nx2d[i], Ny2d[i], d2d[i], k2d[i]);
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	return ret;
}

int crease::calcTN2d(int Xcnt0, double *Xx2d0, double *Xy2d0,
					 double *Tx2d0, double *Ty2d0, double *d2d0,
					 double *Nx2d0, double *Ny2d0, double *k2d0 )
{
	int i, ret=0;
	double d2d1[MAX_SPCNT], dx[MAX_SPCNT], dy[MAX_SPCNT];
	memset( Tx2d0, 0, sizeof(double)*Xcnt0 );
	memset( Ty2d0, 0, sizeof(double)*Xcnt0 );
	memset( Nx2d0, 0, sizeof(double)*Xcnt0 );
	memset( Ny2d0, 0, sizeof(double)*Xcnt0 );
	memset(  d2d0, 0, sizeof(double)*Xcnt0 );
	memset(  k2d0, 0, sizeof(double)*Xcnt0 );

	for( i=0; i<Xcnt0-1; i++ ){
		dx[i] = Xx2d0[i+1] - Xx2d0[i];
		dy[i] = Xy2d0[i+1] - Xy2d0[i];
		d2d1[i] = sqrt( dx[i]*dx[i] + dy[i]*dy[i] );
	}
	for( i=1; i<Xcnt0-1; i++ ){
		if( d2d1[i-1]<=0.0 && d2d1[i]<=0.0 ){
			Tx2d0[i] = Ty2d0[i] = d2d0[i] = 0.0;
			continue;
		}
		if( d2d1[i-1]>0.0 && d2d1[i]>0.0 ){
		Tx2d0[i] = (dx[i-1]/d2d1[i-1] + dx[i]/d2d1[i])*0.5;
		Ty2d0[i] = (dy[i-1]/d2d1[i-1] + dy[i]/d2d1[i])*0.5;
		d2d0[i] = (d2d1[i-1]+d2d1[i])*0.5;
		} else if( d2d1[i-1]>0.0 ) {
			Tx2d0[i] = dx[i-1];
			Ty2d0[i] = dy[i-1];
			d2d0[i] = d2d1[i-1];
		} else {
			Tx2d0[i] = dx[i];
			Ty2d0[i] = dy[i];
			d2d0[i] = d2d1[i];
		}
		normalize_v2( &Tx2d0[i], &Ty2d0[i] );
	}
	for( i=2; i<Xcnt0-2; i++ ){
		double d0 = d2d1[i];
		double d1 = d2d1[i-1];
		double nx0 = Tx2d0[i+1] - Tx2d0[i];
		double nx1 = Tx2d0[i]   - Tx2d0[i-1];
		double ny0 = Ty2d0[i+1] - Ty2d0[i];
		double ny1 = Ty2d0[i]   - Ty2d0[i-1];

		if( d0==0.0 && d1==0.0 ){
			Nx2d0[i] = Ny2d0[i] = k2d0[i] = 0.0;
			continue;
		}
		if( d0>0.0 && d1>0.0 ){
		Nx2d0[i] = (nx0/d0 + nx1/d1)*0.5;
		Ny2d0[i] = (ny0/d0 + ny1/d1)*0.5;
		} else if( d0>0.0 ){
			Nx2d0[i] = nx0;
			Ny2d0[i] = ny0;
		} else {
			Nx2d0[i] = nx1;
			Ny2d0[i] = ny1;
		}
		k2d0[i] = sqrt( Nx2d0[i]*Nx2d0[i] + Ny2d0[i]*Ny2d0[i] );
#if 1
		if( Tx2d0[i]*Ny2d0[i] - Ty2d0[i]*Nx2d0[i] < 0 ){
			k2d0[i] = -k2d0[i];
		}
#else
		if((i==2 && Nx2d0[i]*Nx2d0[0] + Ny2d0[i]*Ny2d0[0] < 0 )
			|| (i>2 && Nx2d0[i]*Nx2d0[i-1] + Ny2d0[i]*Ny2d0[i-1] < 0 )){ // îΩì]
				k2d0[i] = -k2d0[i];
		}
#endif
		if( fabs(k2d0[i]) > 0.0/*ZERO*/ ){
			Nx2d0[i] /= k2d0[i];
			Ny2d0[i] /= k2d0[i];
		}
	}
end:
	return ret;
}

// TNBéZèo
//	double Xx[MAX_SPCNT], Xy[MAX_SPCNT], Xz[MAX_SPCNT];
//	double Tx[MAX_SPCNT], Ty[MAX_SPCNT], Tz[MAX_SPCNT];
//	double Nx[MAX_SPCNT], Ny[MAX_SPCNT], Nz[MAX_SPCNT], kv[MAX_SPCNT];
//	double Bx[MAX_SPCNT], By[MAX_SPCNT], Bz[MAX_SPCNT], tr[MAX_SPCNT];
int crease::calcTNB()
{
	int ret=0;

	calcTNB( Xcnt, Xx, Xy, Xz, Tx, Ty, Tz, dx, Nx, Ny, Nz, kv, Bx, By, Bz, tr );
#if 0
	{
		FILE *fp=fopen("dump_TNB.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,Xz,-,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f\n",
					Xx[i],Xy[i],Xz[i], Tx[i],Ty[i],Tz[i], Nx[i],Ny[i],Nz[i], Bx[i],By[i],Bz[i], dx[i],kv[i],tr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}

int crease::calcTNB( int Xsidx0, int Xeidx0 )
{
	int ret=0;
	int size=Xeidx0-Xsidx0+1;

	calcTNB( size, &(Xx[Xsidx0]), &(Xy[Xsidx0]), &(Xz[Xsidx0]),
		&(Tx[Xsidx0]), &(Ty[Xsidx0]), &(Tz[Xsidx0]), &(dx[Xsidx0]),
		&(Nx[Xsidx0]), &(Ny[Xsidx0]), &(Nz[Xsidx0]), &(kv[Xsidx0]),
		&(Bx[Xsidx0]), &(By[Xsidx0]), &(Bz[Xsidx0]), &(tr[Xsidx0]) );
	//Xsidx = Xsidx0;
	//Xeidx = Xeidx0;
#if 0
	{
		FILE *fp=fopen("output/dump_TNB.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,Xz,-,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f\n",
					Xx[i],Xy[i],Xz[i], Tx[i],Ty[i],Tz[i], Nx[i],Ny[i],Nz[i], Bx[i],By[i],Bz[i], dx[i],kv[i],tr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}

int crease::calcTNB( int Xcnt0, double *Xx0, double *Xy0, double *Xz0,
					double *Tx0, double *Ty0, double *Tz0,	double *dx0,
					double *Nx0, double *Ny0, double *Nz0, double *kv0,
					double *Bx0, double *By0, double *Bz0, double *tr0 )
{
	int i, sign=1, ret=0;
	double dx1[MAX_SPCNT], dx[MAX_SPCNT], dy[MAX_SPCNT], dz[MAX_SPCNT];
	memset( Tx0, 0, sizeof(double)*Xcnt0 );
	memset( Ty0, 0, sizeof(double)*Xcnt0 );
	memset( Tz0, 0, sizeof(double)*Xcnt0 );
	memset( Nx0, 0, sizeof(double)*Xcnt0 );
	memset( Ny0, 0, sizeof(double)*Xcnt0 );
	memset( Nz0, 0, sizeof(double)*Xcnt0 );
	memset( Bx0, 0, sizeof(double)*Xcnt0 );
	memset( By0, 0, sizeof(double)*Xcnt0 );
	memset( Bz0, 0, sizeof(double)*Xcnt0 );
	memset( dx0, 0, sizeof(double)*Xcnt0 );
	memset( kv0, 0, sizeof(double)*Xcnt0 );
	memset( tr0, 0, sizeof(double)*Xcnt0 );

	for( i=0; i<Xcnt0-1; i++ ){
		dx[i] = Xx0[i+1] - Xx0[i];
		dy[i] = Xy0[i+1] - Xy0[i];
		dz[i] = Xz0[i+1] - Xz0[i];
		dx1[i] = sqrt( dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i] );
	}

	for( i=1; i<Xcnt0-1; i++ ){
		if( dx1[i-1]==0.0 && dx1[i]==0.0 ){
			Tx0[i] = Ty0[i] = Tz0[i] = dx0[i] = 0.0;
			continue;
		}
		if( dx1[i-1]>0.0 && dx1[i]>0.0 ){
		Tx0[i] = (dx[i-1]/dx1[i-1] + dx[i]/dx1[i])*0.5;
		Ty0[i] = (dy[i-1]/dx1[i-1] + dy[i]/dx1[i])*0.5;
		Tz0[i] = (dz[i-1]/dx1[i-1] + dz[i]/dx1[i])*0.5;
		dx0[i] = (dx1[i-1]+dx1[i])*0.5;
		} else if( dx1[i-1]>0.0 ){
			Tx0[i] = dx[i-1];
			Ty0[i] = dy[i-1];
			Tz0[i] = dz[i-1];
			dx0[i] = dx1[i-1];
		} else{
			Tx0[i] = dx[i];
			Ty0[i] = dy[i];
			Tz0[i] = dz[i];
			dx0[i] = dx1[i];
		}
		normalize_v3( &Tx0[i], &Ty0[i], &Tz0[i] );
	}
	for( i=2; i<Xcnt0-2; i++ ){
		double d0 = dx1[i];
		double d1 = dx1[i-1];
		double nx0 = Tx0[i+1] - Tx0[i];
		double nx1 = Tx0[i]   - Tx0[i-1];
		double ny0 = Ty0[i+1] - Ty0[i];
		double ny1 = Ty0[i]   - Ty0[i-1];
		double nz0 = Tz0[i+1] - Tz0[i];
		double nz1 = Tz0[i]   - Tz0[i-1];

		if( d0<=0.0 && d1<=0.0 ){
			Nx0[i] = Ny0[i] = Nz0[i] = kv0[i] = 0.0;
			continue;
		}
		if( d0>0.0 && d1>0.0 ){
		Nx0[i] = (nx0/d0 + nx1/d1)*0.5;
		Ny0[i] = (ny0/d0 + ny1/d1)*0.5;
		Nz0[i] = (nz0/d0 + nz1/d1)*0.5;
		} else if( d0>0.0 ){
			Nx0[i] = nx0;
			Ny0[i] = ny0;
			Nz0[i] = nz0;
		} else {
			Nx0[i] = nx1;
			Ny0[i] = ny1;
			Nz0[i] = nz1;
		}
		kv0[i] = sqrt( Nx0[i]*Nx0[i] + Ny0[i]*Ny0[i] + Nz0[i]*Nz0[i] );
		if((i==2 && Nx0[i]*Nx0[0] + Ny0[i]*Ny0[0] + Nz0[i]*Nz0[0] < 0)
			|| (i>2 && Nx0[i]*Nx0[i-1] + Ny0[i]*Ny0[i-1] + Nz0[i]*Nz0[i-1] < 0))
		{
			kv0[i] = -kv0[i];
		}
		if( fabs(kv0[i]) > 0.0 /*ZERO*/ ){	// kv=0ïtãﬂÇ≈ïsà¿íËÇ…Ç»ÇÈÅB
			Nx0[i] /= kv0[i];
			Ny0[i] /= kv0[i];
			Nz0[i] /= kv0[i];
		}
	}
	for( i=2; i<Xcnt0-2; i++ ){
		// B=TÅ~N
		Bx0[i] = Ty0[i]*Nz0[i] - Tz0[i]*Ny0[i];
		By0[i] = Tz0[i]*Nx0[i] - Tx0[i]*Nz0[i];
		Bz0[i] = Tx0[i]*Ny0[i] - Ty0[i]*Nx0[i];
	}

	calcTau( Xcnt0, Xx0, Xy0, Xz0, Tx0, Ty0, Tz0, dx0, Nx0, Ny0, Nz0, kv0, Bx0, By0, Bz0, tr0 );
#if 0
	{
		FILE *fp=fopen("output/dump_TNB.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,Xz,-,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f\n",
					Xx[i],Xy[i],Xz[i], Tx[i],Ty[i],Tz[i], Nx[i],Ny[i],Nz[i], Bx[i],By[i],Bz[i], dx[i],kv[i],tr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}



int crease::calcTau()
{
	return calcTau( Xcnt, Xx, Xy, Xz, Tx, Ty, Tz, dx, Nx, Ny, Nz, kv, Bx, By, Bz, tr );
}

int crease::calcTau( int Xcnt0, double *Xx0, double *Xy0, double *Xz0,
					double *Tx0, double *Ty0, double *Tz0, double *dx0_,
					double *Nx0, double *Ny0, double *Nz0, double *kv0,
					double *Bx0, double *By0, double *Bz0, double *tr0 )
{
	int i, ret=0;
	memset( tr0, 0, sizeof(double)*Xcnt0 );

	FILE *fp=NULL;
	//fp=fopen("output/dump_Tau.csv","w");
	if( fp ){
		fprintf( fp, "d0,d1,,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t,-,dTx,dTy,dTz,-,dNx,dNy,dNz,-,dBx,dBy,dBz" );
		fprintf( fp, ",,(dN+kT)x,(dN+kT)y,(dN+kT)z,,len(dN+kT),(dN+kT).B,,'-dBx,'-dBy,'-dBz,,len(-dB),'-dB.N" );
		fprintf( fp, "\n" );
	}

	// N' = -k*T + t*B -> t = (N'+kT)/B
	// B' = -t*N -> t = -B'/N
	for( i=2; i<Xcnt0-2; i++ ){
		double dTx=0,dTy=0,dTz=0, dNx=0,dNy=0,dNz=0, dBx=0,dBy=0,dBz=0;
		double nktx=0,nkty=0,nktz=0, bx=0,by=0,bz=0, t0=0,t1=0, ip0=0,ip1=0;
		double dx0,dx1,dy0,dy1,dz0,dz1,d0,d1;
		dx0 = Xx0[i+1] - Xx0[i];
		dx1 = Xx0[i]   - Xx0[i-1];
		dy0 = Xy0[i+1] - Xy0[i];
		dy1 = Xy0[i]   - Xy0[i-1];
		dz0 = Xz0[i+1] - Xz0[i];
		dz1 = Xz0[i]   - Xz0[i-1];
		d0 = sqrt( dx0*dx0 + dy0*dy0 + dz0*dz0 );
		d1 = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1 );
		if( d0==0.0 || d1==0.0 ){
			if( fp ){
				fprintf( fp, "%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f\n",
					d0,d1,Tx0[i],Ty0[i],Tz0[i], Nx0[i],Ny0[i],Nz0[i], Bx0[i],By0[i],Bz0[i] );
			}
			continue;
		}

		dx0 = Tx0[i+1] - Tx0[i];
		dx1 = Tx0[i]   - Tx0[i-1];
		dy0 = Ty0[i+1] - Ty0[i];
		dy1 = Ty0[i]   - Ty0[i-1];
		dz0 = Tz0[i+1] - Tz0[i];
		dz1 = Tz0[i]   - Tz0[i-1];
		dTx = (dx0/d0 + dx1/d1)*0.5;
		dTy = (dy0/d0 + dy1/d1)*0.5;
		dTz = (dz0/d0 + dz1/d1)*0.5;

		dx0 = Nx0[i+1] - Nx0[i];
		dx1 = Nx0[i]   - Nx0[i-1];
		dy0 = Ny0[i+1] - Ny0[i];
		dy1 = Ny0[i]   - Ny0[i-1];
		dz0 = Nz0[i+1] - Nz0[i];
		dz1 = Nz0[i]   - Nz0[i-1];
		dNx = (dx0/d0 + dx1/d1)*0.5;
		dNy = (dy0/d0 + dy1/d1)*0.5;
		dNz = (dz0/d0 + dz1/d1)*0.5;

		dx0 = Bx0[i+1] - Bx0[i];
		dx1 = Bx0[i]   - Bx0[i-1];
		dy0 = By0[i+1] - By0[i];
		dy1 = By0[i]   - By0[i-1];
		dz0 = Bz0[i+1] - Bz0[i];
		dz1 = Bz0[i]   - Bz0[i-1];
		dBx = (dx0/d0 + dx1/d1)*0.5;
		dBy = (dy0/d0 + dy1/d1)*0.5;
		dBz = (dz0/d0 + dz1/d1)*0.5;

		// N' = -k*T + t*B -> t = (N'+kT)/B
		nktx = dNx+kv0[i]*Tx0[i];
		nkty = dNy+kv0[i]*Ty0[i];
		nktz = dNz+kv0[i]*Tz0[i];
		t0 = sqrt( nktx*nktx + nkty*nkty + nktz*nktz );
		if( t0!=0 ){
			nktx /= t0;
			nkty /= t0;
			nktz /= t0;
		} else {
			nktx = nkty = nktz = 0.0;
		}
		ip0 = nktx*Bx0[i] + nkty*By0[i] + nktz*Bz0[i];

		// B' = -t*N -> t = -B'/N
		bx = -dBx;
		by = -dBy;
		bz = -dBz;
		t1 = sqrt( bx*bx + by*by + bz*bz );
		if( t1!=0 ){
			bx /= t1;
			by /= t1;
			bz /= t1;
		} else {
			bx = by = bz = 0.0;
		}
		ip1 = bx*Nx0[i] + by*Ny0[i] + bz*Nz0[i];

		if( fabs(ip0) >= 0.99 && fabs(ip0) >= fabs(ip1) ){
			if( ip0>0 ){
				tr0[i] = t0;
			} else {
				tr0[i] = -t0;
			}
		} else if( fabs(ip1) >= 0.99 && fabs(ip1) >= fabs(ip0) ){
			if( ip1>0 ){
				tr0[i] = t1;
			} else {
				tr0[i] = -t1;
			}
		} else {
			tr0[i] = 0.0;
		}

		if( fp ){
			fprintf( fp, "%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f,,%f,%f,%f",
				d0,d1,Tx0[i],Ty0[i],Tz0[i], Nx0[i],Ny0[i],Nz0[i], Bx0[i],By0[i],Bz0[i], dx0_[i],kv0[i],tr0[i],
				dTx,dTy,dTz, dNx,dNy,dNz, dBx,dBy,dBz);
			fprintf( fp, ",,%f,%f,%f,,%f,%f,,%f,%f,%f,,%f,%f",
				nktx,nkty,nktz, t0,ip0, bx,by,bz, t1,ip1 );
			fprintf( fp, "\n" );
		}
	} // i

	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}

// kv,tr -> X*, T*, N*, B*
int crease::calcXTNB( double m3[16] )
{
	int i,ret=0;
	double len;

	memset( Xx, 0, sizeof(double)*MAX_SPCNT );
	memset( Xy, 0, sizeof(double)*MAX_SPCNT );
	memset( Xz, 0, sizeof(double)*MAX_SPCNT );
	memset( Tx, 0, sizeof(double)*MAX_SPCNT );
	memset( Ty, 0, sizeof(double)*MAX_SPCNT );
	memset( Tz, 0, sizeof(double)*MAX_SPCNT );
	memset( Nx, 0, sizeof(double)*MAX_SPCNT );
	memset( Ny, 0, sizeof(double)*MAX_SPCNT );
	memset( Nz, 0, sizeof(double)*MAX_SPCNT );
	memset( Bx, 0, sizeof(double)*MAX_SPCNT );
	memset( By, 0, sizeof(double)*MAX_SPCNT );
	memset( Bz, 0, sizeof(double)*MAX_SPCNT );

	if( kv[0]==0.0 && kv[1]==0.0 && kv[Xcnt-2]==0.0 && kv[Xcnt-1]==0.0 ){
		kv[0] = kv[1] = kv[2];
		kv[Xcnt-1] = kv[Xcnt-2] = kv[Xcnt-3];
	}
	if( tr[0]==0.0 && tr[1]==0.0 && tr[2]==0.0
		&& tr[Xcnt-1]==0.0 && tr[Xcnt-2]==0.0 && tr[Xcnt-3]==0.0 ){
			tr[0] = tr[1] = tr[2] = tr[3];
			tr[Xcnt-1] = tr[Xcnt-2] = tr[Xcnt-3] = tr[Xcnt-4];
	}
	for( i=0; i<Xcnt; i++ ){
		dx[i] = XSPC;
	}

	Xx[0] = m3[12]; Xy[0] = m3[13]; Xz[0] = m3[14];
	Tx[0] = m3[0];  Ty[0] = m3[1];  Tz[0] = m3[2];
	Nx[0] = m3[4];  Ny[0] = m3[5];  Nz[0] = m3[6]; 
	Bx[0] = m3[8];  By[0] = m3[9];  Bz[0] = m3[10];

	// T[1] = kv[0]*dx[0]*N[0]
	// N[1] = - kv[0]*dx[0]*T[0] + tr[0]*dx[0]*B[0]
	// B=TÅ~N
	Tx[1] = Tx[0] + kv[0]*dx[0]*Nx[0];
	Ty[1] = Ty[0] + kv[0]*dx[0]*Ny[0];
	Tz[1] = Tz[0] + kv[0]*dx[0]*Nz[0];
	len = sqrt(Tx[1]*Tx[1] + Ty[1]*Ty[1] + Tz[1]*Tz[1]);
	Tx[1] /= len;	Ty[1] /= len;	Tz[1] /= len;

	// (T[0]+T[1])/2
	double tx01 = Tx[0]+Tx[1];
	double ty01 = Ty[0]+Ty[1];
	double tz01 = Tz[0]+Tz[1];
	len = sqrt( tx01*tx01 + ty01*ty01 + tz01*tz01 );
	tx01 /= len;	ty01 /= len;	tz01 /= len;
	Xx[1] = Xx[0] + tx01*dx[0];
	Xy[1] = Xy[0] + ty01*dx[0];
	Xz[1] = Xz[0] + tz01*dx[0];

	Nx[1] = Nx[0] - kv[0]*dx[0]*Tx[0] + tr[0]*dx[0]*Bx[0];
	Ny[1] = Ny[0] - kv[0]*dx[0]*Ty[0] + tr[0]*dx[0]*By[0];
	Nz[1] = Nz[0] - kv[0]*dx[0]*Tz[0] + tr[0]*dx[0]*Bz[0];
	len = sqrt(Nx[1]*Nx[1] + Ny[1]*Ny[1] + Nz[1]*Nz[1]);
	Nx[1] /= len;	Ny[1] /= len;	Nz[1] /= len;

	Bx[1] = Ty[1]*Nz[1] - Tz[1]*Ny[1];
	By[1] = Tz[1]*Nx[1] - Tx[1]*Nz[1];
	Bz[1] = Tx[1]*Ny[1] - Ty[1]*Nx[1];
	len = sqrt(Bx[1]*Bx[1] + By[1]*By[1] + Bz[1]*Bz[1]);
	Bx[1] /= len;	By[1] /= len;	Bz[1] /= len;

	for( i=2; i<Xcnt; i++ ){
		double pdx=Xx[i-1]-Xx[i-2], pdy=Xy[i-1]-Xy[i-2], pdz=Xz[i-1]-Xz[i-2], pdl=sqrt(pdx*pdx+pdy*pdy+pdz*pdz);
		double ip = (Xx[i-1]-Xx[i-2])*Tx[i-1] + (Xy[i-1]-Xy[i-2])*Ty[i-1] + (Xz[i-1]-Xz[i-2])*Tz[i-1];
		Xx[i] = Xx[i-2] + 2*ip*Tx[i-1];
		Xy[i] = Xy[i-2] + 2*ip*Ty[i-1];
		Xz[i] = Xz[i-2] + 2*ip*Tz[i-1];

		Tx[i] = Tx[i-2] + 2*kv[i-1]*dx[i-1]*Nx[i-1];
		Ty[i] = Ty[i-2] + 2*kv[i-1]*dx[i-1]*Ny[i-1];
		Tz[i] = Tz[i-2] + 2*kv[i-1]*dx[i-1]*Nz[i-1];
		len = sqrt(Tx[i]*Tx[i] + Ty[i]*Ty[i] + Tz[i]*Tz[i]);
		Tx[i] /= len;	Ty[i] /= len;	Tz[i] /= len;

		Nx[i] = Nx[i-2] - 2*kv[i-1]*dx[i-1]*Tx[i-1] + 2*tr[i-1]*dx[i-1]*Bx[i-1];
		Ny[i] = Ny[i-2] - 2*kv[i-1]*dx[i-1]*Ty[i-1] + 2*tr[i-1]*dx[i-1]*By[i-1];
		Nz[i] = Nz[i-2] - 2*kv[i-1]*dx[i-1]*Tz[i-1] + 2*tr[i-1]*dx[i-1]*Bz[i-1];
		len = sqrt(Nx[i]*Nx[i] + Ny[i]*Ny[i] + Nz[i]*Nz[i]);
		Nx[i] /= len;	Ny[i] /= len;	Nz[i] /= len;

		Bx[i] = Ty[i]*Nz[i] - Tz[i]*Ny[i];
		By[i] = Tz[i]*Nx[i] - Tx[i]*Nz[i];
		Bz[i] = Tx[i]*Ny[i] - Ty[i]*Nx[i];
		len = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);
		Bx[i] /= len;	By[i] /= len;	Bz[i] /= len;
	}

#if 1
	double eTx[MAX_SPCNT], eTy[MAX_SPCNT], eTz[MAX_SPCNT];
	for( i=0; i<Xcnt-1; i++ ){
		eTx[i] = Tx[i] + Tx[i+1];
		eTy[i] = Ty[i] + Ty[i+1];
		eTz[i] = Tz[i] + Tz[i+1];
		len = sqrt( eTx[i]*eTx[i] + eTy[i]*eTy[i] + eTz[i]*eTz[i] );
		eTx[i] /= len;	eTy[i] /= len;	eTz[i] /= len;
	}
	for( i=1; i<Xcnt; i++ ){
		Xx[i] = Xx[i-1] + eTx[i-1]*dx[i-1];
		Xy[i] = Xy[i-1] + eTy[i-1]*dx[i-1];
		Xz[i] = Xz[i-1] + eTz[i-1]*dx[i-1];
	}
#endif

#if 0
	{
		FILE *fp=fopen("output/dump_calcXTNB0.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,Xz,-,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					Xx[i],Xy[i],Xz[i], sqrt( Xx[i]*Xx[i] + Xy[i]*Xy[i] + Xz[i]*Xz[i] ),
					Tx[i],Ty[i],Tz[i], sqrt( Tx[i]*Tx[i] + Ty[i]*Ty[i] + Tz[i]*Tz[i] ),
					Nx[i],Ny[i],Nz[i], sqrt( Nx[i]*Nx[i] + Ny[i]*Ny[i] + Nz[i]*Nz[i] ),
					Bx[i],By[i],Bz[i], sqrt( Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i] ),
					dx[i],kv[i],tr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
	// T,N,B,dx,kv,tr ÇãÅÇﬂÇ»Ç®Ç∑
	calcTNB();

#if 0
	{
		FILE *fp=fopen("output/dump_calcXTNB1.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,Xz,-,Tx,Ty,Tz,-,Nx,Ny,Nz,-,Bx,By,Bz,-,d,k,t\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					Xx[i],Xy[i],Xz[i], sqrt( Xx[i]*Xx[i] + Xy[i]*Xy[i] + Xz[i]*Xz[i] ),
					Tx[i],Ty[i],Tz[i], sqrt( Tx[i]*Tx[i] + Ty[i]*Ty[i] + Tz[i]*Tz[i] ),
					Nx[i],Ny[i],Nz[i], sqrt( Nx[i]*Nx[i] + Ny[i]*Ny[i] + Nz[i]*Nz[i] ),
					Bx[i],By[i],Bz[i], sqrt( Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i] ),
					dx[i],kv[i],tr[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	return ret;
}

int crease::calcXTN2d( double m2[9] )
{
	int i,ret=0;
	double len;
	memset( Xx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Xy2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Tx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Ty2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Nx2d, 0, sizeof(double)*MAX_SPCNT );
	memset( Ny2d, 0, sizeof(double)*MAX_SPCNT );

	if( k2d[0]==0.0 && k2d[1]==0.0 && k2d[Xcnt-2]==0.0 && k2d[Xcnt-1]==0.0 ){
		k2d[0] = k2d[1] = k2d[2];
		k2d[Xcnt-1] = k2d[Xcnt-2] = k2d[Xcnt-3];
	}
	for( i=0; i<Xcnt; i++ ){
		d2d[i] = XSPC;
	}

	Xx2d[0] = m2[6]; Xy2d[0] = m2[7];
	Tx2d[0] = m2[0]; Ty2d[0] = m2[1];
	Nx2d[0] = m2[3]; Ny2d[0] = m2[4];

	Tx2d[1] = Tx2d[0] + k2d[0]*d2d[0]*Nx2d[0];
	Ty2d[1] = Ty2d[0] + k2d[0]*d2d[0]*Ny2d[0];
	len = sqrt(Tx2d[1]*Tx2d[1] + Ty2d[1]*Ty2d[1]);
	Tx2d[1] /= len;	Ty2d[1] /= len;

	double tx01 = Tx2d[0]+Tx2d[1];
	double ty01 = Ty2d[0]+Ty2d[1];
	len = sqrt(tx01*tx01+ty01*ty01);	tx01 /= len;	ty01 /= len;
	Xx2d[1] = Xx2d[0] + tx01*d2d[0];
	Xy2d[1] = Xy2d[0] + ty01*d2d[0];
	Nx2d[1] = -Ty2d[1];
	Ny2d[1] = Tx2d[1];

	for( i=2; i<Xcnt; i++ ){
		double ip = (Xx2d[i-1]-Xx2d[i-2])*Tx2d[i-1] + (Xy2d[i-1]-Xy2d[i-2])*Ty2d[i-1];
		Xx2d[i] = Xx2d[i-2] + 2*ip*Tx2d[i-1];
		Xy2d[i] = Xy2d[i-2] + 2*ip*Ty2d[i-1];

		Tx2d[i] = Tx2d[i-2] + 2*k2d[i-1]*d2d[i-1]*Nx2d[i-1];
		Ty2d[i] = Ty2d[i-2] + 2*k2d[i-1]*d2d[i-1]*Ny2d[i-1];
		len = sqrt(Tx2d[i]*Tx2d[i] + Ty2d[i]*Ty2d[i]);
		Tx2d[i] /= len;	Ty2d[i] /= len;

		// îΩéûåvâÒÇËÇ…90ìxâÒì]
		Nx2d[i] = -Ty2d[i];
		Ny2d[i] = Tx2d[i];
	}
#if 1
	double eTx[MAX_SPCNT], eTy[MAX_SPCNT];
	for( i=0; i<Xcnt-1; i++ ){
		eTx[i] = Tx2d[i] + Tx2d[i+1];
		eTy[i] = Ty2d[i] + Ty2d[i+1];
		len = sqrt( eTx[i]*eTx[i] + eTy[i]*eTy[i] );
		eTx[i] /= len;	eTy[i] /= len;
	}
	for( i=1; i<Xcnt; i++ ){
		Xx2d[i] = Xx2d[i-1] + eTx[i-1]*d2d[i-1];
		Xy2d[i] = Xy2d[i-1] + eTy[i-1]*d2d[i-1];
	}
#endif

	// T,N,d2d,k2d ÇãÅÇﬂÇ»Ç®Ç∑
	calcTN2d();

#if 0
	{
		FILE *fp=fopen("output/dump_calcXTN2d.csv","w");
		if( fp ){
			fprintf( fp, "Xx,Xy,-,Tx,Ty,-,Nx,Ny,-,d,k2d\n" );
			for( i=0; i<Xcnt; i++ ){
				fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					Xx2d[i],Xy2d[i], sqrt( Xx2d[i]*Xx2d[i] + Xy2d[i]*Xy2d[i] ),
					Tx2d[i],Ty2d[i], sqrt( Tx2d[i]*Tx2d[i] + Ty2d[i]*Ty2d[i] ),
					Nx2d[i],Ny2d[i], sqrt( Nx2d[i]*Nx2d[i] + Ny2d[i]*Ny2d[i] ),
					d2d[i],k2d[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif
end:
	return ret;
}

