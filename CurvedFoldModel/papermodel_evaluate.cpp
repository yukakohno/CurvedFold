#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "papermodel.h"
#include "util.h"

// 点(x0,y0)と線分(x1,y1)-(x2,y2)の距離
double min_d2( double x0, double y0, double z0, double x1, double y1, double z1,
			  double x2, double y2, double z2, double *x1_, double *y1_, double *z1_ )
{
	double x12 = x2-x1, x01 = x1-x0, x02 = x2-x0;
	double y12 = y2-y1, y01 = y1-y0, y02 = y2-y0;
	double z12 = z2-z1, z01 = z1-z0, z02 = z2-z0;
	double l2_12 = x12*x12 + y12*y12 + z12*z12;
	double l2_01 = x01*x01 + y01*y01 + z01*z01;
	double l2_02 = x02*x02 + y02*y02 + z02*z02;
	double tt = -(x12*x01 + y12*y01 + z12*z01);
	if( tt < 0 ) {
		if( x1_ && y1_ && z1_ ){ *x1_ = x1; *y1_ = y1; *z1_ = z1; }
		//return -1;
		return x01*x01 + y01*y01 + z01*z01;
	}
	if( tt > l2_12 ) { // tt/len12 > len12
		if( x1_ && y1_ && z1_ ){ *x1_ = x2; *y1_ = y2; *z1_ = z2; }
		//return -1;
		return x02*x02 + y02*y02 + z02*z02;
	}
	if( x1_ && y1_ && z1_ ){
		*x1_ = x1 + tt*x12/l2_12; // tt/len12 * x12/len12
		*y1_ = y1 + tt*y12/l2_12;
		*z1_ = z1 + tt*z12/l2_12;
	}
	//double f1 = a*(y1-y0)-b*(x1-x0);
	//AB、APを外積して求められたベクトルの長さが、平行四辺形Dの面積になる
	double cvx = (y12 * z01) - (z12 * y01);
	double cvy = (z12 * x01) - (x12 * z01);
	double cvz = (x12 * y01) - (y12 * x01);
	double f2 = cvx*cvx + cvy*cvy + cvz*cvz;
	return f2/l2_12;
}

int papermodel::checkGap( double *errdata, int row, int col, int divnum )
{
	int ret=0, spcnt, p3cnt0, p3cnt1;
	double space = 1, p3len0[MAX_STP_CNT], p3len1[MAX_STP_CNT];
	double stx0[MAX_STP_CNT], sty0[MAX_STP_CNT], stz0[MAX_STP_CNT];
	double stx1[MAX_STP_CNT], sty1[MAX_STP_CNT], stz1[MAX_STP_CNT];
	double stx1_[MAX_STP_CNT], sty1_[MAX_STP_CNT], stz1_[MAX_STP_CNT];
	double dif[MAX_STP_CNT];

	if( lcrcnt<=1 && rcrcnt<=1 ){
		return -1;
	}

	p3cnt0 = lcrs[1]->getSamplePoints3D( space, stx0, sty0, stz0, p3len0 );
	//p3cnt1 = rcrs[1]->getSamplePoints3D( space, stx1, sty1, stz1, p3len1 );
	p3cnt1 = rcrs[1]->getVertexPoints3D( space, stx1, sty1, stz1, p3len1 );

	//spcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
	printf( "p3cnt0 = %d, p3cnt1 = %d, p3len0 = %f, p3len1 = %f\n",
		p3cnt0, p3cnt1, p3len0[p3cnt0-1], p3len1[p3cnt1-1] );

	// rotation & mid points
	double axis[3]={0.0,0.0,1.0}, quat[4], rmat[16];
	axis_ang_quat( axis, 2.0*M_PI/(double)divnum, quat );
	quat_mat( quat, rmat );

	for( int i=0; i<p3cnt1; i++ )
	{
		double tmpx, tmpy, tmpz;
		tmpx = rmat[0]*stx1[i] + rmat[4]*sty1[i] + rmat[8]*stz1[i] + rmat[12];
		tmpy = rmat[1]*stx1[i] + rmat[5]*sty1[i] + rmat[9]*stz1[i] + rmat[13];
		tmpz = rmat[2]*stx1[i] + rmat[6]*sty1[i] + rmat[10]*stz1[i] + rmat[14];
		stx1[i] = tmpx;
		sty1[i] = tmpy;
		stz1[i] = tmpz;
	}

	spcnt = p3cnt0;
	double d0 = 0.0;
	for( int i=0; i<p3cnt0; i++ )
	{
		dif[i] = 10000;	// max
		if( p3len0[i] > p3len1[p3cnt1-1] ){
			spcnt = i;
			break;
		}
		int minj=-1;
		for( int j=0; j<p3cnt1-1; j++ )
		{
			double diff, minx1, miny1, minz1;
			diff = min_d2( stx0[i], sty0[i], stz0[i], stx1[j], sty1[j], stz1[j],
				stx1[j+1], sty1[j+1], stz1[j+1], &minx1, &miny1, &minz1 );
			//dif[i] = dif[i] < diff ? dif[i] : diff;	// MIN( dif[i], diff )
			if( diff > 0 && dif[i] > diff ){
				dif[i] = diff;
				minj = j;
				stx1_[i] = minx1;
				sty1_[i] = miny1;
				stz1_[i] = minz1;
			}
		}
		if( dif[i] == 10000 ){
			errdata[i*col+0] = dif[i] = -1;
			//spcnt = i;
			//break;
		} else {
			errdata[i*col+0] =  sqrt(dif[i]);
		}
		errdata[i*col+1] =  stx0[i];
		errdata[i*col+2] =  sty0[i];
		errdata[i*col+3] =  stz0[i];
		errdata[i*col+4] =  stx1_[i];
		errdata[i*col+5] =  sty1_[i];
		errdata[i*col+6] =  stz1_[i];
	}

#if 0
	{
		FILE *fp = fopen("output/checkgap.csv", "w");
		if( fp ){
			fprintf( fp, "stx0,sty0,stz0,len0,stx1,sty1,stz1,len1,diff,stx1_,sty1_,stz1_\n" );
			for( int i=0; i<spcnt; i++ ){
				double diff = -1;
				if( dif[i]>0 ){ diff = sqrt(dif[i]); }
				fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
					stx0[i], sty0[i], stz0[i], p3len0[i], stx1[i], sty1[i], stz1[i], p3len1[i],
					diff, stx1_[i], sty1_[i], stz1_[i] );
			}
			int mincnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
			for( int i=spcnt; i<mincnt; i++ ){
				fprintf( fp, "%f,%f,%f,%f,%f,%f,%f,%f\n",
					stx0[i], sty0[i], stz0[i], p3len0[i], stx1[i], sty1[i], stz1[i], p3len1[i] );
			}
			if( mincnt==p3cnt1 ){
				for( int i=mincnt; i<p3cnt0; i++ ){
					fprintf( fp, "%f,%f,%f,%f,,,,\n", stx0[i], sty0[i], stz0[i], p3len0[i] );
				}
			} else {
				for( int i=mincnt; i<p3cnt1; i++ ){
					fprintf( fp, ",,,,%f,%f,%f,%f\n", stx1[i], sty1[i], stz1[i], p3len1[i] );
				}
			}
			fclose( fp );	fp=NULL;
		}
	}
#endif

	return spcnt;
}

/* http://marupeke296.com/COL_3D_No21_TriTri.html */
int collisiontest( double x00, double y00, double z00, double x01, double y01, double z01,
				  double x02, double y02, double z02, double nx0, double ny0, double nz0,
				  double x10, double y10, double z10, double x11, double y11, double z11,
				  double x12, double y12, double z12, double nx1, double ny1, double nz1 )
{
	int ret=0, mcnt0=0, mcnt1=0;
	double vx00,vy00,vz00, vx01,vy01,vz01, vx02,vy02,vz02, vx10,vy10,vz10, vx11,vy11,vz11, vx12,vy12,vz12;
	double ip00, ip01, ip02, ip10, ip11, ip12, d00, d01, d02, d10, d11, d12;
	double mx0[2],my0[2],mz0[2], mx1[2],my1[2],mz1[2], min0,max0,min1,max1; 
	double zero=0.00001;

	// test
	vx01 = x01-x00; vy01 = y01-y00; vz01 = z01-z00;
	vx02 = x02-x00; vy02 = y02-y00; vz02 = z02-z00;
	ip01 = vx01*nx0 + vy01*ny0 + vz01*nz0;
	ip02 = vx02*nx0 + vy02*ny0 + vz02*nz0;
	vx11 = x11-x10; vy11 = y11-y10; vz11 = z11-z10;
	vx12 = x12-x10; vy12 = y12-y10; vz12 = z12-z10;
	ip11 = vx11*nx1 + vy11*ny1 + vz11*nz1;
	ip12 = vx12*nx1 + vy12*ny1 + vz12*nz1;
	//printf("ip = %f, %f, %f, %f\n", ip01, ip02, ip11, ip12 );


	// すべての頂点が片側にいるか？
	vx00 = x00-x10; vy00 = y00-y10; vz00 = z00-z10;
	vx01 = x01-x10; vy01 = y01-y10; vz01 = z01-z10;
	vx02 = x02-x10; vy02 = y02-y10; vz02 = z02-z10;
	ip00 = vx00*nx1 + vy00*ny1 + vz00*nz1;	d00 = fabs(ip00);
	ip01 = vx01*nx1 + vy01*ny1 + vz01*nz1;	d01 = fabs(ip01);
	ip02 = vx02*nx1 + vy02*ny1 + vz02*nz1;	d02 = fabs(ip02);
	if( ip00>=zero && ip01>=zero && ip02>=zero || ip00<=-zero && ip01<=-zero && ip02<=-zero ){
		goto end;
	}

	vx10 = x10-x00; vy10 = y10-y00; vz10 = z10-z00;
	vx11 = x11-x00; vy11 = y11-y00; vz11 = z11-z00;
	vx12 = x12-x00; vy12 = y12-y00; vz12 = z12-z00;
	ip10 = vx10*nx0 + vy10*ny0 + vz10*nz0;	d10 = fabs(ip10);
	ip11 = vx11*nx0 + vy11*ny0 + vz11*nz0;	d11 = fabs(ip11);
	ip12 = vx12*nx0 + vy12*ny0 + vz12*nz0;	d12 = fabs(ip12);
	if( ip10>=zero && ip11>=zero && ip12>=zero || ip10<=-zero && ip11<=-zero && ip12<=-zero ){
		goto end;
	}

	if( d00<=zero && mcnt0<2 ){	mx0[mcnt0]=x00; my0[mcnt0]=y00; mz0[mcnt0]=z00; mcnt0++; }
	if( d01<=zero && mcnt0<2 ){	mx0[mcnt0]=x01; my0[mcnt0]=y01; mz0[mcnt0]=z01; mcnt0++; }
	if( d02<=zero && mcnt0<2 ){	mx0[mcnt0]=x02; my0[mcnt0]=y02; mz0[mcnt0]=z02; mcnt0++; }
	if(( ip00>zero && ip01<-zero || ip00<-zero && ip01>zero ) && mcnt0<2 ){
		double t = d00/(d00+d01);
		mx0[mcnt0] = vx00 + t*(vx01-vx00) +x10;
		my0[mcnt0] = vy00 + t*(vy01-vy00) +y10;
		mz0[mcnt0] = vz00 + t*(vz01-vz00) +z10;
		mcnt0++;
	}
	if(( ip01>zero && ip02<-zero || ip01<-zero && ip02>zero ) && mcnt0<2 ){
		double t = d01/(d01+d02);
		mx0[mcnt0] = vx01 + t*(vx02-vx01) +x10;
		my0[mcnt0] = vy01 + t*(vy02-vy01) +y10;
		mz0[mcnt0] = vz01 + t*(vz02-vz01) +z10;
		mcnt0++;
	}
	if(( ip02>zero && ip00<-zero || ip02<-zero && ip00>zero ) && mcnt0<2 ){
		double t = d02/(d02+d00);
		mx0[mcnt0] = vx02 + t*(vx00-vx02) +x10;
		my0[mcnt0] = vy02 + t*(vy00-vy02) +y10;
		mz0[mcnt0] = vz02 + t*(vz00-vz02) +z10;
		mcnt0++;
	}
	if( mcnt0<2 ){
		goto end;
	}

	if( d10<=zero && mcnt1<2 ){	mx1[mcnt1]=x10; my1[mcnt1]=y10; mz1[mcnt1]=z10; mcnt1++; }
	if( d11<=zero && mcnt1<2 ){	mx1[mcnt1]=x11; my1[mcnt1]=y11; mz1[mcnt1]=z11; mcnt1++; }
	if( d12<=zero && mcnt1<2 ){	mx1[mcnt1]=x12; my1[mcnt1]=y12; mz1[mcnt1]=z12; mcnt1++; }
	if(( ip10>zero && ip11<-zero || ip10<-zero && ip11>zero ) && mcnt1<2 ){
		double t = d10/(d10+d11);
		mx1[mcnt1] = vx10 + t*(vx11-vx10) +x00;
		my1[mcnt1] = vy10 + t*(vy11-vy10) +y00;
		mz1[mcnt1] = vz10 + t*(vz11-vz10) +z00;
		mcnt1++;
	}
	if(( ip11>zero && ip12<-zero || ip11<-zero && ip12>zero ) && mcnt1<2 ){
		double t = d11/(d11+d12);
		mx1[mcnt1] = vx11 + t*(vx12-vx11) +x00;
		my1[mcnt1] = vy11 + t*(vy12-vy11) +y00;
		mz1[mcnt1] = vz11 + t*(vz12-vz11) +z00;
		mcnt1++;
	}
	if(( ip12>zero && ip10<-zero || ip12<-zero && ip10>zero ) && mcnt1<2 ){
		double t = d12/(d12+d10);
		mx1[mcnt1] = vx12 + t*(vx10-vx12) +x00;
		my1[mcnt1] = vy12 + t*(vy10-vy12) +y00;
		mz1[mcnt1] = vz12 + t*(vz10-vz12) +z00;
		mcnt1++;
	}
	if( mcnt1<2 ){
		goto end;
	}
	//printf( "mcnt0=%d, mcnt1=%d\n", mcnt0, mcnt1 );
	// test
	{
		double mx001, my001, mz001, mx101, my101, mz101, len0, len1, ip;
		mx001 = mx0[1]-mx0[0];
		my001 = my0[1]-my0[0];
		mz001 = mz0[1]-mz0[0];
		len0 = sqrt(mx001*mx001 + my001*my001 + mz001*mz001);
		mx001 /= len0;	my001 /= len0;	mz001 /= len0;
		mx101 = mx1[1]-mx1[0];
		my101 = my1[1]-my1[0];
		mz101 = mz1[1]-mz1[0];
		len1 = sqrt(mx101*mx101 + my101*my101 + mz101*mz101);
		mx101 /= len1;	my101 /= len1;	mz101 /= len1;
		ip = fabs(mx001*mx101 + my001*my101 + mz001*mz101);
		if( ip<0.99 ){
			printf( "m0 = ( %0.2f, %0.2f, %0.2f ), m1 = ( %0.2f, %0.2f, %0.2f )\n",
				mx001, my001, mz001, mx101, my101, mz101 );
		}
		goto end;
	}

	if( mx0[0] != mx0[1] && mx1[0] != mx1[1] ){
		min0 = mx0[0] < mx0[1] ? mx0[0] : mx0[1];
		max0 = mx0[0] > mx0[1] ? mx0[0] : mx0[1];
		min1 = mx1[0] < mx1[1] ? mx1[0] : mx1[1];
		max1 = mx1[0] > mx1[1] ? mx1[0] : mx1[1];
	} else if( my0[0] != my0[1] && my1[0] != my1[1] ){
		min0 = my0[0] < my0[1] ? my0[0] : my0[1];
		max0 = my0[0] > my0[1] ? my0[0] : my0[1];
		min1 = my1[0] < my1[1] ? my1[0] : my1[1];
		max1 = my1[0] > my1[1] ? my1[0] : my1[1];
	} else if( mz0[0] != mz0[1] && mz1[0] != mz1[1] ){
		min0 = mz0[0] < mz0[1] ? mz0[0] : mz0[1];
		max0 = mz0[0] > mz0[1] ? mz0[0] : mz0[1];
		min1 = mz1[0] < mz1[1] ? mz1[0] : mz1[1];
		max1 = mz1[0] > mz1[1] ? mz1[0] : mz1[1];
	} else {
		printf( "error: m0[0]==m0[1] || m1[0]==m1[1]\n" );
		goto end;
	}

	if( max0 < min1 || max1 < min0 ){
		goto end;
	}

	ret = 1;
end:
	return ret;
}

int nvec( double x0, double y0, double z0, double x1, double y1, double z1,
		 double x2, double y2, double z2, double *_nx, double *_ny, double *_nz )
{
	double dx1,dy1,dz1,len1, dx2,dy2,dz2,len2, nx,ny,nz,len;

	dx1 = x1-x0;
	dy1 = y1-y0;
	dz1 = z1-z0;
	len1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
	if( len1==0.0 ){ return -1; }
	dx1/=len1;	dy1/=len1;	dz1/=len1;

	dx2 = x2-x0;
	dy2 = y2-y0;
	dz2 = z2-z0;
	len2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
	if( len2==0.0 ){ return -1; }
	dx2/=len2;	dy2/=len2;	dz2/=len2;

	nx = dy1*dz2-dz1*dy2;
	ny = dz1*dx2-dx1*dz2;
	nz = dx1*dy2-dy1*dx2;
	len = sqrt(nx*nx+ny*ny+nz*nz);
	if( len==0.0 ){ return -1; }
	nx/=len;	ny/=len;	nz/=len;

	if( _nx && _ny && _nz ){
		*_nx = nx;
		*_ny = ny;
		*_nz = nz;
	}
}

int collisiontest( int vcnt0, double *x0, double *y0, double *z0, double nx0, double ny0, double nz0,
				  int vcnt1, double *x1, double *y1, double *z1, double nx1, double ny1, double nz1 )
{
	int ret = 0;

	if( vcnt0==3 && vcnt1==3 ){	// 012-012
		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx0,ny0,nz0,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx1,ny1,nz1 );
	} else if( vcnt0==3 && vcnt1==4 ){	// 012-012, 012-230
		nvec( x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], &nx1,&ny1,&nz1 );
		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx0,ny0,nz0,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx1,ny1,nz1 );
		if( ret ){ goto end; }
		nvec( x1[2],y1[2],z1[2], x1[3],y1[3],z1[3], x1[0],y1[0],z1[0], &nx1,&ny1,&nz1 );
		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx0,ny0,nz0,
							 x1[2],y1[2],z1[2], x1[3],y1[3],z1[3], x1[0],y1[0],z1[0], nx1,ny1,nz1 );
	} else if( vcnt0==4 && vcnt1==3 ){	// 012-012, 230-012
		nvec( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], &nx0,&ny0,&nz0 );
		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx0,ny0,nz0,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx1,ny1,nz1 );
		if( ret ){ goto end; }
		nvec( x0[2],y0[2],z0[2], x0[3],y0[3],z0[3], x0[0],y0[0],z0[0], &nx0,&ny0,&nz0 );
		ret = collisiontest( x0[2],y0[2],z0[2], x0[3],y0[3],z0[3], x0[0],y0[0],z0[0], nx0,ny0,nz0,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx1,ny1,nz1 );
	} else if( vcnt0==4 && vcnt1==4 ){	// 012-012, 012-230, 230-012, 230-230
		double nx00,ny00,nz00, nx01,ny01,nz01, nx10,ny10,nz10, nx11,ny11,nz11; 
		nvec( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], &nx00,&ny00,&nz00 );
		nvec( x0[2],y0[2],z0[2], x0[3],y0[3],z0[3], x0[0],y0[0],z0[0], &nx01,&ny01,&nz01 );
		nvec( x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], &nx10,&ny10,&nz10 );
		nvec( x1[2],y1[2],z1[2], x1[3],y1[3],z1[3], x1[0],y1[0],z1[0], &nx11,&ny11,&nz11 );

		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx00,ny00,nz00,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx10,ny10,nz10 );
		if( ret ){ goto end; }
		ret = collisiontest( x0[0],y0[0],z0[0], x0[1],y0[1],z0[1], x0[2],y0[2],z0[2], nx00,ny00,nz00,
							 x1[2],y1[2],z1[2], x1[3],y1[3],z1[3], x1[0],y1[0],z1[0], nx11,ny11,nz11 );
		if( ret ){ goto end; }
		ret = collisiontest( x0[2],y0[2],z0[2], x0[3],y0[3],z0[3], x0[0],y0[0],z0[0], nx01,ny01,nz01,
							 x1[0],y1[0],z1[0], x1[1],y1[1],z1[1], x1[2],y1[2],z1[2], nx10,ny10,nz10 );
		if( ret ){ goto end; }
		ret = collisiontest( x0[2],y0[2],z0[2], x0[3],y0[3],z0[3], x0[0],y0[0],z0[0], nx01,ny01,nz01,
							 x1[2],y1[2],z1[2], x1[3],y1[3],z1[3], x1[0],y1[0],z1[0], nx11,ny11,nz11 );
	}
end:
	return ret;
}

int papermodel::checkCollision( int *errdata, int row, int col, int divnum )
{
	int ret=0;
	memset( errdata, 0, sizeof(int)*row*col );

	// papermodel
	// int plcnt, plvcnt[MAX_PL_FCNT];
	// double plx[MAX_PL_FCNT*4], ply[MAX_PL_FCNT*4], plz[MAX_PL_FCNT*4];
	double plnx[MAX_PL_FCNT], plny[MAX_PL_FCNT], plnz[MAX_PL_FCNT];
	for( int i=0; i<plcnt; i++ ){
		plnx[i] = plmat[i*16+8];
		plny[i] = plmat[i*16+9];
		plnz[i] = plmat[i*16+10];
	}

	for( int i=0; i<plcnt; i++ ){
		if( pl_cridx[i]>0 ){
			continue;
		}
		for( int j=i+1; j<plcnt; j++ ){
			if( pl_cridx[j]>0 ){
				continue;
			}
			ret = collisiontest( plvcnt[i], &(plx[i*4]), &(ply[i*4]), &(plz[i*4]), plnx[i],plny[i],plnz[i],
								 plvcnt[j], &(plx[j*4]), &(ply[j*4]), &(plz[j*4]), plnx[j],plny[j],plnz[j] );
			if( ret ){
				errdata[i] = errdata[j] = 1;
			}
#if 0
			printf( "i=%d,j=%d: ", i, j );
			for( int k=0; k<plcnt; k++ ){
				printf( "%d ", errdata[k] );
			}
			printf( "\n" );
#endif
		}
	}

	// TODO: ↑を3個作成、総当たりで collision detection

end:
	return ret;
}

