#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include "papermodel.h"
#include "util.h"
#include "Curvature.h"

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y

int papermodel::getMat_stitch( int divnum, int maxpcnt, double *mObj )
{
	int ret=0;
	int p2cnt0, p2cnt1, p3cnt0, p3cnt1;
	double space = 10;
	double stx0[MAX_STP_CNT*2], sty0[MAX_STP_CNT*2], stz0[MAX_STP_CNT*2];
	double stx1[MAX_STP_CNT*2], sty1[MAX_STP_CNT*2], stz1[MAX_STP_CNT*2];

	if( lcrcnt<=1 && rcrcnt<=1 ){
		return -1;
	}

	//p2cnt0 = lcrs[1]->getSamplePoints2D( space, stx_cp[0], sty_cp[0] );
	//p2cnt1 = rcrs[1]->getSamplePoints2D( space, stx_cp[1], sty_cp[1] );
	p3cnt0 = lcrs[1]->getSamplePoints3D( space, stx[0], sty[0], stz[0], NULL );
	p3cnt1 = rcrs[1]->getSamplePoints3D( space, stx[1], sty[1], stz[1], NULL );
	stpcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
	if( maxpcnt>0 ){
		stpcnt = stpcnt < maxpcnt ? stpcnt : maxpcnt;
	}

	// rotation & mid points
	double axis[3]={0.0,0.0,1.0}, quat[4], rtmat0[16], rtmat1[16];
	axis_ang_quat( axis, 2.0*M_PI/(double)divnum, quat );	quat_mat( quat, rtmat0 );
	axis_ang_quat( axis, -2.0*M_PI/(double)divnum, quat );	quat_mat( quat, rtmat1 );

	for( int i=0; i<stpcnt; i++ ){
		int i2 = stpcnt+i;
		this->stx1[0][i] = stx1[i] = rtmat0[0]*stx[1][i] + rtmat0[4]*sty[1][i] + rtmat0[8]*stz[1][i] + rtmat0[12];
		this->sty1[0][i] = sty1[i] = rtmat0[1]*stx[1][i] + rtmat0[5]*sty[1][i] + rtmat0[9]*stz[1][i] + rtmat0[13];
		this->stz1[0][i] = stz1[i] = rtmat0[2]*stx[1][i] + rtmat0[6]*sty[1][i] + rtmat0[10]*stz[1][i] + rtmat0[14];
		this->stx1[1][i] = stx1[i2] = rtmat1[0]*stx[0][i] + rtmat1[4]*sty[0][i] + rtmat1[8]*stz[0][i] + rtmat1[12];
		this->sty1[1][i] = sty1[i2] = rtmat1[1]*stx[0][i] + rtmat1[5]*sty[0][i] + rtmat1[9]*stz[0][i] + rtmat1[13];
		this->stz1[1][i] = stz1[i2] = rtmat1[2]*stx[0][i] + rtmat1[6]*sty[0][i] + rtmat1[10]*stz[0][i] + rtmat1[14];
	}
	for( int i=0; i<stpcnt; i++ ){
		int i2 = stpcnt+i;
		this->stx0[0][i] = stx0[i] = stx[0][i];
		this->sty0[0][i] = sty0[i] = sty[0][i];
		this->stz0[0][i] = stz0[i] = stz[0][i];
		this->stx0[1][i] = stx0[i2] = stx[1][i];
		this->sty0[1][i] = sty0[i2] = sty[1][i];
		this->stz0[1][i] = stz0[i2] = stz[1][i];
	}
	for( int i=0; i<stpcnt*2; i++ ){
		this->stxm[i/stpcnt][i%stpcnt] = stx1[i] = (stx0[i] + stx1[i]) *0.5;
		this->stym[i/stpcnt][i%stpcnt] = sty1[i] = (sty0[i] + sty1[i]) *0.5;
		this->stzm[i/stpcnt][i%stpcnt] = stz1[i] = (stz0[i] + stz1[i]) *0.5;
	}
	if(0){
		FILE *fp=fopen("output/dumpstp.csv", "w");
		for( int i=0; i<stpcnt*2; i++ ){
			fprintf( fp, "%d,%f,%f,%f,%f,%f,%f\n", i,stx0[i],sty0[i],stz0[i],stx1[i],sty1[i],stz1[i]);
		}
		fclose(fp);
	}

	// getMat
	double mat[16];
	getMatRot( stpcnt*2, stx0, sty0, stz0, stx1, sty1, stz1, mat );
	mObj[0] = mat[0];	mObj[1] = mat[4];	mObj[2] = mat[8];	mObj[3] = mat[12];
	mObj[4] = mat[1];	mObj[5] = mat[5];	mObj[6] = mat[9];	mObj[7] = mat[13];
	mObj[8] = mat[2];	mObj[9] = mat[6];	mObj[10] = mat[10];	mObj[11] = mat[14];
	mObj[12] = mat[3];	mObj[13] = mat[7];	mObj[14] = mat[11];	mObj[15] = mat[15];

	return ret;
}

int papermodel::getMat_stitch1( int divnum, int maxpcnt, double *mObj )
{
	int ret=0;
	int p2cnt0, p2cnt1, p3cnt0, p3cnt1;
	double space = 10;
	double stx0[MAX_STP_CNT*2], sty0[MAX_STP_CNT*2], stz0[MAX_STP_CNT*2];
	double stx1[MAX_STP_CNT*2], sty1[MAX_STP_CNT*2], stz1[MAX_STP_CNT*2];
	double mat0[16], mat0i[16], mat1[16];

	if( lcrcnt<=1 && rcrcnt<=1 ){
		return -1;
	}

	//p2cnt0 = lcrs[1]->getSamplePoints2D( space, stx_cp[0], sty_cp[0] );
	//p2cnt1 = rcrs[1]->getSamplePoints2D( space, stx_cp[1], sty_cp[1] );
	p3cnt0 = lcrs[1]->getSamplePoints3D( space, stx[0], sty[0], stz[0], NULL );
	p3cnt1 = rcrs[1]->getSamplePoints3D( space, stx[1], sty[1], stz[1], NULL );
	stpcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
	if( maxpcnt>0 ){
		stpcnt = stpcnt < maxpcnt ? stpcnt : maxpcnt;
	}

	// approximate by 3D line -> evec[0][](left), evec[1][](right)
	double cen[3] = { 0.0, 0.0, 0.0 }, eval0[3], evec0[16], eval1[3], evec1[16];
	double tmpx,tmpy,tmpz;
	ret = PCA3( stx[0], sty[0], stz[0], p3cnt0, cen, eval0, evec0 );
	tmpx = stx[0][p3cnt0-1] - stx[0][0];
	tmpy = sty[0][p3cnt0-1] - sty[0][0];
	tmpz = stz[0][p3cnt0-1] - stz[0][0];
	normalize_v3( &tmpx, &tmpy, &tmpz ); 
	if( ret<0 ){
		evec0[0]=tmpx; evec0[1]=tmpy; evec0[2]=tmpz;
	} else if( evec0[0]*tmpx + evec0[1]*tmpy + evec0[2]*tmpz < 0 ){
		for( int i=0; i<16; i++ ){ evec0[i] = -evec0[i]; }
	}
	ret = PCA3( stx[1], sty[1], stz[1], p3cnt1, cen, eval1, evec1 );
	tmpx = stx[1][p3cnt1-1] - stx[1][0];
	tmpy = sty[1][p3cnt1-1] - sty[1][0];
	tmpz = stz[1][p3cnt1-1] - stz[1][0];
	normalize_v3( &tmpx, &tmpy, &tmpz ); 
	if( ret<0 ){
		evec1[0]=tmpx; evec1[1]=tmpy; evec1[2]=tmpz;
	} else if( evec1[0]*tmpx + evec1[1]*tmpy + evec1[2]*tmpz < 0 ){
		for( int i=0; i<16; i++ ){ evec1[i] = -evec1[i]; }
	}
	printf( "evec0: %.2f, %.2f, %.2f\n", evec0[0], evec0[1], evec0[2] );
	printf( "evec1: %.2f, %.2f, %.2f\n", evec1[0], evec1[1], evec1[2] );

	double cosa = evec0[0]*evec1[0] + evec0[1]*evec1[1] + evec0[2]*evec1[2];
	double cosap = cos( 2.0*M_PI/(double)divnum );
	printf( "cosa = %f\n", cosa );
	printf( "cosap = %f\n", cosap );
	if( cosa <= cosap ){
		unit_m44( mat0 );	// start from original pose
	} else {
		double sint = sqrt( (1.0-cosa)/(1.0-cosap) );
		double cost = sqrt( (cosa-cosap)/(1.0-cosap) );
		printf( "sint = %f\n", sint );
		printf( "cost = %f\n", cost );

		double vxO, vyO, vxA, vyA, vxB, vyB;
		double sinp = sin( M_PI/(double)divnum ), cosp = cos( M_PI/(double)divnum );
		double vx0[3], vy0[3], vz0[3], vx1[3], vy1[3], vz1[3];

		// original
		vx0[0] = vy0[0] = vz0[0] = 0.0;
		vx0[1] = evec0[0]; vy0[1] = evec0[1]; vz0[1] = evec0[2];
		vx0[2] = evec1[0]; vy0[2] = evec1[1]; vz0[2] = evec1[2];
		printf( "vx0: %.2f, %.2f, %.2f\n", vx0[0], vy0[0], vz0[0] );
		printf( "vx0: %.2f, %.2f, %.2f\n", vx0[1], vy0[1], vz0[1] );
		printf( "vx0: %.2f, %.2f, %.2f\n", vx0[2], vy0[2], vz0[2] );
		printf( "cos(vx0): %f\n", vx0[1]*vx0[2] + vy0[1]*vy0[2] + vz0[1]*vz0[2] );

		// target
		vx1[0] = vy1[0] = vz1[0] = 0.0;
		vxO=evec0[0]+evec1[0];	vyO=evec0[1]+evec1[1];	normalize_v2(&vxO,&vyO);
		vxA=cosp*vxO-sinp*vyO;	vyA=sinp*vxO+cosp*vyO;
		vxB=cosp*vxO+sinp*vyO;	vyB=-sinp*vxO+cosp*vyO;
		printf( "v*O: %.2f, %.2f\n", vxO, vyO );
		printf( "v*A: %.2f, %.2f\n", vxA, vyA );
		printf( "v*B: %.2f, %.2f\n", vxB, vyB );
		printf( "cos(AOB): %f\n", vxA * vxB + vyA * vyB );
		vx1[1] = sint*vxA;	vy1[1] = sint*vyA;	vz1[1] = -cost;
		vx1[2] = sint*vxB;	vy1[2] = sint*vyB;	vz1[2] = -cost;
		printf( "vx1: %.2f, %.2f, %.2f\n", vx1[0], vy1[0], vz1[0] );
		printf( "vx1: %.2f, %.2f, %.2f\n", vx1[1], vy1[1], vz1[1] );
		printf( "vx1: %.2f, %.2f, %.2f\n", vx1[2], vy1[2], vz1[2] );
		printf( "cos(vx1): %f\n", vx1[1]*vx1[2] + vy1[1]*vy1[2] + vz1[1]*vz1[2] );

		getMat( 3, vx0, vy0, vz0, vx1, vy1, vz1, mat0 );
		printf( "mat0: %.2f, %.2f, %.2f, %.2f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
		printf( "      %.2f, %.2f, %.2f, %.2f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
		printf( "      %.2f, %.2f, %.2f, %.2f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
		printf( "      %.2f, %.2f, %.2f, %.2f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);

		// rotate inv(mat0) -> new original points st*[0][], st*[1][]
		for( int j=0; j<2; j++ ){
			for( int i=0; i<stpcnt; i++ ){
				tmpx = mat0[0]*stx[j][i] + mat0[1]*sty[j][i] + mat0[2]*stz[j][i];
				tmpy = mat0[4]*stx[j][i] + mat0[5]*sty[j][i] + mat0[6]*stz[j][i];
				tmpz = mat0[8]*stx[j][i] + mat0[9]*sty[j][i] + mat0[10]*stz[j][i];
				stx[j][i] = tmpx;
				sty[j][i] = tmpy;
				stz[j][i] = tmpz;
			}
		}
	}

	// original
	for( int i=0; i<stpcnt; i++ ){
		int i2 = stpcnt+i;
		this->stx0[0][i] = stx0[i] = stx[0][i];
		this->sty0[0][i] = sty0[i] = sty[0][i];
		this->stz0[0][i] = stz0[i] = stz[0][i];
		this->stx0[1][i] = stx0[i2] = stx[1][i];
		this->sty0[1][i] = sty0[i2] = sty[1][i];
		this->stz0[1][i] = stz0[i2] = stz[1][i];
	}

	// rotated
	double axis[3]={0.0,0.0,1.0}, quat[4], rtmat0[16], rtmat1[16];
	axis_ang_quat( axis, 2.0*M_PI/(double)divnum, quat );	quat_mat( quat, rtmat0 );
	printf( "quat0: %.3f, %.3f, %.3f, %.3f\n", quat[0], quat[1], quat[2], quat[3] );
	printf( "norm = %f\n", sqrt( quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3] ));
	axis_ang_quat( axis, -2.0*M_PI/(double)divnum, quat );	quat_mat( quat, rtmat1 );
	printf( "quat1: %.3f, %.3f, %.3f, %.3f\n", quat[0], quat[1], quat[2], quat[3] );
	printf( "norm = %f\n", sqrt( quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3] ));
	printf( "rt0: %.3f, %.3f, %.3f, %.3f,\n", rtmat0[0], rtmat0[1], rtmat0[2], rtmat0[3]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[4], rtmat0[5], rtmat0[6], rtmat0[7]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[8], rtmat0[9], rtmat0[10], rtmat0[11]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[12], rtmat0[13], rtmat0[14], rtmat0[15]);
	printf( "rt1: %.3f, %.3f, %.3f, %.3f,\n", rtmat1[0], rtmat1[1], rtmat1[2], rtmat1[3]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[4], rtmat1[5], rtmat1[6], rtmat1[7]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[8], rtmat1[9], rtmat1[10], rtmat1[11]);
	printf( "     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[12], rtmat1[13], rtmat1[14], rtmat1[15]);

	for( int i=0; i<stpcnt; i++ ){
		int i2 = stpcnt+i;
		// right edge -> rotate to match left edge
		this->stx1[0][i] = stx1[i] = rtmat0[0]*stx[1][i] + rtmat0[4]*sty[1][i] + rtmat0[8]*stz[1][i] + rtmat0[12];
		this->sty1[0][i] = sty1[i] = rtmat0[1]*stx[1][i] + rtmat0[5]*sty[1][i] + rtmat0[9]*stz[1][i] + rtmat0[13];
		this->stz1[0][i] = stz1[i] = rtmat0[2]*stx[1][i] + rtmat0[6]*sty[1][i] + rtmat0[10]*stz[1][i] + rtmat0[14];
		// left edge -> rotate to match right edge
		this->stx1[1][i] = stx1[i2] = rtmat1[0]*stx[0][i] + rtmat1[4]*sty[0][i] + rtmat1[8]*stz[0][i] + rtmat1[12];
		this->sty1[1][i] = sty1[i2] = rtmat1[1]*stx[0][i] + rtmat1[5]*sty[0][i] + rtmat1[9]*stz[0][i] + rtmat1[13];
		this->stz1[1][i] = stz1[i2] = rtmat1[2]*stx[0][i] + rtmat1[6]*sty[0][i] + rtmat1[10]*stz[0][i] + rtmat1[14];
	}
	if(0){
		FILE *fp=fopen("output/dumpstp.csv", "w");
		for( int i=0; i<stpcnt*2; i++ ){
			fprintf( fp, "%d,%f,%f,%f,,%f,%f,%f\n",
				i, stx0[i],sty0[i],stz0[i],stx1[i],sty1[i],stz1[i]);
		}
		fclose(fp);
	}

	// targets
	for( int i=0; i<stpcnt*2; i++ ){
#if 0	// tartet points
		this->stxm[i/stpcnt][i%stpcnt] = stx1[i];
		this->stym[i/stpcnt][i%stpcnt] = sty1[i];
		this->stzm[i/stpcnt][i%stpcnt] = stz1[i];
#else	// mid points
		this->stxm[i/stpcnt][i%stpcnt] = stx1[i] = (stx0[i] + stx1[i]) *0.5;
		this->stym[i/stpcnt][i%stpcnt] = sty1[i] = (sty0[i] + sty1[i]) *0.5;
		this->stzm[i/stpcnt][i%stpcnt] = stz1[i] = (stz0[i] + stz1[i]) *0.5;
#endif
	}

	// 戻す（表示用）
	for( int j=0; j<2; j++ ){
		for( int i=0; i<stpcnt; i++ ){
			tmpx = mat0[0]*stx[j][i] + mat0[4]*sty[j][i] + mat0[8]*stz[j][i];
			tmpy = mat0[1]*stx[j][i] + mat0[5]*sty[j][i] + mat0[9]*stz[j][i];
			tmpz = mat0[2]*stx[j][i] + mat0[6]*sty[j][i] + mat0[10]*stz[j][i];
			stx[j][i] = tmpx;
			sty[j][i] = tmpy;
			stz[j][i] = tmpz;
		}
	}

	// getMat
	getMatRot( stpcnt*2, stx0, sty0, stz0, stx1, sty1, stz1, mat1 );
	printf( "m0 : %f, %f, %f, %f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
	printf( "     %f, %f, %f, %f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
	printf( "     %f, %f, %f, %f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
	printf( "     %f, %f, %f, %f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);
	printf( "m1 : %f, %f, %f, %f,\n", mat1[0], mat1[1], mat1[2], mat1[3]);
	printf( "     %f, %f, %f, %f,\n", mat1[4], mat1[5], mat1[6], mat1[7]);
	printf( "     %f, %f, %f, %f,\n", mat1[8], mat1[9], mat1[10], mat1[11]);
	printf( "     %f, %f, %f, %f,\n", mat1[12], mat1[13], mat1[14], mat1[15]);
	//mult_m44_n44( mat0, mat1, 1 ); // n44[16] = n44[16] * m44[16]
	mult_m44_n44( mat1, mat0 ); // m44[16] = n44[16] * m44[16] <- こっちが正しそう
	printf( "mat: %f, %f, %f, %f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
	printf( "     %f, %f, %f, %f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
	printf( "     %f, %f, %f, %f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
	printf( "     %f, %f, %f, %f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);

	mObj[0] = mat0[0];	mObj[1] = mat0[4];	mObj[2] = mat0[8];	mObj[3] = mat0[12];
	mObj[4] = mat0[1];	mObj[5] = mat0[5];	mObj[6] = mat0[9];	mObj[7] = mat0[13];
	mObj[8] = mat0[2];	mObj[9] = mat0[6];	mObj[10] = mat0[10];	mObj[11] = mat0[14];
	mObj[12] = mat0[3];	mObj[13] = mat0[7];	mObj[14] = mat0[11];	mObj[15] = mat0[15];

	return ret;
}

int papermodel::getStitch(int divnum, int maxpcnt, double* mObj)
{
	int ret = 0;
	int p2cnt0, p2cnt1, p3cnt0, p3cnt1;
	double space = 10, p3len0[MAX_STP_CNT], p3len1[MAX_STP_CNT];

	if (lcrcnt <= 1 && rcrcnt <= 1) {
		return -1;
	}

	//p2cnt0 = lcrs[1]->getSamplePoints2D( space, stx_cp[0], sty_cp[0] );
	//p2cnt1 = rcrs[1]->getSamplePoints2D( space, stx_cp[1], sty_cp[1] );
	//p3cnt0 = lcrs[1]->getSamplePoints3D(space, stx[0], sty[0], stz[0], NULL);
	//p3cnt1 = rcrs[1]->getSamplePoints3D(space, stx[1], sty[1], stz[1], NULL);
	p2cnt0 = p3cnt0 = lcrs[1]->getSamplePoints2D3D(space, stx[0], sty[0], stz[0], stx_cp[0], sty_cp[0], p3len0);
	p2cnt1 = p3cnt1 = rcrs[1]->getSamplePoints2D3D(space, stx[1], sty[1], stz[1], stx_cp[1], sty_cp[1], p3len1);
	stpcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
	if (maxpcnt > 0) {
		stpcnt = stpcnt < maxpcnt ? stpcnt : maxpcnt;
	}

	return ret;
}

int papermodel::getMat_stitch2(int divnum, double* mObj)
{
	int ret = 0;
	//int p2cnt0, p2cnt1, p3cnt0, p3cnt1;
	//double space = 10;
	double stx0[MAX_STP_CNT * 2], sty0[MAX_STP_CNT * 2], stz0[MAX_STP_CNT * 2];
	double stx1[MAX_STP_CNT * 2], sty1[MAX_STP_CNT * 2], stz1[MAX_STP_CNT * 2];
	double mat0[16], mat0i[16], mat1[16];

	if (stpcnt==0) {
		return -1;
	}
#if 0
	p2cnt0 = lcrs[1]->getSamplePoints2D(space, stx_cp[0], sty_cp[0]);
	p2cnt1 = rcrs[1]->getSamplePoints2D(space, stx_cp[1], sty_cp[1]);
	p3cnt0 = lcrs[1]->getSamplePoints3D(space, stx[0], sty[0], stz[0], NULL);
	p3cnt1 = rcrs[1]->getSamplePoints3D(space, stx[1], sty[1], stz[1], NULL);
	stpcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;
	if (maxpcnt > 0) {
		stpcnt = stpcnt < maxpcnt ? stpcnt : maxpcnt;
	}
#endif
	// approximate by 3D line -> evec[0][](left), evec[1][](right)
	double cen[3] = { 0.0, 0.0, 0.0 }, eval0[3], evec0[16], eval1[3], evec1[16];
	double tmpx, tmpy, tmpz;
	ret = PCA3(stx[0], sty[0], stz[0], stpcnt, cen, eval0, evec0);
	tmpx = stx[0][stpcnt - 1] - stx[0][0];
	tmpy = sty[0][stpcnt - 1] - sty[0][0];
	tmpz = stz[0][stpcnt - 1] - stz[0][0];
	normalize_v3(&tmpx, &tmpy, &tmpz);
	if (ret < 0) {
		evec0[0] = tmpx; evec0[1] = tmpy; evec0[2] = tmpz;
	}
	else if (evec0[0] * tmpx + evec0[1] * tmpy + evec0[2] * tmpz < 0) {
		for (int i = 0; i < 16; i++) { evec0[i] = -evec0[i]; }
	}
	ret = PCA3(stx[1], sty[1], stz[1], stpcnt, cen, eval1, evec1);
	tmpx = stx[1][stpcnt - 1] - stx[1][0];
	tmpy = sty[1][stpcnt - 1] - sty[1][0];
	tmpz = stz[1][stpcnt - 1] - stz[1][0];
	normalize_v3(&tmpx, &tmpy, &tmpz);
	if (ret < 0) {
		evec1[0] = tmpx; evec1[1] = tmpy; evec1[2] = tmpz;
	}
	else if (evec1[0] * tmpx + evec1[1] * tmpy + evec1[2] * tmpz < 0) {
		for (int i = 0; i < 16; i++) { evec1[i] = -evec1[i]; }
	}
	printf("evec0: %.2f, %.2f, %.2f\n", evec0[0], evec0[1], evec0[2]);
	printf("evec1: %.2f, %.2f, %.2f\n", evec1[0], evec1[1], evec1[2]);

	double cosa = evec0[0] * evec1[0] + evec0[1] * evec1[1] + evec0[2] * evec1[2];
	double cosap = cos(2.0 * M_PI / (double)divnum);
	printf("cosa = %f\n", cosa);
	printf("cosap = %f\n", cosap);
	if (cosa <= cosap) {
		unit_m44(mat0);	// start from original pose
	}
	else {
		double sint = sqrt((1.0 - cosa) / (1.0 - cosap));
		double cost = sqrt((cosa - cosap) / (1.0 - cosap));
		printf("sint = %f\n", sint);
		printf("cost = %f\n", cost);

		double vxO, vyO, vxA, vyA, vxB, vyB;
		double sinp = sin(M_PI / (double)divnum), cosp = cos(M_PI / (double)divnum);
		double vx0[3], vy0[3], vz0[3], vx1[3], vy1[3], vz1[3];

		// original
		vx0[0] = vy0[0] = vz0[0] = 0.0;
		vx0[1] = evec0[0]; vy0[1] = evec0[1]; vz0[1] = evec0[2];
		vx0[2] = evec1[0]; vy0[2] = evec1[1]; vz0[2] = evec1[2];
		printf("vx0: %.2f, %.2f, %.2f\n", vx0[0], vy0[0], vz0[0]);
		printf("vx0: %.2f, %.2f, %.2f\n", vx0[1], vy0[1], vz0[1]);
		printf("vx0: %.2f, %.2f, %.2f\n", vx0[2], vy0[2], vz0[2]);
		printf("cos(vx0): %f\n", vx0[1] * vx0[2] + vy0[1] * vy0[2] + vz0[1] * vz0[2]);

		// target
		vx1[0] = vy1[0] = vz1[0] = 0.0;
		vxO = evec0[0] + evec1[0];	vyO = evec0[1] + evec1[1];	normalize_v2(&vxO, &vyO);
		vxA = cosp * vxO - sinp * vyO;	vyA = sinp * vxO + cosp * vyO;
		vxB = cosp * vxO + sinp * vyO;	vyB = -sinp * vxO + cosp * vyO;
		printf("v*O: %.2f, %.2f\n", vxO, vyO);
		printf("v*A: %.2f, %.2f\n", vxA, vyA);
		printf("v*B: %.2f, %.2f\n", vxB, vyB);
		printf("cos(AOB): %f\n", vxA * vxB + vyA * vyB);
		vx1[1] = sint * vxA;	vy1[1] = sint * vyA;	vz1[1] = -cost;
		vx1[2] = sint * vxB;	vy1[2] = sint * vyB;	vz1[2] = -cost;
		printf("vx1: %.2f, %.2f, %.2f\n", vx1[0], vy1[0], vz1[0]);
		printf("vx1: %.2f, %.2f, %.2f\n", vx1[1], vy1[1], vz1[1]);
		printf("vx1: %.2f, %.2f, %.2f\n", vx1[2], vy1[2], vz1[2]);
		printf("cos(vx1): %f\n", vx1[1] * vx1[2] + vy1[1] * vy1[2] + vz1[1] * vz1[2]);

		getMat(3, vx0, vy0, vz0, vx1, vy1, vz1, mat0);
		printf("mat0: %.2f, %.2f, %.2f, %.2f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
		printf("      %.2f, %.2f, %.2f, %.2f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
		printf("      %.2f, %.2f, %.2f, %.2f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
		printf("      %.2f, %.2f, %.2f, %.2f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);

		// rotate inv(mat0) -> new original points st*[0][], st*[1][]
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < stpcnt; i++) {
				tmpx = mat0[0] * stx[j][i] + mat0[1] * sty[j][i] + mat0[2] * stz[j][i];
				tmpy = mat0[4] * stx[j][i] + mat0[5] * sty[j][i] + mat0[6] * stz[j][i];
				tmpz = mat0[8] * stx[j][i] + mat0[9] * sty[j][i] + mat0[10] * stz[j][i];
				stx[j][i] = tmpx;
				sty[j][i] = tmpy;
				stz[j][i] = tmpz;
			}
		}
	}

	// original
	for (int i = 0; i < stpcnt; i++) {
		int i2 = stpcnt + i;
		this->stx0[0][i] = stx0[i] = stx[0][i];
		this->sty0[0][i] = sty0[i] = sty[0][i];
		this->stz0[0][i] = stz0[i] = stz[0][i];
		this->stx0[1][i] = stx0[i2] = stx[1][i];
		this->sty0[1][i] = sty0[i2] = sty[1][i];
		this->stz0[1][i] = stz0[i2] = stz[1][i];
	}

	// rotated
	double axis[3] = { 0.0,0.0,1.0 }, quat[4], rtmat0[16], rtmat1[16];
	axis_ang_quat(axis, 2.0 * M_PI / (double)divnum, quat);	quat_mat(quat, rtmat0);
	printf("quat0: %.3f, %.3f, %.3f, %.3f\n", quat[0], quat[1], quat[2], quat[3]);
	printf("norm = %f\n", sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]));
	axis_ang_quat(axis, -2.0 * M_PI / (double)divnum, quat);	quat_mat(quat, rtmat1);
	printf("quat1: %.3f, %.3f, %.3f, %.3f\n", quat[0], quat[1], quat[2], quat[3]);
	printf("norm = %f\n", sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]));
	printf("rt0: %.3f, %.3f, %.3f, %.3f,\n", rtmat0[0], rtmat0[1], rtmat0[2], rtmat0[3]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[4], rtmat0[5], rtmat0[6], rtmat0[7]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[8], rtmat0[9], rtmat0[10], rtmat0[11]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat0[12], rtmat0[13], rtmat0[14], rtmat0[15]);
	printf("rt1: %.3f, %.3f, %.3f, %.3f,\n", rtmat1[0], rtmat1[1], rtmat1[2], rtmat1[3]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[4], rtmat1[5], rtmat1[6], rtmat1[7]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[8], rtmat1[9], rtmat1[10], rtmat1[11]);
	printf("     %.3f, %.3f, %.3f, %.3f,\n", rtmat1[12], rtmat1[13], rtmat1[14], rtmat1[15]);

	for (int i = 0; i < stpcnt; i++) {
		int i2 = stpcnt + i;
		// right edge -> rotate to match left edge
		this->stx1[0][i] = stx1[i] = rtmat0[0] * stx[1][i] + rtmat0[4] * sty[1][i] + rtmat0[8] * stz[1][i] + rtmat0[12];
		this->sty1[0][i] = sty1[i] = rtmat0[1] * stx[1][i] + rtmat0[5] * sty[1][i] + rtmat0[9] * stz[1][i] + rtmat0[13];
		this->stz1[0][i] = stz1[i] = rtmat0[2] * stx[1][i] + rtmat0[6] * sty[1][i] + rtmat0[10] * stz[1][i] + rtmat0[14];
		// left edge -> rotate to match right edge
		this->stx1[1][i] = stx1[i2] = rtmat1[0] * stx[0][i] + rtmat1[4] * sty[0][i] + rtmat1[8] * stz[0][i] + rtmat1[12];
		this->sty1[1][i] = sty1[i2] = rtmat1[1] * stx[0][i] + rtmat1[5] * sty[0][i] + rtmat1[9] * stz[0][i] + rtmat1[13];
		this->stz1[1][i] = stz1[i2] = rtmat1[2] * stx[0][i] + rtmat1[6] * sty[0][i] + rtmat1[10] * stz[0][i] + rtmat1[14];
	}
	if (0) {
		FILE* fp = fopen("output/dumpstp.csv", "w");
		for (int i = 0; i < stpcnt * 2; i++) {
			fprintf(fp, "%d,%f,%f,%f,,%f,%f,%f\n",
				i, stx0[i], sty0[i], stz0[i], stx1[i], sty1[i], stz1[i]);
		}
		fclose(fp);
	}

	// targets
	for (int i = 0; i < stpcnt * 2; i++) {
#if 0	// tartet points
		this->stxm[i / stpcnt][i % stpcnt] = stx1[i];
		this->stym[i / stpcnt][i % stpcnt] = sty1[i];
		this->stzm[i / stpcnt][i % stpcnt] = stz1[i];
#else	// mid points
		this->stxm[i / stpcnt][i % stpcnt] = stx1[i] = (stx0[i] + stx1[i]) * 0.5;
		this->stym[i / stpcnt][i % stpcnt] = sty1[i] = (sty0[i] + sty1[i]) * 0.5;
		this->stzm[i / stpcnt][i % stpcnt] = stz1[i] = (stz0[i] + stz1[i]) * 0.5;
#endif
	}

	// 戻す（表示用）
	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < stpcnt; i++) {
			tmpx = mat0[0] * stx[j][i] + mat0[4] * sty[j][i] + mat0[8] * stz[j][i];
			tmpy = mat0[1] * stx[j][i] + mat0[5] * sty[j][i] + mat0[9] * stz[j][i];
			tmpz = mat0[2] * stx[j][i] + mat0[6] * sty[j][i] + mat0[10] * stz[j][i];
			stx[j][i] = tmpx;
			sty[j][i] = tmpy;
			stz[j][i] = tmpz;
		}
	}

	// getMat
	getMatRot(stpcnt * 2, stx0, sty0, stz0, stx1, sty1, stz1, mat1);
	printf("m0 : %f, %f, %f, %f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
	printf("     %f, %f, %f, %f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
	printf("     %f, %f, %f, %f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
	printf("     %f, %f, %f, %f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);
	printf("m1 : %f, %f, %f, %f,\n", mat1[0], mat1[1], mat1[2], mat1[3]);
	printf("     %f, %f, %f, %f,\n", mat1[4], mat1[5], mat1[6], mat1[7]);
	printf("     %f, %f, %f, %f,\n", mat1[8], mat1[9], mat1[10], mat1[11]);
	printf("     %f, %f, %f, %f,\n", mat1[12], mat1[13], mat1[14], mat1[15]);
	//mult_m44_n44( mat0, mat1, 1 ); // n44[16] = n44[16] * m44[16]
	mult_m44_n44(mat1, mat0); // m44[16] = n44[16] * m44[16] <- こっちが正しそう
	printf("mat: %f, %f, %f, %f,\n", mat0[0], mat0[1], mat0[2], mat0[3]);
	printf("     %f, %f, %f, %f,\n", mat0[4], mat0[5], mat0[6], mat0[7]);
	printf("     %f, %f, %f, %f,\n", mat0[8], mat0[9], mat0[10], mat0[11]);
	printf("     %f, %f, %f, %f,\n", mat0[12], mat0[13], mat0[14], mat0[15]);

	mObj[0] = mat0[0];	mObj[1] = mat0[4];	mObj[2] = mat0[8];	mObj[3] = mat0[12];
	mObj[4] = mat0[1];	mObj[5] = mat0[5];	mObj[6] = mat0[9];	mObj[7] = mat0[13];
	mObj[8] = mat0[2];	mObj[9] = mat0[6];	mObj[10] = mat0[10];	mObj[11] = mat0[14];
	mObj[12] = mat0[3];	mObj[13] = mat0[7];	mObj[14] = mat0[11];	mObj[15] = mat0[15];

	return ret;
}

