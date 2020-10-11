#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "papermodel.h"
#include "util.h"

int set3D( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double x3_cp, double y3_cp,
		  double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		  double x2_3d, double y2_3d, double z2_3d, double *x3_3d, double *y3_3d, double *z3_3d )
{
	int ret=0;
	double x2[3], y2[3], z2[3], x3[3], y3[3], z3[3], mat[16];
	x2[0] = x0_cp;	y2[0] = y0_cp;	z2[0] = 0.0;
	x2[1] = x1_cp;	y2[1] = y1_cp;	z2[1] = 0.0;
	x2[2] = x2_cp;	y2[2] = y2_cp;	z2[2] = 0.0;
	x3[0] = x0_3d;	y3[0] = y0_3d;	z3[0] = z0_3d;
	x3[1] = x1_3d;	y3[1] = y1_3d;	z3[1] = z1_3d;
	x3[2] = x2_3d;	y3[2] = y2_3d;	z3[2] = z2_3d;
	ret = getMat( 3, x2,y2,z2, x3,y3,z3, mat);
	*x3_3d = mat[0]*x3_cp + mat[1]*y3_cp + mat[2]*0.0 + mat[3];
	*y3_3d = mat[4]*x3_cp + mat[5]*y3_cp + mat[6]*0.0 + mat[7];
	*z3_3d = mat[8]*x3_cp + mat[9]*y3_cp + mat[10]*0.0 + mat[11];
	return ret;
}

int papermodel::postproc()
{
	if(flg_postproc_type == PPTYPE_UNDEF || flg_postproc_type == PPTYPE_X )
	{
		crs[0].initEnds();
		lcrs[0]=rcrs[0]=&(crs[0]);
		crcnt=lcrcnt=rcrcnt=1;
		return -1;
	}

	if( flg_postproc_type==PPTYPE_OPEN
		|| flg_postproc_type==PPTYPE_PRICURVE || flg_postproc_type==PPTYPE_RTCURVE
		|| flg_postproc_type==PPTYPE_ADDCREASE || flg_postproc_type==PPTYPE_FMOT )
	{
		crs[0].Xsidx = CEMGN;
		crs[0].Xeidx = crs[0].Xcnt-CEMGN-1;

		if( flg_postproc_type==PPTYPE_OPEN ){
			lcrs[0] = rcrs[0] = &(crs[0]);
			crcnt = lcrcnt = rcrcnt = 1;
		}

		CP2DC();	// generate dcurve from CP, add dcurve if neccessary, overwrite if flg_org==1

		// generate fcurve/tcurve from dcurve
		fccnt = tccnt = 0;
		for( int i=0; i<dccnt; i++ ){
			if( dcurve[i].ctype == CTYPE_FOLD || dcurve[i].ctype == CTYPE_FOLD_LEFT || dcurve[i].ctype == CTYPE_FOLD_RIGHT )
			{	// 左右判定のため（だけ）に使う
				fcurve[fccnt].create_drawcurve( &(dcurve[i]), &(crs[0]), psx,psy,pex,pey );
				fccnt++;
			}
			if( dcurve[i].ctype == CTYPE_TRIM )
			{
				//tcurve[dccnt].create_drawcurve( &(dcurve[i]), &(crs[0]), psx,psy,pex,pey );
				dcurve[i].closecv( 5.0 );		// loop draw curve, unit: mm
				tcurve[tccnt].closecv( 5.0 );	// loop trim curve, unit: mm
				//tccnt++;
			}
		}

		// dcurve をもとに crease の順番入れ替え、CP以外のcrease情報リセット
		FC2Crs();

		// FC -> CP if flg_org==0
		for( int j=1; j<lcrcnt; j++ ){
			crease *c = lcrs[j];
			if( c->flg_org==0 ){ // 0:org* 1:CP*
				c->org2CP();
				c->flg_org=1;
			}
			c->flg_xspc_def=0;
		}
		for( int j=1; j<rcrcnt; j++ ){
			crease *c = rcrs[j];
			if( c->flg_org==0 ){ // 0:org* 1:CP*
				c->org2CP();
				c->flg_org=1;
			}
			c->flg_xspc_def=0;
		}

		addCrease( 1, 1, -1.0, 0.0 );	// curve fold -> papermodel
	}

	for( int c=0; c<crcnt; c++ )
	{
		crs[c].setEnds( psx,psy,pex,pey, dccnt,dcurve );
		crs[c].calcRLenP(psx,psy,pex,pey);		// ruling長さ修正（枠線、曲線までの長さ優先）
		crs[c].calcRLenC( dccnt, dcurve );		// ruling長さ修正（トリム）
		if( rp.flg_rectifyR ){
			crs[c].calcRLenHoseiR( rp.rectifyR_kvthres );	// 曲率小さい部分は長さ0に
		}
	}

	// c1-> Xxs,Xys,Xzs, Xxe,Xye,Xze 算出
	for( int i=1; i<lcrcnt; i++ ){
		crease *c0=lcrs[i-1];
		crease *c1=lcrs[i];
		double c0x_cp, c0y_cp, c0x,c0y,c0z;
#if 1
		if( c1->Xsidx-1 < c0->Xsidx ){
			//printf( "lcrs[%d] c0*=X*s\n", i );
			c0x_cp = c0->Xxs2d;	c0y_cp = c0->Xys2d;
			c0x = c0->Xxs;	c0y = c0->Xys;	c0z = c0->Xzs;
		} else {
			//printf( "lcrs[%d] c0*=X*[(c1->Xsidx-1)=%d]\n", i, c1->Xsidx-1 );
			c0x_cp = c0->Xx2d[c1->Xsidx-1];	c0y_cp = c0->Xy2d[c1->Xsidx-1];
			c0x = c0->Xx[c1->Xsidx-1];	c0y = c0->Xy[c1->Xsidx-1];	c0z = c0->Xz[c1->Xsidx-1];
		}
		set3D( c1->Xx2d[c1->Xsidx], c1->Xy2d[c1->Xsidx], c0->Xx2d[c1->Xsidx], c0->Xy2d[c1->Xsidx],
			c0x_cp, c0y_cp, c1->Xxs2d, c1->Xys2d,
			c1->Xx[c1->Xsidx], c1->Xy[c1->Xsidx], c1->Xz[c1->Xsidx],
			c0->Xx[c1->Xsidx], c0->Xy[c1->Xsidx], c0->Xz[c1->Xsidx],
			c0x, c0y, c0z, &(c1->Xxs), &(c1->Xys), &(c1->Xzs) );
#endif
#if 1
		if( c1->Xeidx+1 > c0->Xeidx ){
			//printf( "lcrs[%d] c0*=X*e\n", i );
			c0x_cp = c0->Xxe2d;	c0y_cp = c0->Xye2d;
			c0x = c0->Xxe;	c0y = c0->Xye;	c0z = c0->Xze;
		} else {
			//printf( "lcrs[%d] c0*=X*[(c1->Xsidx-1)=%d]\n", i, c1->Xeidx+1 );
			c0x_cp = c0->Xx2d[c1->Xeidx+1];	c0y_cp = c0->Xy2d[c1->Xeidx+1];
			c0x = c0->Xx[c1->Xeidx+1];	c0y = c0->Xy[c1->Xeidx+1];	c0z = c0->Xz[c1->Xeidx+1];
		}
		set3D( c1->Xx2d[c1->Xeidx], c1->Xy2d[c1->Xeidx],
			c0->Xx2d[c1->Xeidx], c0->Xy2d[c1->Xeidx],
			c0x_cp, c0y_cp, c1->Xxe2d, c1->Xye2d,
			c1->Xx[c1->Xeidx], c1->Xy[c1->Xeidx], c1->Xz[c1->Xeidx],
			c0->Xx[c1->Xeidx], c0->Xy[c1->Xeidx], c0->Xz[c1->Xeidx],
			c0x, c0y, c0z, &(c1->Xxe), &(c1->Xye), &(c1->Xze) );
#endif
	}
	for( int i=1; i<rcrcnt; i++ ){
		crease *c0=rcrs[i-1];
		crease *c1=rcrs[i];
		double c0x_cp, c0y_cp, c0x,c0y,c0z;
#if 1
		if( c1->Xsidx-1 < c0->Xsidx ){
			//printf( "rcrs[%d] c0*=X*s\n", i );
			c0x_cp = c0->Xxs2d;	c0y_cp = c0->Xys2d;
			c0x = c0->Xxs;	c0y = c0->Xys;	c0z = c0->Xzs;
		} else {
			//printf( "rcrs[%d] c0*=X*[(c1->Xsidx-1)=%d]\n", i, c1->Xsidx-1 );
			c0x_cp = c0->Xx2d[c1->Xsidx-1];	c0y_cp = c0->Xy2d[c1->Xsidx-1];
			c0x = c0->Xx[c1->Xsidx-1];	c0y = c0->Xy[c1->Xsidx-1];	c0z = c0->Xz[c1->Xsidx-1];
		}
		set3D( c1->Xx2d[c1->Xsidx], c1->Xy2d[c1->Xsidx], c0->Xx2d[c1->Xsidx], c0->Xy2d[c1->Xsidx],
			c0x_cp, c0y_cp, c1->Xxs2d, c1->Xys2d,
			c1->Xx[c1->Xsidx], c1->Xy[c1->Xsidx], c1->Xz[c1->Xsidx],
			c0->Xx[c1->Xsidx], c0->Xy[c1->Xsidx], c0->Xz[c1->Xsidx],
			c0x, c0y, c0z, &(c1->Xxs), &(c1->Xys), &(c1->Xzs) );
#endif
#if 1
		if( c1->Xeidx+1 > c0->Xeidx ){
			//printf( "rcrs[%d] c0*=X*e\n", i );
			c0x_cp = c0->Xxe2d;	c0y_cp = c0->Xye2d;
			c0x = c0->Xxe;	c0y = c0->Xye;	c0z = c0->Xze;
		} else {
			//printf( "rcrs[%d] c0*=X*[(c1->Xsidx-1)=%d]\n", i, c1->Xeidx+1 );
			c0x_cp = c0->Xx2d[c1->Xeidx+1];	c0y_cp = c0->Xy2d[c1->Xeidx+1];
			c0x = c0->Xx[c1->Xeidx+1];	c0y = c0->Xy[c1->Xeidx+1];	c0z = c0->Xz[c1->Xeidx+1];
		}
		set3D( c1->Xx2d[c1->Xeidx], c1->Xy2d[c1->Xeidx], c0->Xx2d[c1->Xeidx], c0->Xy2d[c1->Xeidx],
			c0x_cp, c0y_cp, c1->Xxe2d, c1->Xye2d,
			c1->Xx[c1->Xeidx], c1->Xy[c1->Xeidx], c1->Xz[c1->Xeidx],
			c0->Xx[c1->Xeidx], c0->Xy[c1->Xeidx], c0->Xz[c1->Xeidx],
			c0x, c0y, c0z, &(c1->Xxe), &(c1->Xye), &(c1->Xze) );
#endif
	}

	calcRulingPly();	// TODO: 2D,3Dポリゴン、変換行列作成
	getTgt2D3D();

	return 0;
}
