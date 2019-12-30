#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "papermodel.h"

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
			loadcv0( "input/cv0.txt" );
			loadBsplineCP( "input/bspline.txt" );
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
		}
		for( int j=1; j<rcrcnt; j++ ){
			crease *c = rcrs[j];
			if( c->flg_org==0 ){ // 0:org* 1:CP*
				c->org2CP();
				c->flg_org=1;
			}
		}

		addCrease( 1, 1, -1.0, 0.0 );	// curve fold -> papermodel
	}
	else if( flg_postproc_type==PPTYPE_MOVECP0 || flg_postproc_type==PPTYPE_MOVECP1 )
	{
		if( hCrs_cridx<0 ){
			crs[-hCrs_cridx].init_rlflg();
			for( int c= -hCrs_cridx+1; c<lcrcnt; c++ ){
				crs[c].init_rlflg();
				crs[c].init_rrflg();
			}
			addCrease( -hCrs_cridx, 0, -1.0, 0.0 ); // 左側 -hCrs_cridx番目以降のみ更新
		} else if( hCrs_cridx>0 ){
			crs[hCrs_cridx].init_rrflg();
			for( int c= hCrs_cridx+1; c<rcrcnt; c++ ){
				crs[c].init_rlflg();
				crs[c].init_rrflg();
			}
			addCrease( 0, hCrs_cridx, -1.0, 0.0 ); // 右側 hCrs_cridx番目以降のみ更新
		}
	}
	else if( flg_postproc_type==PPTYPE_OPTIMIZE )
	{
		if( flg_rectifyCrs<0 ){
			addCrease( -flg_rectifyCrs, 0, -1.0, 0.0 ); // 左側 -hCrs_cridx番目以降のみ更新
		} else if( flg_rectifyCrs>0 ){
			addCrease( 0, flg_rectifyCrs, -1.0, 0.0 ); // 右側 hCrs_cridx番目以降のみ更新
		}

	}

	for( int c=0; c<crcnt; c++ )
	{
		crs[c].setEnds( psx,psy,pex,pey, dccnt,dcurve );
		crs[c].calcRLenP(psx,psy,pex,pey);		// ruling長さ修正（枠線、曲線までの長さ優先）
		crs[c].calcRLenC( dccnt, dcurve );		// ruling長さ修正（トリム）
		if( flg_rectifyR ){
			crs[c].calcRLenHoseiR();	// 曲率小さい部分は長さ0に
		}
	}

	calcRulingPly();	// TODO: 2D,3Dポリゴン、変換行列作成
	return 0;
}
#if 0
int papermodel::resetpostproc( int _cidx )
{
	if( 0 < -_cidx && -_cidx < lcrcnt ){
		int idx = lcrs[-_cidx]->org_c0_idx;
		cvcnt0[idx] = cvcnt1[idx] = 0;
		lcrs[-_cidx]->init();
		// 詰める
		for( int i=-_cidx; i<lcrcnt-1; i++ ){
			lcrs[i] = lcrs[i+1];
		}
		lcrcnt--;
		for( int i=0; i<ccnt1; i++ ){
			if( ctype1[i]==CTYPE_FOLD_LEFT && cridx1[i] >= -_cidx ){
				cridx1[i]--;
			}
		}
	} else if( 0 < _cidx && _cidx < rcrcnt ){
		int idx = rcrs[_cidx]->org_c0_idx;
		cvcnt0[idx] = cvcnt1[idx] = 0;
		rcrs[_cidx]->init();
		// 詰める
		for( int i=_cidx; i<rcrcnt-1; i++ ){
			rcrs[i] = rcrs[i+1];
		}
		rcrcnt--;
		for( int i=0; i<ccnt1; i++ ){
			if( ctype1[i]==CTYPE_FOLD_RIGHT && cridx1[i] >= _cidx ){
				cridx1[i]--;
			}
		}
	} else {
		initC0();
	}

#if 1
	postproc();
#else
	//dumpcv0( "cv0.txt" );
	//dumpBsplineCP( "bspline.txt" );
	cv0to1();
	calcCurveEnd2();	// 折り線の始点終点を指定
	for( int i=0; i<crcnt; i++ ){
		crs[i].flg_org=0;
	}
	cv1toCreaseOrg();	// curve fold -> ruling
	addCrease(-1.0,0.0);	// curve fold -> ruling
	calcRLenP();		// 枠線 → ruling長さ修正
	calcRLenC();		// curve trim → ruling長さ修正
	calcRulingPly();	// TODO: 2D,3Dポリゴン、変換行列作成
#endif
	return 0;
}
#endif



