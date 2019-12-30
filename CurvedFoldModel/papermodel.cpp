#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "papermodel.h"
#include "Bspline.h"

papermodel::papermodel()
{
	init();
	int Scnt = 10000;
	ht1 = new double[Scnt];
	hSx = new double[Scnt];
	hSy = new double[Scnt];
	hL = new double[Scnt];

	Bspline bs;
	int ret = bs.calcT( Scnt-1, 3, 0, ht1); ht1[Scnt-1]=3.0;
}

papermodel::~papermodel()
{
	delete[] ht1;
	delete[] hSx;
	delete[] hSy;
	delete[] hL;
}

void papermodel::init()
{
	flg_postproc_type = PPTYPE_UNDEF;
	flg_interpolate = 1;
	flg_rectifyA = 0;
	flg_rectifyT = 0;
	flg_rectifyR= 0;

	pw = PPWIDTH;
	ph = PPHEIGHT;
	psx = 0;
	psy = 0;
	pex = pw;
	pey = ph;

	crcnt = lcrcnt = rcrcnt = 0;
	for( int i=0; i<MAX_CRS_CNT; i++ ){
		crs[i].init();
	}
	memset(lcrs, NULL, sizeof(crease*)*MAX_CRS_CNT);
	memset(rcrs, NULL, sizeof(crease*)*MAX_CRS_CNT);

	for( int i=0; i<dccnt; i++ ){
		dcurve[i].init();
	}
	for( int i=0; i<fccnt; i++ ){
		fcurve[i].init();
	}
	for( int i=0; i<tccnt; i++ ){
		tcurve[i].init();
	}

	flg_rectifyCrs = 0;
	hCrs_evaltype = EVTYPE_RULCROSS;
	hCrs_cridx = 0;
	hCrs_cpidx = -1;

	initHCrs();
	// これ以下になったら最適化終了, コスト増加の許容値, cverr[]にかける表示用
	// tau
	hCrs_thresdiff[EVTYPE_TORSION] = 0.002;
	hCrs_mgndiff[EVTYPE_TORSION] = 0.0001;
	hCrs_cvmaxerr[EVTYPE_TORSION] = 0.05;
	// betal-betar
	hCrs_thresdiff[EVTYPE_BETALR] = 5.0/180.0*M_PI;
	hCrs_mgndiff[EVTYPE_BETALR] = 3.0/180.0*M_PI;
	hCrs_cvmaxerr[EVTYPE_BETALR] = M_PI;
	// ruling crossing area
	hCrs_thresdiff[EVTYPE_RULCROSS] = 0.1;
	hCrs_mgndiff[EVTYPE_RULCROSS] = 0.005;
	hCrs_cvmaxerr[EVTYPE_RULCROSS] = 100.0;

	initX();
}

void papermodel::initX()
{
	plcnt = plecnt = 0;
	memset( plvcnt, 0, sizeof(double)*MAX_PL_FCNT );
	memset( plx, 0, sizeof(double)*MAX_PL_FCNT*4 );
	memset( ply, 0, sizeof(double)*MAX_PL_FCNT*4 );
	memset( plz, 0, sizeof(double)*MAX_PL_FCNT*4 );
	memset( plx_cp, 0, sizeof(double)*MAX_PL_FCNT*4 );
	memset( ply_cp, 0, sizeof(double)*MAX_PL_FCNT*4 );
	memset( plmat, 0, sizeof(double)*MAX_PL_FCNT*16 );
	memset( plex, 0, sizeof(double)*MAX_PL_ECNT*2 );
	memset( pley, 0, sizeof(double)*MAX_PL_ECNT*2 );
	memset( plez, 0, sizeof(double)*MAX_PL_ECNT*2 );

	o_vcnt = o_fcnt = 0;
	memset( o_vx, 0, sizeof(double)*MAX_PL_FCNT*4);
	memset( o_vy, 0, sizeof(double)*MAX_PL_FCNT*4);
	memset( o_vz, 0, sizeof(double)*MAX_PL_FCNT*4);
	memset( o_fvcnt, 0, sizeof(int)*MAX_PL_FCNT);
	memset( o_fi, 0, sizeof(int)*MAX_PL_FCNT*4);

	memset( cverr, 0, sizeof(double)*MAX_SPCNT );
	memset( cvminerr, 0, sizeof(double)*MAX_SPCNT );
}

void papermodel::initHCrs()
{
	hCrs_cpdx = hCrs_cpdy = 0;
	for( int i=0; i<MAX_CRS_OTYPE; i++ ){
		hCrs_mindiff[i] = 10000.0;	// （暫定）最小値
	}
	memset(cverr, 0, sizeof(double)*MAX_SPCNT);
	memset(cvminerr, 0, sizeof(double)*MAX_SPCNT);
}


void papermodel::set_postproc_type( postproc_type ppt )
{
	if( flg_postproc_type != PPTYPE_X )
	{
		flg_postproc_type = ppt;
	}
}

