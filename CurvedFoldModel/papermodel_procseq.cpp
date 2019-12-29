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
		crcnt=1;
		return -1;
	}

	if( flg_postproc_type==PPTYPE_OPEN
		|| flg_postproc_type==PPTYPE_PRICURVE || flg_postproc_type==PPTYPE_RTCURVE
		|| flg_postproc_type==PPTYPE_ADDCREASE || flg_postproc_type==PPTYPE_FMOT )
	{
		crs[0].Xsidx = CEMGN;
		crs[0].Xeidx = crs[0].Xcnt-CEMGN-1;

		if( flg_postproc_type==PPTYPE_OPEN ){
			crcnt = 1;
			loadcv0( "input/cv0.txt" );
		}

		// generate fcurve/tcurve from dcurve
		tccnt = 0;
		for( int i=0; i<dccnt; i++ ){
			if( dcurve[i].ctype == CTYPE_TRIM )
			{
				//tcurve[dccnt].create_drawcurve( &(dcurve[i]), &(crs[0]), psx,psy,pex,pey );
				dcurve[i].closecv( 5.0 );		// loop draw curve, unit: mm
				tcurve[tccnt].closecv( 5.0 );	// loop trim curve, unit: mm
				//tccnt++;
			}
		}
	}

	for( int c=0; c<crcnt; c++ )
	{
		crs[c].setEnds( psx,psy,pex,pey, dccnt,dcurve );
		crs[c].calcRLenP(psx,psy,pex,pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj
		crs[c].calcRLenC( dccnt, dcurve );		// ruling’·‚³C³iƒgƒŠƒ€j
		if( flg_rectifyR ){
			crs[c].calcRLenHoseiR();	// ‹È—¦¬‚³‚¢•”•ª‚Í’·‚³0‚É
		}
	}

	calcRulingPly();	// TODO: 2D,3Dƒ|ƒŠƒSƒ“A•ÏŠ·s—ñì¬
	return 0;
}



