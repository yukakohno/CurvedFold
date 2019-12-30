#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"

// a) 3D曲線と折り角度から、折り線を求める
int crease::calcXA_CP( int flg_interpolate, rectify_params *rp )
{
	if( flg_interpolate ){
		initX();	Xcnt = XcntOrg;
		interpolate_spline( Px, Pcnt, kv, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
		interpolate_spline( Pa, Pcnt, alpha, Xcnt );
	}
	if( rp->flg_rectifyT ){ rectifyTauBezier2( &rp->rectT ); }
	calcXTNB( m3 );				// kv,tr -> X -> calcTNB();
	if( rp->flg_rectifyA ){ rectifyAlphaBezier( &rp->rectA ); }
	calcAK_K2D( );				// alpha, kv -> sina, cosa, tana, da, k2d
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	calcRuling( rp->flg_rectifyR, rp->rectifyR_kvthres );

	//setP_k(-1);
	//setP_t(-1);
	setP_k2(-1);
	//setP_a(-1);

	return 0;
}

// b) 折り線と折り角度から、3D曲線を求める
int crease::calcCPA_X( int flg_interpolate, rectify_params *rp )
{
	if( flg_interpolate ){
		initX();	Xcnt = XcntOrg;
		interpolate_spline( Px2d, Pcnt, k2d, Xcnt );
		interpolate_spline( Pa, Pcnt, alpha, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
	}
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	if( rp->flg_rectifyA ){ rectifyAlphaBezier( &rp->rectA ); }
	calcAK2D_K();				// alpha, k2d -> kv
	if( rp->flg_rectifyT ){ rectifyTauBezier2( &rp->rectT ); }
	calcXTNB( m3 );				// kv,tr -> X -> calcTNB();
	calcDA();
	calcRuling( rp->flg_rectifyR, rp->rectifyR_kvthres );

	setP_k(-1);
	//setP_t(-1);
	//setP_k2(-1);
	//setP_a(-1);

	return 0;
}


// c) 折り線と3D曲線から、折り角度を求める
int crease::calcCPX_A( int flg_interpolate, rectify_params *rp )
{
	if( flg_interpolate ){
		initX();	Xcnt = XcntOrg;
		interpolate_spline( Px, Pcnt, kv, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
		interpolate_spline( Px2d, Pcnt, k2d, Xcnt );
	}
	if( rp->flg_rectifyT ){ rectifyTauBezier2( &rp->rectT ); }
	calcXTNB( m3 );				// kv,tr -> X -> calcTNB();
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	calcAlpha( rp->flg_rectifyA, &rp->rectA );	// kv,k2d -> (rectifyAlpha) -> alpha, sina, cosa, tana, da
	calcRuling( rp->flg_rectifyR, rp->rectifyR_kvthres );

	//setP_k(-1);
	//setP_t(-1);
	//setP_k2(-1);
	setP_a(-1);

	return 0;
}
