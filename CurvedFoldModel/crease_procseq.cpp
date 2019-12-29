#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"
#include "Spline.h"

// a) 3D‹Èü‚ÆÜ‚èŠp“x‚©‚çAÜ‚èü‚ð‹‚ß‚é
int crease::calcXA_CP( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR )
{
	if( flg_interpolate ){
		initX();	Xcnt = XCNT;
		interpolate_spline( Px, Pcnt, kv, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
		interpolate_spline( Pa, Pcnt, alpha, Xcnt );
	}
	flg_src_s1e1 = 1;			// 0:k2d, 1:kv*cos(alpha)
	calcXTNB( flg_rectifyT, m3 );	// kv,tr -> (rectifyTau) -> X -> calcTNB();
	calcAK_K2D( flg_rectifyA );	// alpha, kv -> (rectifyAlpha) -> sina, cosa, tana, da, k2d
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	flg_src_s1e1 = 0;
	calcRuling( flg_rectifyR );

	//setP_k(-1);
	//setP_t(-1);
	setP_k2(-1);
	//setP_a(-1);

	return 0;
}

// b) Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
int crease::calcCPA_X( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR )
{
	if( flg_interpolate ){
		initX();	Xcnt = XCNT;
		interpolate_spline( Px2d, Pcnt, k2d, Xcnt );
		interpolate_spline( Pa, Pcnt, alpha, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
	}
	flg_src_s1e1 = 0;			// 0:k2d, 1:kv*cos(alpha)
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	calcAK2D_K( flg_rectifyA );	// alpha, k2d -> (rectifyAlpha) -> kv
	calcXTNB( flg_rectifyT, m3 );	// kv,tr -> (rectifyTau) -> X -> calcTNB();
	calcDA();
	calcRuling( flg_rectifyR );

	setP_k(-1);
	//setP_t(-1);
	//setP_k2(-1);
	//setP_a(-1);

	return 0;
}


// c) Ü‚èü‚Æ3D‹Èü‚©‚çAÜ‚èŠp“x‚ð‹‚ß‚é
int crease::calcCPX_A( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR )
{
	if( flg_interpolate ){
		initX();	Xcnt = XCNT;
		interpolate_spline( Px, Pcnt, kv, Xcnt );
		interpolate_spline( Py, Pcnt, tr, Xcnt );
		interpolate_spline( Px2d, Pcnt, k2d, Xcnt );
	}
	flg_src_s1e1 = 0;			// 0:k2d, 1:kv*cos(alpha)
	calcXTNB( flg_rectifyT, m3 );	// kv,tr -> (rectifyTau) -> X -> calcTNB();
	calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();
	calcAlpha( flg_rectifyA );	// kv,k2d -> (rectifyAlpha) -> alpha, sina, cosa, tana, da
	calcRuling( flg_rectifyR );

	//setP_k(-1);
	//setP_t(-1);
	//setP_k2(-1);
	setP_a(-1);

	return 0;
}

