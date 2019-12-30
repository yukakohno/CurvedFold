#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"

// a) 3D�Ȑ��Ɛ܂�p�x����A�܂�������߂�
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
	setP_Bl(-1);
	setP_Br(-1);

	return 0;
}

// b) �܂���Ɛ܂�p�x����A3D�Ȑ������߂�
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
	setP_Bl(-1);
	setP_Br(-1);

	return 0;
}


// c) �܂����3D�Ȑ�����A�܂�p�x�����߂�
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
	setP_Bl(-1);
	setP_Br(-1);

	return 0;
}

int crease::calcR_TA( int flg_interpolate, rectify_params *rp, int mini, int maxi, double malpha )
{
	int ret = 0;
	if( flg_interpolate ){
		initX();	Xcnt = XcntOrg;
		interpolate_spline( Px2d, Pcnt, k2d, Xcnt );
		interpolate_spline( Pa, Pcnt, alpha, Xcnt );
		interpolate_spline_RulAngle( Pbl, Pcnt,	betal, Xcnt, cotbl, cosbl, sinbl );
		interpolate_spline_RulAngle( Pbr, Pcnt,	betar, Xcnt, cotbr, cosbr, sinbr );
	}
	if( flg_interpolate ){
		calcXTN2d( m2 );		// k2d -> X2d -> calcTN2d();
	}
	if( -M_PI < malpha && malpha < M_PI ){
		if( mini < 0 ){ mini = 0; }
		if( maxi < 0 ){ maxi = Xcnt; }
		int m = (maxi+mini)/2;
		alpha[m]=malpha;
	}
	ret = calcRul2TA( 1, mini, maxi );	// Br, Bl -> tr,alpha
	if( rp->flg_rectifyA ){ rectifyAlphaBezier( &rp->rectA ); }
	calcAK2D_K();				// alpha, k2d -> kv
	if( rp->flg_rectifyT ){ rectifyTauBezier2( &rp->rectT ); }
	calcXTNB( m3 );				// kv,tr -> X -> calcTNB();
	calcDA();
	calcRuling( rp->flg_rectifyR, rp->rectifyR_kvthres );

	setP_k(-1);
	setP_t(-1);
	//setP_k2(-1);
	setP_a(-1);
	//setP_Bl(-1);
	//setP_Br(-1);

	return ret;
}

int crease::calcR_TA2( rectify_params *rp,
					  int mini, int maxi )
{
#if 0
	memset(fileioflg,0,sizeof(int)*30);
	fileioflg[D2] = 1;
	fileioflg[K2] = 1;
	fileioflg[D3] = 1;
	fileioflg[K3] = 1;
	fileioflg[T3] = 1;
	fileioflg[ALPHA] = 1;
	fileioflg[DA] = 1;
	fileioflg[BETA] = 1;
	dumpAll( "output/0.csv" );
#endif
	//calcXTN2d( m2 );			// k2d -> X2d -> calcTN2d();	// �܂���͕ς��Ȃ�
	calcRul2TA( 1, mini, maxi );	// Br, Bl -> tr,alpha
	//dumpAll( "output/1.csv" );
	if( rp->flg_rectifyA ){ rectifyAlphaBezier( &rp->rectA ); }
	calcAK2D_K();				// alpha, k2d -> kv
	//dumpAll( "output/2.csv" );
	if( rp->flg_rectifyT ){ rectifyTauBezier2( &rp->rectT ); }
	calcXTNB( m3 );				// kv,tr -> X -> calcTNB();
	//dumpAll( "output/3.csv" );
	calcDA();
	//dumpAll( "output/4.csv" );
	calcRuling( rp->flg_rectifyR, rp->rectifyR_kvthres );
	//dumpAll( "output/5.csv" );

	setP_k(-1);
	setP_t(-1);
	//setP_k2(-1);
	setP_a(-1);
	//setP_Bl(-1);
	//setP_Br(-1);

	return 0;
}
