#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "crease.h"
#include "Bezier.h"

int crease::gets1e1( int Xcnt, int Xs, int Xe, double *val, double thres_val, int mgn,
					int *_s1, int *_e1, int semax )
{
	int i, j, s0=-1,e0=-1, scnt=0, scnt0=0;
	int flg[MAX_SPCNT]; memset( flg, 0, sizeof(int)*MAX_SPCNT );

	// curv=0 に近い範囲を検出 -> s0,e0
	for( i=Xs; i<=Xe; i++ ){
		double aval = fabs(val[i]);
		if( s0>-1 ){
			flg[i] = 1;
		}
		if( s0<0 && aval<thres_val ){
			_s1[scnt] = i-1;
			flg[i-1] = flg[i] = 1;
			s0 = 1;
		} else if( s0>0 && aval>thres_val ){
			_e1[scnt] = i;
			scnt++;
			s0 = -1;
			if( scnt >= semax){
				printf("scnt overflow");
				scnt = -1; goto end;
			}
		}
	}
	if( s0>0 ){
		_e1[scnt] = i;
		scnt++;
		s0 = -1;
	}

	// + → -、- → + を追加で検出（勾配が大きい場合） -> s0,e0
	for( i=Xs; i<Xe; i++ ){
		if( flg[i] || flg[i+1] ){
			continue;
		}
		if( val[i] < 0 && val[i+1] > 0 || val[i] > 0 && val[i+1] < 0 ){

			if( scnt==0 ){
				_s1[scnt] = i;
				_e1[scnt] = i+1;
				scnt++;
				continue;
			}

			// insert
			scnt++;
			if( scnt > semax){
				printf("scnt overflow");
				scnt = -1; goto end;
			}
			int ins=scnt-1;
			for( j=0; j<scnt-1; j++ ){
				if( i < _s1[j] ){
					ins = j;
					break;
				}
			}
			for( j=scnt-2; j>=ins; j-- ){
				_s1[j+1] = _s1[j];
				_e1[j+1] = _e1[j];
			}
			_s1[ins] = i;
			_e1[ins] = i+1;
		}
	}
	if( scnt==0 ){
		goto end;
	}

	//	for( i=0; i<scnt; i++ ){
	//		printf( "%d: s0=%d, e0=%d\n", i, _s1[i], _e1[i] );
	//	}

	// curv=0 前後の補正範囲を設定 -> s1,e1
	for( i=0; i<scnt; i++ ){
		_s1[i] = _s1[i]-mgn;
		if( _s1[i]<0 ){ _s1[i]=0; }
		_e1[i] = _e1[i]+mgn;
		if( _e1[i]>=Xcnt ){ _e1[i]=Xcnt-1; }
	}

	// 隣とくっついた場合は統合
	for( i=1; i<scnt; i++ ){
		if( _s1[i] < _e1[i-1] ){
			_s1[i] = _s1[i-1];
			_s1[i-1] = _e1[i-1] = -1;
		}
	}
	for( i=0; i<scnt; i++ ){
		if( _s1[i] > -1 && _e1[i] > -1 ){
			_s1[scnt0] = _s1[i];
			_e1[scnt0] = _e1[i];
			scnt0++;
		}
	}
	scnt = scnt0;

	//	for( i=0; i<scnt; i++ ){
	//		printf( "%d: s1=%d, e1=%d\n", i, _s1[i], _e1[i] );
	//	}
end:
	return scnt;
}

// abs(k2d[i])>k0thres の範囲、符号が変わる点、をmgnだけ広げた範囲を求める
int crease::gets1e1( int *_s1, int *_e1, int semax, int mgn, int flg_src_s1e1 )
{
	if( flg_src_s1e1==0 ){
		return gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, MIN_CURV, mgn, _s1, _e1, semax );
	} else {
		double ak2d[MAX_SPCNT] = {0.0};
		for( int i=CEMGN; i<Xcnt-CEMGN; i++ ){
			ak2d[i] = kv[i]*cos(alpha[i]);
		}
		return gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, ak2d, MIN_CURV, mgn, _s1, _e1, semax );
	}
}
#if 1

int crease::rectifyTau( int Xcnt, int *s1, int *e1, int scnt, double *tr, double *kv )
{
	int i, ret=0;

	if( scnt==0 ){
		goto end;
	}

	// tau 補正
	double t_k[MAX_SPCNT], di, diunit;
	for( i=0; i<Xcnt; i++ ){
		if( kv[i] != 0.0 ){
			t_k[i] = tr[i]/kv[i];
		}
	}
#if 1
	for( int k=0; k<scnt; k++ ){
		diunit=20.0/(e1[k]-s1[k]);
		for( i=s1[k], di=-10.0; i<=e1[k]; i++, di+=diunit ){
			t_k[i] = 1.0/(1.0+exp(-di*0.5)) * (t_k[e1[k]]-t_k[s1[k]]) + t_k[s1[k]];
			tr[i] = t_k[i] * kv[i];
		}
	}
#else
	diunit=20.0/(s0-s1);
	for( i=s1, di=10.0; i<=s0; i++, di-=diunit ){
		t_k[i] = 1.0/(1.0+exp(-di*0.5)) * t_k[s1];
		tr[i] = t_k[i] * kv[i];
	}
	diunit=20.0/(e1-e0);
	for( i=e1, di=10.0; i>=e0; i--, di-=diunit ){
		t_k[i] = 1.0/(1.0+exp(-di*0.5)) * t_k[e1];
		tr[i] = t_k[i] * kv[i];
	}
	for( i=s0; i<=e0; i++ ){
		tr[i] = 0.0;
	}
#endif
end:
	return ret;
}

#else
int crease::rectifyTau( int Xcnt, double *kv, double *alpha, double *tr )
{
	int i, s1[10],e1[10], scnt=0, ret=0;

	scnt = gets1e1(Xcnt, kv, alpha, 10, s1, e1, 10 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	// tau 補正
	double t_k[MAX_SPCNT], di, diunit;
	for( i=0; i<Xcnt; i++ ){
		if( kv[i] != 0.0 ){
			t_k[i] = tr[i]/kv[i];
		}
	}
#if 1
	for( int k=0; k<scnt; k++ ){
		diunit=20.0/(e1[k]-s1[k]);
		for( i=s1[k], di=-10.0; i<=e1[k]; i++, di+=diunit ){
			t_k[i] = 1.0/(1.0+exp(-di*0.5)) * (t_k[e1[k]]-t_k[s1[k]]) + t_k[s1[k]];
			tr[i] = t_k[i] * kv[i];
		}
	}
#else
	diunit=20.0/(s0-s1);
	for( i=s1, di=10.0; i<=s0; i++, di-=diunit ){
		t_k[i] = 1.0/(1.0+exp(-di*0.5)) * t_k[s1];
		tr[i] = t_k[i] * kv[i];
	}
	diunit=20.0/(e1-e0);
	for( i=e1, di=10.0; i>=e0; i--, di-=diunit ){
		t_k[i] = 1.0/(1.0+exp(-di*0.5)) * t_k[e1];
		tr[i] = t_k[i] * kv[i];
	}
	for( i=s0; i<=e0; i++ ){
		tr[i] = 0.0;
	}
#endif
end:
	return ret;
}
#endif

int crease::rectifyTau()
{
	int s1[10],e1[10], scnt=0, ret=0;
	scnt = gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, MIN_CURV, 10, s1, e1, 10 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }
	ret = rectifyTau( Xcnt, s1, e1, scnt, tr, kv );
end:
	return ret;
}

int crease::rectifyTau2()
{
	int i, s1[10],e1[10], s2[10],e2[10], scnt=0, ret=0;
	double tr0[MAX_SPCNT], cumtr0_pos=0.0, cumtr0_neg=0.0, cumtr1_pos=0.0, cumtr1_neg=0.0;
	Bezier b0, b1;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifytau.csv", "w" );

	for( i=0; i<Xcnt; i++ ){
		if( tr[i]>0.0 ){
			cumtr0_pos += tr[i];
		} else {
			cumtr0_neg += tr[i];
		}
	}

	memcpy( tr0, tr, sizeof(double)*Xcnt );

	// τ=0 にする範囲 -> s1,e1
	scnt = gets1e1( s1, e1, 10, 1, flg_src_s1e1 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	// τを補正する範囲 -> s2,e2
	int m[10];
	for( i=0; i<scnt-1; i++ ){
		m[i] = (int)( ( e1[i] + s1[i+1] ) /2.0 );
	}
	for( i=0; i<scnt; i++ ){
		if( i==0 ){
			s2[i] = CEMGN;
		} else {
			s2[i] = m[i-1]+1;
		}
		if( i==scnt-1 ){
			e2[i] = Xcnt-CEMGN;
		} else {
			e2[i] = m[i]-1;
		}
	}

	//for( i=0; i<scnt; i++ ){
	//	printf( "%d: s2=%d, e2=%d\n", i, s2[i], e2[i] );
	//}

	// tau 補正
	for( int k=0; k<scnt; k++ ){
		double ave, dx,dy;
		ave = 0.0;

		if( s2[k]>CEMGN ){
			dx = 2;
			dy = tr[ s2[k]+1 ] - tr[ s2[k]-1 ];
		} else {
			dx = 1;
			dy = tr[ s2[k]+1 ] - tr[ s2[k] ];
		}
		b0.set( s2[k], tr[s2[k]], dx, dy, (s1[k]-s2[k])*0.5,
			s1[k], ave, -1.0, 0.0, (s1[k]-s2[k])*0.5, NULL /*"output/b0.csv"*/ );

		if( e2[k]<Xcnt-CEMGN ){
			dx = -2;
			dy = tr[ e2[k]-1 ] - tr[ e2[k]+1 ];
		} else {
			dx = -1;
			dy = tr[ e2[k]-1 ] - tr[ e2[k] ];
		}
		b1.set( e1[k], ave, 1.0, 0.0, (e2[k]-e1[k])*0.5,
			e2[k], tr[e2[k]], dx, dy, (e2[k]-e1[k])*0.5, NULL /*"ouput/b1.csv"*/ );

		for( i=s2[k]+1; i<s1[k]; i++ ){
			tr0[i] = b0.gety( (double)i );
		}
		for( i=s1[k]; i<=e1[k]; i++ ){
			tr0[i] = ave;
		}
		for( i=e1[k]+1; i<e2[k]; i++ ){
			tr0[i] = b1.gety( (double)i );
		}

		if( fp ){
			fprintf( fp, "iter %d", k );
			for( i=0; i<Xcnt; i++){
				fprintf( fp, ",%f", tr0[i] );
			}
			fprintf( fp, "\n" );
		}

	}
	memcpy( tr, tr0, sizeof(double)*Xcnt );

	for( i=0; i<Xcnt; i++ ){
		if( tr[i]==0.0 ){
		} else if( tr[i]>0.0 ){
			cumtr1_pos += tr[i];
		} else {
			cumtr1_neg += tr[i];
		}
	}
	if( cumtr1_pos!=0.0 ){ cumtr0_pos /= cumtr1_pos; }
	if( cumtr1_neg!=0.0 ){ cumtr0_neg /= cumtr1_neg; }

	for( i=0; i<Xcnt; i++ ){
		if( tr[i]==0.0 ){
		} else if( tr[i]>0.0 ){
			tr[i] *= cumtr0_pos;
		} else {
			tr[i] *= cumtr0_neg;
		}
	}

end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}

int crease::rectifyTauBezier()
{
	int i, s1[10],e1[10], scnt=0, ret=0;
	Bezier b;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifytau.csv", "w" );

	scnt = gets1e1( s1, e1, 10, 10, flg_src_s1e1 ); // arraysize, mgn
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	// tau 補正
	double t_k[MAX_SPCNT];
	for( i=0; i<Xcnt; i++ ){
		if( kv[i] != 0.0 ){
			t_k[i] = tr[i]/kv[i];
		}
	}

	if( fp ){
		fprintf( fp, "original tr" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original kv" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original t/k" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

	for( int k=0; k<scnt; k++ ){

		double dx0,dy0,dx1,dy1;
		if( s1[k]>0 ){
			dx0 = 2;
			dy0 = t_k[ s1[k]+1 ] - t_k[ s1[k]-1 ];
		} else {
			dx0 = 1;
			dy0 = t_k[ s1[k]+1 ] - t_k[ s1[k] ];
		}
		if( e1[k]<Xcnt-1 ){
			dx1 = -2;
			dy1 = t_k[ e1[k]-1 ] - t_k[ e1[k]+1 ];
		} else {
			dx1 = -1;
			dy1 = t_k[ e1[k]-1 ] - t_k[ e1[k] ];
		}

		b.set( s1[k], t_k[s1[k]], dx0, dy0, 10,
			e1[k], t_k[e1[k]], dx1, dy1, 10, NULL/*"b.csv"*/ );

		for( i=s1[k]+1; i<e1[k]; i++ ){
			t_k[i] = b.gety( (double)i );
			tr[i] = t_k[i] * kv[i];
		}
	}

	if( fp ){
		fprintf( fp, "rectify tr" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify kv" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify t/k" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}

int crease::rectifyTauBezier2()
{
	int i, s1[10],e1[10], s2[10],e2[10], scnt=0, ret=0;
	double tr0[MAX_SPCNT], t_k[MAX_SPCNT], cumtr0_pos=0.0, cumtr0_neg=0.0, cumtr1_pos=0.0, cumtr1_neg=0.0;
	Bezier b0, b1;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifytau.csv", "w" );

	memcpy( tr0, tr, sizeof(double)*Xcnt );
	for( i=0; i<Xcnt; i++ ){
		if( kv[i] != 0.0 ){
			t_k[i] = tr[i]/kv[i];
		}
	}
	if( fp ){
		fprintf( fp, "original tr" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original kv" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original t/k" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

	scnt = gets1e1( s1, e1, 10, 0, flg_src_s1e1 ); // arraysize, mgn
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	int m[10];
	for( i=0; i<scnt-1; i++ ){
		m[i] = (int)( ( e1[i] + s1[i+1] ) /2.0 );
	}
	for( i=0; i<scnt; i++ ){
		// まず最大まで広げる
		if( i==0 ){
			s2[i] = CEMGN;
		} else {
			s2[i] = m[i-1]+1;
		}
		if( i==scnt-1 ){
			e2[i] = Xcnt-CEMGN-1;
		} else {
			e2[i] = m[i]-1;
		}
		// 最も近い極大／極小まで狭める
		double sgn = k2d[ s1[i]-1 ] - k2d[ s1[i] ];
		for( int j=s1[i]; j>s2[i]; j-- ){
			if( sgn>=0 && k2d[j-1]<k2d[j] || sgn<0 && k2d[j-1]>k2d[j] ){
				s2[i] = j;
				break;
			}
		}
		sgn = k2d[ s1[i]+1 ] - k2d[ s1[i] ];
		for( int j=e1[i]; j<e2[i]; j++ ){
			if( sgn>=0 && k2d[j+1]<k2d[j] || sgn<0 && k2d[j+1]>k2d[j] ){
				e2[i] = j;
				break;
			}
		}
	}

	for( int k=0; k<scnt; k++ ){
		// [s2,s1], [e1,e2] を補間
		double dx,dy, len0, len1, maxlen=3.0, ave=0.0; // ★TODO: ave=0.0で問題ないか再考

		if( s2[k]>0 ){
			dx = 2;
			dy = t_k[ s2[k]+1 ] - t_k[ s2[k]-1 ];
		} else {
			dx = 1;
			dy = t_k[ s2[k]+1 ] - t_k[ s2[k] ];
		}
		len0 = len1 = (s1[k]-s2[k])*0.5;
		len0 = len0 < maxlen ? len0 : maxlen;
		len1 = len1 < maxlen ? len1 : maxlen;
		b0.set( s2[k], t_k[s2[k]], dx, dy, len0,
			s1[k], ave, -1.0, 0.0, len1, NULL );

		if( e2[k]<Xcnt-1 ){
			dx = -2;
			dy = t_k[ e2[k]-1 ] - t_k[ e2[k]+1 ];
		} else {
			dx = -1;
			dy = t_k[ e2[k]-1 ] - t_k[ e2[k] ];
		}
		len0 = len1 = (e2[k]-e1[k])*0.5;
		len0 = len0 < maxlen ? len0 : maxlen;
		len1 = len1 < maxlen ? len1 : maxlen;
		b1.set( e1[k], ave, 1.0, 0.0, len0,
			e2[k], t_k[e2[k]], dx, dy, len1, NULL );

		for( i=s2[k]; i<s1[k]; i++ ){
			t_k[i] = b0.gety(i);
			tr0[i] = t_k[i] * kv[i];
		}
		for( i=s1[k]; i<=e1[k]; i++ ){
			t_k[i] = tr0[i] = 0.0;
		}
		for( i=e1[k]+1; i<=e2[k]; i++ ){
			t_k[i] = b1.gety(i);
			tr0[i] = t_k[i] * kv[i];
		}
	}

#if 1 // tr 総和が同じになるように
	for( i=0; i<Xcnt; i++ ){
		if( tr[i]>0.0 ){
			cumtr0_pos += tr[i];
		} else if( tr[i]<0.0 ){
			cumtr0_neg += tr[i];
		}
		if( tr0[i]>0.0 ){
			cumtr1_pos += tr0[i];
		} else if( tr0[i]<0.0 ){
			cumtr1_neg += tr0[i];
		}
	}
	if( cumtr1_pos!=0.0 ){ cumtr0_pos /= cumtr1_pos; }
	if( cumtr1_neg!=0.0 ){ cumtr0_neg /= cumtr1_neg; }

	for( i=0; i<Xcnt; i++ ){
		if( tr0[i]==0.0 ){
		} else if( tr0[i]>0.0 ){
			tr0[i] *= cumtr0_pos;
		} else {
			tr0[i] *= cumtr0_neg;
		}
	}
#endif
	memcpy( tr, tr0, sizeof(double)*Xcnt );

	if( fp ){
		fprintf( fp, "rectify tr" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify kv" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify t/k" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}

int crease::rectifyAlpha()
{
	int i, s1[10],e1[10], s2[10],e2[10], scnt=0, ret=0;
	double alpha0[MAX_SPCNT];
	memcpy( alpha0, alpha, sizeof(double)*Xcnt );

	scnt = gets1e1( s1, e1, 10, 10, flg_src_s1e1 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret=-1; goto end; }
	for( i=0; i<scnt; i++){
		s2[i] = s1[i]-10;
		if( s2[i]<0 ){ s2[i]=0; }
		e2[i] = e1[i]+10;
		if( e2[i]>=Xcnt ){ e2[i]=Xcnt-1; }
	}
	// alpha 補正（両端の間で一定の値にして、外側を補間）
	for( int k=0; k<scnt; k++ ){
		double ave = (alpha[s1[k]]+alpha[e1[k]])*0.5;
		double di, diunit;
		for( i=s1[k]; i<=e1[k]; i++ ){
			alpha0[i] = ave;
		}
		diunit=20.0/(s1[k]-s2[k]);
		for( i=s2[k], di=-10.0; i<s1[k]; i++, di+=diunit ){
			double tmp = 1.0/(1.0+exp(-di*0.5)) * (ave - alpha[s2[k]]) + alpha[s2[k]];
			//alpha0[i] = alpha0[i] < tmp ? alpha0[i] : tmp;
			alpha0[i] = tmp;
		}
		diunit=20.0/(e2[k]-e1[k]);
		for( i=e2[k], di=-10.0; i>e1[k]; i--, di+=diunit ){
			double tmp = 1.0/(1.0+exp(-di*0.5)) * (ave - alpha[e2[k]]) + alpha[e2[k]];
			//alpha0[i] = alpha0[i] < tmp ? alpha0[i] : tmp;
			alpha0[i] = tmp;
		}
	}
	memcpy( alpha, alpha0, sizeof(double)*Xcnt );
end:
	return ret;
}

int crease::rectifyAlphaBezier( int Xcnt, int *s1, int *e1, int scnt, double *alpha )
{
	int i, s2[10],e2[10], mgn2=5, ret=0;
	double alpha0[MAX_SPCNT], maxlen=3; // XCNT*0.05; // Xcnt=100 -> maxlen=5
	Bezier b0, b1;

	FILE *fp = NULL;
	//fp = fopen( "output/rectifyalpha.csv", "w" );
	if( fp ){
		fprintf( fp, "original" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha[i] );
		}
		fprintf( fp, "\n" );
	}

	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret=-1; goto end; }

	memcpy( alpha0, alpha, sizeof(double)*Xcnt );

#if 0	// 指定量広げる
	for( i=0; i<scnt; i++){
		s2[i] = s1[i]-mgn2;
		if( s2[i]<0 ){ s2[i]=0; }
		e2[i] = e1[i]+mgn2;
		if( e2[i]>=Xcnt ){ e2[i]=Xcnt-1; }
	}
#else
	int m[10];
	for( i=0; i<scnt-1; i++ ){
		m[i] = (int)( ( e1[i] + s1[i+1] ) /2.0 );
	}
	for( i=0; i<scnt; i++ ){
		// まず最大まで広げる
		if( i==0 ){
			s2[i] = CEMGN;
		} else {
			s2[i] = m[i-1]+1;
		}
		if( i==scnt-1 ){
			e2[i] = Xcnt-CEMGN-1;
		} else {
			e2[i] = m[i]-1;
		}
		// 最も近い極大／極小まで狭める
		double sgn = alpha[ s1[i]-1 ] - alpha[ s1[i] ];
		for( int j=s1[i]; j>s2[i]; j-- ){
			if( sgn>=0 && alpha[j-1]<alpha[j] || sgn<0 && alpha[j-1]>alpha[j] ){
				s2[i] = j;
				break;
			}
		}
		sgn = alpha[ s1[i]+1 ] - alpha[ s1[i] ];
		for( int j=e1[i]; j<e2[i]; j++ ){
			if( sgn>=0 && alpha[j+1]<alpha[j] || sgn<0 && alpha[j+1]>alpha[j] ){
				e2[i] = j;
				break;
			}
		}
	}
#endif
	//for( i=0; i<scnt; i++ ){
	//	printf( "%d: s2=%d, e2=%d\n", i, s2[i], e2[i] );
	//}

#if 1	// alpha 補正(linear)
	for( int k=0; k<scnt; k++ ){
		double ave, da_s,da_e;
		ave = (alpha[s1[k]]+alpha[e1[k]])*0.5;
		da_s = (ave - alpha[s1[k]]) / (double)(s1[k]-s2[k]+1);
		da_e = (ave - alpha[e1[k]]) / (double)(e2[k]-e1[k]+1);
		for( i=s2[k]+1; i<=s1[k]+1; i++ ){	// i<=s1[k]+1: da への影響を考えて内側にはみ出させる
			alpha0[i] = alpha[i] + da_s*(double)(i-s2[k]);
		}
		for( i=s1[k]+2; i<e1[k]-1; i++ ){
			alpha0[i] = ave;
		}
		for( i=e1[k]-1; i<e2[k]; i++ ){	// i=e1[k]-1: da への影響を考えて内側にはみ出させる
			alpha0[i] = alpha[i] + da_e*(double)(e2[k]-i);
		}

		if( fp ){
			fprintf( fp, "iter %d", k );
			for( i=0; i<Xcnt; i++){
				fprintf( fp, ",%f", alpha0[i] );
			}
			fprintf( fp, "\n" );
		}
	}
#else	// alpha 補正(Bezier)
	for( int k=0; k<scnt; k++ ){
		double ave, dx,dy, len0, len1;
		ave = (alpha[s1[k]]+alpha[e1[k]])*0.5;

		if( s2[k]>0 ){
			dx = 2;
			dy = alpha[ s2[k]+1 ] - alpha[ s2[k]-1 ];
		} else {
			dx = 1;
			dy = alpha[ s2[k]+1 ] - alpha[ s2[k] ];
		}
		len0 = len1 = (s1[k]-s2[k])*0.5;
		//len0 = len0 < maxlen ? len0 : maxlen;
		//len1 = len1 < maxlen ? len1 : maxlen;
		b0.set( s2[k], alpha[s2[k]], dx, dy, len0,
			s1[k], ave, -1.0, 0.0, len1, NULL );

		if( e2[k]<Xcnt-1 ){
			dx = -2;
			dy = alpha[ e2[k]-1 ] - alpha[ e2[k]+1 ];
		} else {
			dx = -1;
			dy = alpha[ e2[k]-1 ] - alpha[ e2[k] ];
		}
		len0 = len1 = (e2[k]-e1[k])*0.5;
		//len0 = len0 < maxlen ? len0 : maxlen;
		//len1 = len1 < maxlen ? len1 : maxlen;
		b1.set( e1[k], ave, 1.0, 0.0, len0,
			e2[k], alpha[e2[k]], dx, dy, len1, NULL );

		for( i=s2[k]+1; i<s1[k]; i++ ){
			alpha0[i] = b0.gety( (double)i );
		}
		for( i=s1[k]; i<=e1[k]; i++ ){
			alpha0[i] = ave;
		}
		for( i=e1[k]+1; i<e2[k]; i++ ){
			alpha0[i] = b1.gety( (double)i );
		}

		if( fp ){
			fprintf( fp, "iter %d", k );
			for( i=0; i<Xcnt; i++){
				fprintf( fp, ",%f", alpha0[i] );
			}
			fprintf( fp, "\n" );
		}
	}
#endif
	memcpy( alpha, alpha0, sizeof(double)*Xcnt );

	if( fp ){
		fprintf( fp, "final" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha[i] );
		}
		fprintf( fp, "\n" );
	}
end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}


int crease::rectifyAlphaBezier( int flg_src_s1e1 )
{
#if 0
	int s1[10],e1[10], scnt=0, ret=0;
	if( src_s1e1 ){
		double k2d0[MAX_SPCNT]; memset(k2d0, 0, sizeof(double)*MAX_SPCNT);
		for( int i=CEMGN; i<Xcnt-CEMGN; i++ ){
			k2d0[i] = kv[i]*cos(alpha[i]);
		}
		scnt = gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d0, MIN_CURV, 10, s1, e1, 10 );
	} else {
		scnt = gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, MIN_CURV, 10, s1, e1, 10 );
	}
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	return rectifyAlphaBezier( Xcnt, s1, e1, scnt, alpha );
#else
	int i, s1[10],e1[10], s2[10],e2[10], scnt=0, mgn2=5, ret=0;
	double alpha0[MAX_SPCNT], maxlen=3; // XCNT*0.05; // Xcnt=100 -> maxlen=5
	Bezier b0, b1;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifyalpha.csv", "w" );

	if( fp ){
		fprintf( fp, "original" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha[i] );
		}
		fprintf( fp, "\n" );
	}

	memcpy( alpha0, alpha, sizeof(double)*Xcnt );

	scnt = gets1e1( s1, e1, 10, 0, flg_src_s1e1 ); // arraysize, mgn
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret=-1; goto end; }

#if 0	// 指定量広げる
	for( i=0; i<scnt; i++){
		s2[i] = s1[i]-mgn2;
		if( s2[i]<0 ){ s2[i]=0; }
		e2[i] = e1[i]+mgn2;
		if( e2[i]>=Xcnt ){ e2[i]=Xcnt-1; }
	}
#else
	int m[10];
	for( i=0; i<scnt-1; i++ ){
		m[i] = (int)( ( e1[i] + s1[i+1] ) /2.0 );
	}
	for( i=0; i<scnt; i++ ){
		// まず最大まで広げる
		if( i==0 ){
			s2[i] = CEMGN;
		} else {
			s2[i] = m[i-1]+1;
		}
		if( i==scnt-1 ){
			e2[i] = Xcnt-CEMGN-1;
		} else {
			e2[i] = m[i]-1;
		}
		// 最も近い極大／極小まで狭める
		double sgn = alpha[ s1[i]-1 ] - alpha[ s1[i] ];
		for( int j=s1[i]; j>s2[i]; j-- ){
			if( sgn>=0 && alpha[j-1]<alpha[j] || sgn<0 && alpha[j-1]>alpha[j] ){
				s2[i] = j;
				break;
			}
		}
		sgn = alpha[ s1[i]+1 ] - alpha[ s1[i] ];
		for( int j=e1[i]; j<e2[i]; j++ ){
			if( sgn>=0 && alpha[j+1]<alpha[j] || sgn<0 && alpha[j+1]>alpha[j] ){
				e2[i] = j;
				break;
			}
		}
	}
#endif
	//for( i=0; i<scnt; i++ ){
	//	printf( "%d: s2=%d, e2=%d\n", i, s2[i], e2[i] );
	//}

#if 0	// alpha 補正(linear)
	for( int k=0; k<scnt; k++ ){
		double ave, da_s,da_e, dx,dy, len0, len1;
		ave = (alpha[s1[k]]+alpha[e1[k]])*0.5;
		da_s = (ave - alpha[s1[k]]) / (double)(s1[k]-s2[k]+1);
		da_e = (ave - alpha[e1[k]]) / (double)(e2[k]-e1[k]+1);
		for( i=s2[k]+1; i<=s1[k]+1; i++ ){	// i<=s1[k]+1: da への影響を考えて内側にはみ出させる
			alpha0[i] = alpha[i] + da_s*(double)(i-s2[k]);
		}
		for( i=s1[k]+2; i<e1[k]-1; i++ ){
			alpha0[i] = ave;
		}
		for( i=e1[k]-1; i<e2[k]; i++ ){	// i=e1[k]-1: da への影響を考えて内側にはみ出させる
			alpha0[i] = alpha[i] + da_e*(double)(e2[k]-i);
		}

		if( fp ){
			fprintf( fp, "iter %d", k );
			for( i=0; i<Xcnt; i++){
				fprintf( fp, ",%f", alpha0[i] );
			}
			fprintf( fp, "\n" );
		}
	}
#else	// alpha 補正(Bezier)
	for( int k=0; k<scnt; k++ ){
		double ave, dx,dy, len0, len1;
		ave = (alpha[s1[k]]+alpha[e1[k]])*0.5;

		if( s2[k]>0 ){
			dx = 2;
			dy = alpha[ s2[k]+1 ] - alpha[ s2[k]-1 ];
		} else {
			dx = 1;
			dy = alpha[ s2[k]+1 ] - alpha[ s2[k] ];
		}
		len0 = len1 = (s1[k]-s2[k])*0.5;
		//len0 = len0 < maxlen ? len0 : maxlen;
		//len1 = len1 < maxlen ? len1 : maxlen;
		b0.set( s2[k], alpha[s2[k]], dx, dy, len0,
			s1[k], ave, -1.0, 0.0, len1, NULL );

		if( e2[k]<Xcnt-1 ){
			dx = -2;
			dy = alpha[ e2[k]-1 ] - alpha[ e2[k]+1 ];
		} else {
			dx = -1;
			dy = alpha[ e2[k]-1 ] - alpha[ e2[k] ];
		}
		len0 = len1 = (e2[k]-e1[k])*0.5;
		//len0 = len0 < maxlen ? len0 : maxlen;
		//len1 = len1 < maxlen ? len1 : maxlen;
		b1.set( e1[k], ave, 1.0, 0.0, len0,
			e2[k], alpha[e2[k]], dx, dy, len1, NULL );

		for( i=s2[k]+1; i<s1[k]; i++ ){
			alpha0[i] = b0.gety( (double)i );
		}
		for( i=s1[k]; i<=e1[k]; i++ ){
			alpha0[i] = ave;
		}
		for( i=e1[k]+1; i<e2[k]; i++ ){
			alpha0[i] = b1.gety( (double)i );
		}

		if( fp ){
			fprintf( fp, "iter %d", k );
			for( i=0; i<Xcnt; i++){
				fprintf( fp, ",%f", alpha0[i] );
			}
			fprintf( fp, "\n" );
		}
	}
#endif
	memcpy( alpha, alpha0, sizeof(double)*Xcnt );

	if( fp ){
		fprintf( fp, "final" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha[i] );
		}
		fprintf( fp, "\n" );
	}
#endif
end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}
