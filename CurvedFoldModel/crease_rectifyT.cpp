#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "crease.h"
#include "Bezier.h"

#if 1

int crease::rectifyTau( int Xcnt, int *s1, int *e1, int scnt, double *tr, double *kv )
{
	int i, ret=0;

	if( scnt==0 ){
		goto end;
	}

	// tau ï‚ê≥
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

	// tau ï‚ê≥
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

int crease::rectifyTau( double rectifyT_kvthres, int rectifyT_mgn )
{
	int s1[10],e1[10], scnt=0, ret=0;
	//scnt = gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, MIN_CURV, 10, s1, e1, 10 );
	scnt = gets1e1( Xcnt, CEMGN, Xcnt-CEMGN-1, k2d, rectifyT_kvthres, rectifyT_mgn, s1, e1, 10 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }
	ret = rectifyTau( Xcnt, s1, e1, scnt, tr, kv );
end:
	return ret;
}

int crease::rectifyTau2( int src_s1e1, double rectifyT_kvthres, int rectifyT_mgn )
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

	// É—=0 Ç…Ç∑ÇÈîÕàÕ -> s1,e1
	//scnt = gets1e1( flg_src_s1e1, MIN_CURV, 1, s1, e1, 10 );
	scnt = gets1e1( src_s1e1, rectifyT_kvthres, rectifyT_mgn, s1, e1, 10 );
	if( scnt==0 ) goto end;
	if( scnt<0 ){ ret = -1; goto end; }

	// É—Çï‚ê≥Ç∑ÇÈîÕàÕ -> s2,e2
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

	// tau ï‚ê≥
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

int crease::rectifyTauBezier( rectify_param *rp )
{
	int i, /*s1[10],e1[10], scnt=0,*/ ret=0;
	Bezier b;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifytau.csv", "w" );

	if( rp->src==0 ){
		//scnt = gets1e1( s1, e1, 10, 10, flg_src_s1e1 ); // arraysize, mgn
		rp->scnt = gets1e1( rp->flg_src_s1e1, rp->kvthres, rp->s1mgn, rp->s1, rp->e1, RECT_ASIZE );
	}
	if( rp->scnt==0 ) goto end;
	if( rp->scnt<0 ){ ret = -1; goto end; }

	// tau ï‚ê≥
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

	for( int k=0; k<rp->scnt; k++ ){
		int s1=rp->s1[k], e1=rp->e1[k], s2=rp->s2[k], e2=rp->e2[k];
		double dx0,dy0,dx1,dy1;
		if( s1>0 ){
			dx0 = 2;
			dy0 = t_k[ s1+1 ] - t_k[ s1-1 ];
		} else {
			dx0 = 1;
			dy0 = t_k[ s1+1 ] - t_k[ s1 ];
		}
		if( e1<Xcnt-1 ){
			dx1 = -2;
			dy1 = t_k[ e1-1 ] - t_k[ e1+1 ];
		} else {
			dx1 = -1;
			dy1 = t_k[ e1-1 ] - t_k[ e1 ];
		}

		b.set( s1, t_k[s1], dx0, dy0, rp->svlen0[k], e1, t_k[e1], dx1, dy1, rp->evlen0[k], NULL/*"b.csv"*/ );

		for( i=s1+1; i<e1; i++ ){
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

int crease::rectifyTauBezier2( rectify_param *rp )
{
	int /*i, s1[10],e1[10], s2[10],e2[10], scnt=0,*/ ret=0;
	double tr0[MAX_SPCNT], t_k[MAX_SPCNT], cumtr0_pos=0.0, cumtr0_neg=0.0, cumtr1_pos=0.0, cumtr1_neg=0.0;
	Bezier b0, b1;
	FILE *fp = NULL;
	//fp = fopen( "output/rectifytau.csv", "w" );

	if( tr[0]==0.0 && tr[1]==0.0 && tr[2]==0.0
		&& tr[Xcnt-1]==0.0 && tr[Xcnt-2]==0.0 && tr[Xcnt-3]==0.0 ){
			tr[0] = tr[1] = tr[2] = tr[3];
			tr[Xcnt-1] = tr[Xcnt-2] = tr[Xcnt-3] = tr[Xcnt-4];
	}

	memcpy( tr0, tr, sizeof(double)*Xcnt );
	for( int i=0; i<Xcnt; i++ ){
		if( kv[i] != 0.0 ){
			t_k[i] = tr[i]/kv[i];
		}
	}
	if( fp ){
		fprintf( fp, "original tr" );
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original kv" );
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "original t/k" );
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

	if( rp->src==0 ){
		//scnt = gets1e1( s1, e1, 10, 0, flg_src_s1e1 ); // arraysize, mgn
		rp->scnt = gets1e1( rp->flg_src_s1e1, rp->kvthres, rp->s1mgn, rp->s1, rp->e1, RECT_ASIZE );
	}
	if( rp->scnt==0 ) goto end;
	if( rp->scnt<0 ){ ret = -1; goto end; }

	// get s2[], e2[]
	if( rp->src==0 || rp->src==1 ){
		ret = gets2e2( Xcnt, rp->s1, rp->e1, rp->scnt, k2d, rp->s2, rp->e2 );
	}

	for( int k=0; k<rp->scnt; k++ ){
		int s1=rp->s1[k], e1=rp->e1[k], s2=rp->s2[k], e2=rp->e2[k];
		// [s2,s1], [e1,e2] Çï‚ä‘
		double dx,dy, len0, len1, ave=0.0;

		if( s2>0 ){
			dx = 2;
			dy = t_k[ s2+1 ] - t_k[ s2-1 ];
		} else {
			dx = 1;
			dy = t_k[ s2+1 ] - t_k[ s2 ];
		}
		len0 = len1 = (s1-s2)*0.5;
		if( len0 > rp->svlen0[k] ){ len0 = rp->svlen0[k]; }
		if( len1 > rp->svlen1[k] ){ len1 = rp->svlen1[k]; }
		b0.set( s2, t_k[s2], dx, dy, len0, s1, ave, -1.0, 0.0, len1, NULL );

		if( e2<Xcnt-1 ){
			dx = -2;
			dy = t_k[ e2-1 ] - t_k[ e2+1 ];
		} else {
			dx = -1;
			dy = t_k[ e2-1 ] - t_k[ e2 ];
		}
		len0 = len1 = (e2-e1)*0.5;
		if( len0 > rp->evlen0[k] ){ len0 = rp->evlen0[k]; }
		if( len1 > rp->evlen1[k] ){ len1 = rp->evlen1[k]; }
		b1.set( e1, ave, 1.0, 0.0, len0, e2, t_k[e2], dx, dy, len1, NULL );

		for( int i=s2; i<s1; i++ ){
			t_k[i] = b0.gety(i);
			tr0[i] = t_k[i] * kv[i];
		}
		for( int i=s1; i<=e1; i++ ){
			t_k[i] = tr0[i] = 0.0;
		}
		for( int i=e1+1; i<=e2; i++ ){
			t_k[i] = b1.gety(i);
			tr0[i] = t_k[i] * kv[i];
		}
	}

#if 1 // tr ëçòaÇ™ìØÇ∂Ç…Ç»ÇÈÇÊÇ§Ç…
	for( int i=0; i<Xcnt; i++ ){
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

	for( int i=0; i<Xcnt; i++ ){
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
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", tr[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify kv" );
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", kv[i] );
		}
		fprintf( fp, "\n" );
		fprintf( fp, "rectify t/k" );
		for( int i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", t_k[i] );
		}
		fprintf( fp, "\n" );
	}

end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}
