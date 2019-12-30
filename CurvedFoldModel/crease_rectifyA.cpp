#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "crease.h"
#include "Bezier.h"

int crease::rectifyAlpha( rectify_param *rp )
{
	int i, ret=0;
	double alpha0[MAX_SPCNT];
	memcpy( alpha0, alpha, sizeof(double)*Xcnt );

	rp->scnt = gets1e1( rp->flg_src_s1e1, rp->kvthres, rp->s1mgn, rp->s1, rp->e1, RECT_ASIZE );
	if( rp->scnt==0 ) goto end;
	if( rp->scnt<0 ){ ret=-1; goto end; }

	for( i=0; i<rp->scnt; i++){
		rp->s2[i] = rp->s1[i] - rp->s2mgn;
		if( rp->s2[i]<0 ){ rp->s2[i]=0; }
		rp->e2[i] = rp->e1[i] + rp->s2mgn;
		if( rp->e2[i]>=Xcnt ){ rp->e2[i]=Xcnt-1; }
	}
	// alpha 補正（両端の間で一定の値にして、外側を補間）
	for( int k=0; k<rp->scnt; k++ ){
		double ave = (alpha[rp->s1[k]]+alpha[rp->e1[k]])*0.5;
		double di, diunit;
		for( i=rp->s1[k]; i<=rp->e1[k]; i++ ){
			alpha0[i] = ave;
		}
		diunit = 20.0 / (rp->s1[k] - rp->s2[k]);
		for( i=rp->s2[k], di=-10.0; i<rp->s1[k]; i++, di+=diunit ){
			double tmp = 1.0/(1.0+exp(-di*0.5)) * (ave - alpha[rp->s2[k]]) + alpha[rp->s2[k]];
			//alpha0[i] = alpha0[i] < tmp ? alpha0[i] : tmp;
			alpha0[i] = tmp;
		}
		diunit = 20.0 / (rp->e2[k] - rp->e1[k]);
		for( i=rp->e2[k], di=-10.0; i>rp->e1[k]; i--, di+=diunit ){
			double tmp = 1.0/(1.0+exp(-di*0.5)) * (ave - alpha[rp->e2[k]]) + alpha[rp->e2[k]];
			//alpha0[i] = alpha0[i] < tmp ? alpha0[i] : tmp;
			alpha0[i] = tmp;
		}
	}
	memcpy( alpha, alpha0, sizeof(double)*Xcnt );
end:
	return ret;
}

int crease::rectifyAlphaBezier( rectify_param *rp )
{
	int ret=0;
	if( rp->src==0 ){
		rp->scnt = gets1e1( rp->flg_src_s1e1, rp->kvthres, rp->s1mgn, rp->s1, rp->e1, RECT_ASIZE );
	}
	if( rp->scnt==0 ) goto end;
	if( rp->scnt<0 ){ ret=-1; goto end; }

	ret = rectifyAlphaBezier( Xcnt, alpha, rp );
end:
	return ret;
}

int crease::rectifyAlphaBezier( int Xcnt, double *alpha0, rectify_param *rp )
{
	int i, ret=0;
	double alpha1[MAX_SPCNT];

	FILE *fp = NULL;
	//fp = fopen( "output/rectifyalpha.csv", "w" );
	if( fp ){
		fprintf( fp, "original" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha0[i] );
		}
		fprintf( fp, "\n" );
	}

	if( rp->scnt==0 ) goto end;
	if( rp->scnt<0 ){ ret=-1; goto end; }

	//
	// get s2[], e2[]
	//
	if( rp->src==0 || rp->src==1 ){
		if( rp->s2mgn == 0 ){
			for( i=0; i<rp->scnt; i++){
				rp->s2[i] = rp->s1[i];
				rp->e2[i] = rp->e1[i];
			}
		} else if ( rp->s2mgn > 0 ){
			// 指定量広げる
			for( i=0; i<rp->scnt; i++){
				rp->s2[i] = rp->s1[i] - rp->s2mgn;
				if( rp->s2[i]<0 ){ rp->s2[i]=0; }
				rp->e2[i] = rp->e1[i] + rp->s2mgn;
				if( rp->e2[i]>=Xcnt ){ rp->e2[i]=Xcnt-1; }
			}
			// 隣とくっついた場合？

		} else {
			ret = gets2e2( Xcnt, rp->s1, rp->e1, rp->scnt, alpha0, rp->s2, rp->e2 );
		}
	}
	//for( i=0; i<rp->scnt; i++ ){
	//	printf( "%d: s2=%d, e2=%d\n", i, rp->s2[i], rp->e2[i] );
	//}

	//
	// alpha0 補正
	//
	memcpy( alpha1, alpha0, sizeof(double)*Xcnt );

	if( rp->method==1 ) // linear
	{
		for( int k=0; k<rp->scnt; k++ ){
			int s1=rp->s1[k], e1=rp->e1[k], s2=rp->s2[k], e2=rp->e2[k];
			double ave, da_s,da_e;
			ave = (alpha0[s1]+alpha0[e1])*0.5;
			da_s = (ave - alpha0[s1]) / (double)(s1-s2+1);
			da_e = (ave - alpha0[e1]) / (double)(e2-e1+1);
			for( i=s2+1; i<=s1+1; i++ ){	// i<=s1[k]+1: da への影響を考えて内側にはみ出させる
				alpha1[i] = alpha0[i] + da_s*(double)(i-s2);
			}
			for( i=s1+2; i<e1-1; i++ ){
				alpha1[i] = ave;
			}
			for( i=e1-1; i<e2; i++ ){	// i=e1[k]-1: da への影響を考えて内側にはみ出させる
				alpha1[i] = alpha0[i] + da_e*(double)(e2-i);
			}

			if( fp ){
				fprintf( fp, "iter %d", k );
				for( i=0; i<Xcnt; i++){
					fprintf( fp, ",%f", alpha1[i] );
				}
				fprintf( fp, "\n" );
			}
		}
	}
	else if( rp->method==2 ) // Bezier
	{
		Bezier b0, b1;
		for( int k=0; k<rp->scnt; k++ ){
			int s1=rp->s1[k], e1=rp->e1[k], s2=rp->s2[k], e2=rp->e2[k];
			double ave, dx,dy, len0, len1;
			ave = ( alpha0[s1] + alpha0[e1] )*0.5;

			if( s2>0 ){
				dx = 2;
				dy = alpha0[ s2+1 ] - alpha0[ s2-1 ];
			} else {
				dx = 1;
				dy = alpha0[ s2+1 ] - alpha0[ s2 ];
			}
			len0 = len1 = (s1-s2)*0.5;
			if( len0 > rp->svlen0[k] ){ len0 = rp->svlen0[k]; }
			if( len1 > rp->svlen1[k] ){ len1 = rp->svlen1[k]; }
			b0.set( s2, alpha0[s2], dx, dy, len0, s1, ave, -1.0, 0.0, len1, NULL );

			if( e2<Xcnt-1 ){
				dx = -2;
				dy = alpha0[ e2-1 ] - alpha0[ e2+1 ];
			} else {
				dx = -1;
				dy = alpha0[ e2-1 ] - alpha0[ e2 ];
			}
			len0 = len1 = (e2-e1)*0.5;
			if( len0 > rp->evlen0[k] ){ len0 = rp->evlen0[k]; }
			if( len1 > rp->evlen1[k] ){ len1 = rp->evlen1[k]; }
			b1.set( e1, ave, 1.0, 0.0, len0, e2, alpha0[e2], dx, dy, len1, NULL );

			for( i=s2+1; i<s1; i++ ){
				alpha1[i] = b0.gety( (double)i );
			}
			for( i=s1; i<=e1; i++ ){
				alpha1[i] = ave;
			}
			for( i=e1+1; i<e2; i++ ){
				alpha1[i] = b1.gety( (double)i );
			}

			if( fp ){
				fprintf( fp, "iter %d", k );
				for( i=0; i<Xcnt; i++){
					fprintf( fp, ",%f", alpha1[i] );
				}
				fprintf( fp, "\n" );
			}
		}
	}
	memcpy( alpha0, alpha1, sizeof(double)*Xcnt );

	if( fp ){
		fprintf( fp, "final" );
		for( i=0; i<Xcnt; i++){
			fprintf( fp, ",%f", alpha0[i] );
		}
		fprintf( fp, "\n" );
	}

end:
	if( fp ){ fclose(fp); fp=NULL; }
	return ret;
}



