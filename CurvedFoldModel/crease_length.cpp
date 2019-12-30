#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "crease.h"
#include "curvedraw.h"
#include "util.h"

void crease::setRLen( double len )
{
	for( int i=Xsidx; i<Xeidx; i++ ){
		rllen[i] = rrlen[i] = len;
	}
}

int crease::calcRLen0()
{
	int i, j, ret=0;
	double rllen0[MAX_SPCNT], rrlen0[MAX_SPCNT];

	// CP上の交差判定をもとにcrease長さを決める
	memcpy( rllen0, rllen, sizeof(double)*MAX_SPCNT );
	memcpy( rrlen0, rrlen, sizeof(double)*MAX_SPCNT );
	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		for( j=1+1; j<Xcnt-CEMGN; j++ ){

			double rx0,ry0, rx1,ry1, lx0,ly0, lx1,ly1, denom, dx,dy, l0,l1;
			dx = Xx2d[j] - Xx2d[i];
			dy = Xy2d[j] - Xy2d[i];

			rx0 = rrx_cp[i];
			ry0 = rry_cp[i];
			rx1 = rrx_cp[j];
			ry1 = rry_cp[j];
			denom = rx0*ry1 - rx1*ry0;
			l0 = (ry1*dx - rx1*dy) / denom;
			l1 = (ry0*dx - rx0*dy) / denom;
			if( l0>0 && rrlen[i] > l0 && l1>0 && rrlen[i+1] > l1){
				if( rrlen0[i] > l0){ rrlen0[i] = l0; }
				if( rrlen0[j] > l1){ rrlen0[j] = l1; }
			}

			lx0 = rlx_cp[i];
			ly0 = rly_cp[i];
			lx1 = rlx_cp[j];
			ly1 = rly_cp[j];
			denom = lx0*ly1 - lx1*ly0;
			l0 = (ly1*dx - lx1*dy) / denom;
			l1 = (ly0*dx - lx0*dy) / denom;
			if( l0>0 && rllen[i] > l0 && l1>0 && rllen[i+1] > l1){
				if( rllen0[i] > l0){ rllen0[i] = l0; }
				if( rllen0[j] > l1){ rllen0[j] = l1; }
			}
		}
	}
	memcpy( rllen, rllen0, sizeof(double)*MAX_SPCNT );
	memcpy( rrlen, rrlen0, sizeof(double)*MAX_SPCNT );
end:
	return ret;
}

int crease::calcRLen1_( double *rx_cp, double *ry_cp, double *rlen, char *logfile=NULL )
{
	int i, j, ret=0;
	bool flg[MAX_SPCNT], fmat[MAX_SPCNT*MAX_SPCNT];
	double *lmat=NULL, *len0=NULL, *len1=NULL, dl=0.5, proclen=0.0;
	FILE *fp=NULL;

	if( logfile ){
		fp = fopen( logfile, "w" );
	}

	// CP上の交差判定をもとにcrease長さを決める
	lmat = new double[Xcnt*Xcnt];
	len0 = new double[Xcnt];
	len1 = new double[Xcnt];
	memset( flg, 0, sizeof(bool)*MAX_SPCNT );
	memset( fmat, 0, sizeof(bool)*MAX_SPCNT*MAX_SPCNT );
	memset( lmat, 0, sizeof(double)*Xcnt*Xcnt );
	memcpy( len0, rlen, sizeof(double)*Xcnt );

	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		for( j=1+1; j<Xcnt-CEMGN; j++ ){

			double rx0,ry0, rx1,ry1, denom, dx,dy, l0,l1;
			dx = Xx2d[j] - Xx2d[i];
			dy = Xy2d[j] - Xy2d[i];

			rx0 = rx_cp[i];
			ry0 = ry_cp[i];
			rx1 = rx_cp[j];
			ry1 = ry_cp[j];
			denom = rx0*ry1 - rx1*ry0;
			l0 = lmat[j+i*Xcnt] = (ry1*dx - rx1*dy) / denom;
			l1 = lmat[i+j*Xcnt] = (ry0*dx - rx0*dy) / denom;
			if( l0>0 && rrlen[i] > l0 && l1>0 && rrlen[i+1] > l1){
				if( len0[i] > l0){ len0[i] = l0; }
				if( len0[j] > l1){ len0[j] = l1; }
				fmat[i+j*MAX_SPCNT] = 1;
			}
		}
	}
	if( fp ){
		fprintf(fp, "lmat\n" );
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			for( j=CEMGN; j<Xcnt-CEMGN; j++ ){
				fprintf(fp, ",%f", lmat[j+i*Xcnt] );
			}
			fprintf(fp, "\n" );
		}
		fprintf(fp, "fmat\n" );
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			for( j=CEMGN; j<Xcnt-CEMGN; j++ ){
				fprintf(fp, ",%d", (int)fmat[i+j*MAX_SPCNT] );
			}
			fprintf(fp, "\n" );
		}
	}

	// 他と交差しないrulingは長さ確定 flg[i]=0
	proclen = maxrlen;
	for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
		int flgi=0;
		for( j=CEMGN; j<Xcnt-CEMGN; j++ ){
			flgi += (int)fmat[i+j*MAX_SPCNT];
		}
		if( flgi > 0 ){
			flg[i] = 1; // 交差あり
		}
		if( proclen > len0[i] ){
			proclen = len0[i];
		}
	}

	// 他のrulingも長さ制限された場合に交差しない最大長さ
	for( int k=0; k<Xcnt && proclen < maxrlen; k++, proclen += dl ){

		// すべて長さ確定？ -> break
		int flgsum=0;
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			flgsum += (int)flg[i];
		}
		if( flgsum==0 ){
			//fprintf(fp, "break k=%d\n", k);
			break; // k
		}

		// 長さ確定していないrulingを少し伸ばす
		memcpy( len1, len0, sizeof(double)*Xcnt );
		//fprintf(fp, "%d,%.2f,", k, proclen);
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			if( flg[i] && len0[i] < maxrlen){
				len1[i] = len0[i] > proclen ? len0[i] : proclen;
			}
			//fprintf(fp, "%.2f,", len1[i]);
		}
		//fprintf(fp, "\n");

		// 交差判定、交差したら長さ確定
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			for( j=i+1; j<Xcnt-CEMGN; j++ ){
				int crossing = 0;
				//if( fmat[i+j*MAX_SPCNT] && ( flg[i] || flg[j] )){
				if( fmat[i+j*MAX_SPCNT] ){
					double l0,l1;
					l0 = lmat[j+i*Xcnt];
					l1 = lmat[i+j*Xcnt];
					if( l0>0 && len1[i] > l0 && l1>0 && len1[j] > l1){
						crossing = 1;
					}
				}
				if( crossing ){
					flg[i] = flg[j] = 0;
				}
			}
		}
		for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
			if( flg[i] ){ len0[i] = len1[i]; }
		}
		if( fp ){
			fprintf(fp, "k=%d\n", k );
			for( i=CEMGN; i<Xcnt-CEMGN; i++ ){
				fprintf(fp, ",%f", len0[i] );
			}
			fprintf(fp, "\n" );
		}
	} // k
	memcpy( rlen, len0, sizeof(double)*Xcnt );
	//fclose(fp); fp=NULL;
end:
	if( fp ){ fclose(fp); fp=NULL; }
	if( lmat ){ delete[] lmat; lmat=NULL; }
	if( len0 ){ delete[] len0; len0=NULL; }
	return ret;
}


int crease::calcRLen1()
{
	int ret=0;
	ret = calcRLen1_( rrx_cp, rry_cp, rrlen, NULL/*"output/rrlen.csv"*/ );
	ret = calcRLen1_( rlx_cp, rly_cp, rllen, NULL/*"output/rllen.csv"*/ );
end:
	return ret;
}

int crease::calcRLenP_( double *cvx, double *cvy, double *rx, double *ry, int cvcnt,
					   int psx, int psy, int pex, int pey,
					   double *rlen, ruling_end_type *rflg, char *logfile=NULL )
{
	int ret=0;
	FILE *fp=NULL;
#if 1
	if( logfile ){
		fp = fopen( logfile, "w" );
	}
	if( fp ){
		for( int i=0; i<cvcnt; i++ ){
			fprintf( fp, "%d,%f,%f\n", i, cvx[i], cvy[i] );
		}
		fprintf( fp, ",%d,%d\n", psx,psy );
		fprintf( fp, ",%d,%d\n", psx,pey );
		fprintf( fp, ",%d,%d\n", pex,pey );
		fprintf( fp, ",%d,%d\n", pex,psy );
	}
#endif

	for( int i=0; i<cvcnt; i++ ){
		// 始点(Xx2d[i],Xy2d[i])が範囲内かチェック
		if( cvx[i] < psx || pex <= cvx[i] || cvy[i] < psy || pey < cvy[i] || rx[i]==0.0 || ry[i]==0.0 )
		{
			rlen[i] = 0.0;
			continue;
		}
		// rulingと4辺の交差する点を算出
		// ( Xx + rx*len, Xy + ry*len )
		// Xx + rx*len = psx  ->  len = (psx - Xx) / rx
		// Xy + ry*len = psy  ->  len = (psy - Xy) / ry
		double len[4];
		len[0] = (psy - cvy[i]) / ry[i]; // top
		len[1] = (pex - cvx[i]) / rx[i]; // right
		len[2] = (pey - cvy[i]) / ry[i]; // bottom
		len[3] = (psx - cvx[i]) / rx[i]; // left

		// 長さが正かつ最小のものを選択
		int minj=-1;
		double minlen=1000;
		for( int j=0; j<4; j++ ){
			if( len[j] >= 0 && minlen > len[j] ){
				minj = j;
				minlen = len[j];
			}
		}
		if( rflg[i] == RETYPE_UNDEF || rlen[i] > minlen ){
			rlen[i] = minlen;
			switch(minj){
				case 0:	rflg[i] = RETYPE_EGTOP;		break;
				case 1:	rflg[i] = RETYPE_EGRIGHT;	break;
				case 2:	rflg[i] = RETYPE_EGBOTTOM;	break;
				case 3:	rflg[i] = RETYPE_EGLEFT;	break;
			}
		}
	}
end:
	if( fp ){ fclose(fp); fp==NULL; }
	return ret;
}

int crease::calcRLenP( int psx, int psy, int pex, int pey ) // crease長さを紙の端まで
{
	int ret=0;
	// r*flg, r*len 更新
	ret = calcRLenP_( &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
		&(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), Xeidx-Xsidx+1,
		psx,psy,pex,pey,
		&(rrlen[Xsidx]), &(rrflg[Xsidx]), NULL /*"output/rrlen.csv"*/ );
	ret = calcRLenP_( &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
		&(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), Xeidx-Xsidx+1,
		psx,psy,pex,pey,
		&(rllen[Xsidx]), &(rlflg[Xsidx]), NULL /*"output/rrlen.csv"*/ );
end:
	return ret;
}

int crease::calcRLenC_( double *cvx, double *cvy, double *rx, double *ry, int cvcnt,
					   double *cvx1, double *cvy1, int cvcnt1, double *rlen, ruling_end_type *rflg )
{
	int si=-1, ei=-1, ret=0;

	// 1回しかcrossしないと仮定した簡易版
	for( int i=0; i<cvcnt-1; i++ ){
		double ix,iy,l0,l1;
		for( int j=0; j<cvcnt1-1; j++ ){
			int ret = intersectionOfLine( cvx[i], cvy[i], cvx[i+1], cvy[i+1],
				cvx1[j], cvy1[j], cvx1[j+1], cvy1[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			if( si == -1 ){
				si = i;
			} else {
				ei = i;
			}
			break;
		}
	}
	if(si>cvcnt*0.5){
		ei=si; si=-1;
	}
	if(si>-1){
		for( int i=0; i<=si; i++ ){
			rflg[i] = RETYPE_UNDEF;
			rlen[i] = 0.0;
		}
	}
	if(ei>-1){
		for( int i=ei+1; i<cvcnt; i++ ){
			rflg[i] = RETYPE_UNDEF;
			rlen[i] = 0.0;
		}
	}

	for( int i=0; i<cvcnt; i++ ){
		if( rx[i]==0.0 || ry[i]==0.0 ){
			rlen[i] = 0.0;
			continue;
		}
		double x0, y0, x1, y1, ix,iy,l0,l1;
		x0 = cvx[i];
		y0 = cvy[i];
		x1 = cvx[i]+rx[i]*rlen[i];
		y1 = cvy[i]+ry[i]*rlen[i];
		for( int j=0; j<cvcnt1-1; j++ ){
			int ret = intersectionOfLine( x0, y0, x1, y1, cvx1[j], cvy1[j], cvx1[j+1], cvy1[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			if( rlen[i] > l0 ){
				rflg[i] = RETYPE_TRIM;
				rlen[i] = l0;
			}
		}
	}

end:
	return ret;
}
#if 1
int crease::calcRLenC( int dccnt, void *_dc )
{
	int ret=0;
	curvedraw *dc = (curvedraw *)_dc;
	// TODO: r*flg, r*len 更新
	for( int i=0; i<dccnt; i++ ){
		if( dc[i].ctype!=CTYPE_TRIM ){
			continue;
		}
		ret = calcRLenC_( &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
			&(rrx_cp[Xsidx]), &(rry_cp[Xsidx]), Xeidx-Xsidx+1,
			//cvx0[j], cvy0[j], cvcnt0[j],
			dc[i].cvx, dc[i].cvy, dc[i].cvcnt, 
			&(rrlen[Xsidx]), &(rrflg[Xsidx]) );
		ret = calcRLenC_( &(Xx2d[Xsidx]), &(Xy2d[Xsidx]),
			&(rlx_cp[Xsidx]), &(rly_cp[Xsidx]), Xeidx-Xsidx+1,
			//cvx0[j], cvy0[j], cvcnt0[j],
			dc[i].cvx, dc[i].cvy, dc[i].cvcnt, 
			&(rllen[Xsidx]), &(rlflg[Xsidx]) );
	}
end:
	return ret;
}
#endif
#if 0
int crease::calcRLenC0()
{
	int i, j, ret=0;
	for( i=0; i<ccnt1; i++ ){
		if( ctype1[i] != CTYPE_TRIM ){	// ctype1 = -1:undefined, 0:trim, 1: fold
			continue;
		}
		for( j=0; j<cvcnt1[i]; j++){
			if( cvtype[i][j]==CVTYPE_RUL_LEFT ){	// cvtype=0:rl, 1:rr, 2:曲線始点, 3:曲線終点, 4:紙端
				int idx = cvidx[i][j];
				double len = cvlen[i][j];
				if( rllen[idx] > len ){
					rllen[idx] = len;
					rlflg[idx] = RETYPE_CURVE;
				}
			} else if( cvtype[i][j]==CVTYPE_RUL_RIGHT ) {	// cvtype=0:rl, 1:rr, 2:曲線始点, 3:曲線終点, 4:紙端
				int idx = cvidx[i][j];
				double len = cvlen[i][j];
				if( rrlen[idx] > len ){
					rrlen[idx] = len;
					rrflg[idx] = RETYPE_CURVE;
				}
			}
		}
	}
end:
	return ret;
}
#endif
// 曲率小さい部分は長さ0に
int crease::calcRLenHoseiR( double rectifyR_kvthres )
{
	int ret=0;
	ret = calcRLenHoseiR( Xcnt, k2d, rectifyR_kvthres, rrlen, rrflg );
	ret = calcRLenHoseiR( Xcnt, k2d, rectifyR_kvthres, rllen, rlflg );
	return ret;
}

int crease::calcRLenHoseiR( int Xcnt, double *k2d, double rectifyR_kvthres, double *rlen, ruling_end_type *rflg )
{
	int ret=0;
	for( int i=0; i<Xcnt; i++ ){
		if( fabs(k2d[i]) < rectifyR_kvthres ){
			rlen[i]=0;
		}
	}
	return ret;
}

int crease::setEnds( int psx, int psy, int pex, int pey, int dccnt, void *_dcurve )
{
	crease *c = this;
	curvedraw *dcurve = NULL;
	if( _dcurve ){ dcurve = (curvedraw *)_dcurve; }
	// 初期値
	//c->Xsidx = c->Xsidx > _Xsidx ? c->Xsidx : _Xsidx; // 内側をとる
	//c->Xeidx = c->Xeidx < _Xeidx ? c->Xeidx : _Xeidx; // 内側をとる
	//if( c->rl==0 ){
	c->Xxs2d = c->Xx2d[c->Xsidx];
	c->Xys2d = c->Xy2d[c->Xsidx];
	c->Xxs = c->Xx[c->Xsidx];
	c->Xys = c->Xy[c->Xsidx];
	c->Xzs = c->Xz[c->Xsidx];
	c->Xxe2d = c->Xx2d[c->Xeidx];
	c->Xye2d = c->Xy2d[c->Xeidx];
	c->Xxe = c->Xx[c->Xeidx];
	c->Xye = c->Xy[c->Xeidx];
	c->Xze = c->Xz[c->Xeidx];
	//}
	int Xsidx_P=c->Xsidx, Xeidx_P=c->Xeidx;
	double Xsofs_P=-1, Xeofs_P=-1;
	crease_end_type Xstype_P=CETYPE_UNDEF, Xetype_P=CETYPE_UNDEF;
	double dx2d_Ps, dy2d_Ps, dx3d_Ps, dy3d_Ps, dz3d_Ps;
	double dx2d_Pe, dy2d_Pe, dx3d_Pe, dy3d_Pe, dz3d_Pe;

	int Xsidx_C=c->Xsidx, Xeidx_C=c->Xeidx;
	double Xsofs_C=-1, Xeofs_C=-1;
	double dx2d_Cs, dy2d_Cs, dx3d_Cs, dy3d_Cs, dz3d_Cs;
	double dx2d_Ce, dy2d_Ce, dx3d_Ce, dy3d_Ce, dz3d_Ce;

	// 
	// 紙端との交点（始点）-> Xsidx_P, Xsofs_P, Xstype_P -> Xxs2d, Xys2d
	//
	for( int i=c->Xsidx; i<c->Xeidx; i++ ){

		// 紙の範囲内になるまでたどる
		if( c->Xx2d[i] < psx || pex <= c->Xx2d[i] || c->Xy2d[i] < psy || pey < c->Xy2d[i] ){
			continue;
		}

		// 始点が紙の範囲内
		if( i==c->Xsidx ){
			Xstype_P = CETYPE_EGTOP1; // 仮
			break;
		}

		int minj=-1;
		double d, len[4], minlen;
		dx2d_Ps = c->Xx2d[i-1] - c->Xx2d[i];
		dy2d_Ps = c->Xy2d[i-1] - c->Xy2d[i];
		dx3d_Ps = c->Xx[i-1] - c->Xx[i];
		dy3d_Ps = c->Xy[i-1] - c->Xy[i];
		dz3d_Ps = c->Xz[i-1] - c->Xz[i];
		minlen = d = normalize_v2( &dx2d_Ps, &dy2d_Ps );
		normalize_v3( &dx3d_Ps, &dy3d_Ps, &dz3d_Ps );
		len[0] = (psy - c->Xy2d[i]) / dy2d_Ps;
		len[1] = (pex - c->Xx2d[i]) / dx2d_Ps;
		len[2] = (pey - c->Xy2d[i]) / dy2d_Ps;
		len[3] = (psx - c->Xx2d[i]) / dx2d_Ps;
		for( int j=0; j<4; j++ ){
			if( len[j]<0 ){
				continue;
			}
			if( minlen > len[j] ){
				minj = j;
				minlen = len[j];
			}
		}
		if( minj>-1 ){
			Xsidx_P = i;
			Xsofs_P = minlen;
			switch( minj ){
				case 0:	Xstype_P = CETYPE_EGTOP0;	break;
				case 1:	Xstype_P = CETYPE_EGRIGHT0;	break;
				case 2:	Xstype_P = CETYPE_EGBOTTOM0;	break;
				case 3:	Xstype_P = CETYPE_EGLEFT0;	break;
			}
			break;
		}
		break; // ここまで来ることはないはず
	}
	if( Xstype_P == CETYPE_UNDEF ){ // 曲線全体が紙の外
		goto end;
	}
	if( Xstype_P == CETYPE_EGTOP1 ) { // 始点が紙の範囲内 -> 紙の端との交点 X*s2d, X*s

		// 延長する方向
		dx2d_Ps = c->Xx2d[c->Xsidx] - c->Xx2d[c->Xsidx+1];
		dy2d_Ps = c->Xy2d[c->Xsidx] - c->Xy2d[c->Xsidx+1];
		dx3d_Ps = c->Xx[c->Xsidx] - c->Xx[c->Xsidx+1];
		dy3d_Ps = c->Xy[c->Xsidx] - c->Xy[c->Xsidx+1];
		dz3d_Ps = c->Xz[c->Xsidx] - c->Xz[c->Xsidx+1];
		normalize_v2( &dx2d_Ps, &dy2d_Ps );
		normalize_v3( &dx3d_Ps, &dy3d_Ps, &dz3d_Ps );

		// 延長する長さ
		// ( Xx + edx*len, Xy + edy*len )
		// Xx + edx*len = psx  ->  len = (psx - Xx) / edx
		// Xy + edy*len = psy  ->  len = (psy - Xy) / edy
		double len[4];
		len[0] = (psy - c->Xy2d[c->Xsidx]) / dy2d_Ps;
		len[1] = (pex - c->Xx2d[c->Xsidx]) / dx2d_Ps;
		len[2] = (pey - c->Xy2d[c->Xsidx]) / dy2d_Ps;
		len[3] = (psx - c->Xx2d[c->Xsidx]) / dx2d_Ps;

		// 長さが正かつ最小のものを選択
		int minj=-1;
		double minlen=1000;
		for( int j=0; j<4; j++ ){
			if( len[j] > 0 && minlen > len[j] ){
				minj = j;
				minlen = len[j];
			}
		}
		if( minj>-1 ){
			Xsidx_P = c->Xsidx;
			Xsofs_P = minlen;
			switch( minj ){
				case 0:	Xstype_P = CETYPE_EGTOP1;	break;
				case 1:	Xstype_P = CETYPE_EGRIGHT1;	break;
				case 2:	Xstype_P = CETYPE_EGBOTTOM1;	break;
				case 3:	Xstype_P = CETYPE_EGLEFT1;	break;
			}
		}
	}
	// 始点
	if( c->rl==0 ){
		c->Xxs2d = c->Xx2d[Xsidx_P] + Xsofs_P * dx2d_Ps;
		c->Xys2d = c->Xy2d[Xsidx_P] + Xsofs_P * dy2d_Ps;
		c->Xxs = c->Xx[Xsidx_P] + Xsofs_P * dx3d_Ps;
		c->Xys = c->Xy[Xsidx_P] + Xsofs_P * dy3d_Ps;
		c->Xzs = c->Xz[Xsidx_P] + Xsofs_P * dz3d_Ps;
	}
#if 0
	else if((Xstype_P==CETYPE_EGTOP0 || Xstype_P==CETYPE_EGTOP1) && c->Xys2d0==psy
		|| (Xstype_P==CETYPE_EGBOTTOM0 || Xstype_P==CETYPE_EGBOTTOM1) && c->Xys2d0==pey
		|| (Xstype_P==CETYPE_EGRIGHT0 || Xstype_P==CETYPE_EGRIGHT1) && c->Xxs2d0==pex
		|| (Xstype_P==CETYPE_EGLEFT0 || Xstype_P==CETYPE_EGLEFT1) && c->Xxs2d0==psx)
	{
		c->Xxs2d = c->Xxs2d0;
		c->Xys2d = c->Xys2d0;
		//c->Xxs = c->Xys = c->Xzs = 0;
		c->Xxs = c->Xx[Xsidx_P] + Xsofs_P * dx3d_Ps;
		c->Xys = c->Xy[Xsidx_P] + Xsofs_P * dy3d_Ps;
		c->Xzs = c->Xz[Xsidx_P] + Xsofs_P * dz3d_Ps;
	} else {
		c->Xxs2d = c->Xxe2d0;
		c->Xys2d = c->Xye2d0;
		//c->Xxs = c->Xys = c->Xzs = 0;
		c->Xxs = c->Xx[Xeidx_P] + Xeofs_P * dx3d_Pe;
		c->Xys = c->Xy[Xeidx_P] + Xeofs_P * dy3d_Pe;
		c->Xzs = c->Xz[Xeidx_P] + Xeofs_P * dz3d_Pe;
	}
#endif
	// 
	// 紙端との交点（終点）-> Xeidx_P, Xeofs_P, Xetype_P -> Xxe2d, Xye2d
	//
	for( int i=c->Xeidx; i>c->Xsidx; i-- ){

		// 紙の範囲内になるまでたどる
		if( c->Xx2d[i] < psx || pex <= c->Xx2d[i] || c->Xy2d[i] < psy || pey < c->Xy2d[i] ){
			continue;
		}

		// 終点が紙の範囲内
		if( i==c->Xeidx ){
			Xetype_P = CETYPE_EGTOP1; // 仮
			break;
		}

		int minj=-1;
		double d, len[4], minlen;
		dx2d_Pe = c->Xx2d[i+1] - c->Xx2d[i];
		dy2d_Pe = c->Xy2d[i+1] - c->Xy2d[i];
		dx3d_Pe = c->Xx[i+1] - c->Xx[i];
		dy3d_Pe = c->Xy[i+1] - c->Xy[i];
		dz3d_Pe = c->Xz[i+1] - c->Xz[i];
		minlen = d = normalize_v2( &dx2d_Pe, &dy2d_Pe );
		normalize_v3( &dx3d_Pe, &dy3d_Pe, &dz3d_Pe );
		len[0] = (psy - c->Xy2d[i]) / dy2d_Pe;
		len[1] = (pex - c->Xx2d[i]) / dx2d_Pe;
		len[2] = (pey - c->Xy2d[i]) / dy2d_Pe;
		len[3] = (psx - c->Xx2d[i]) / dx2d_Pe;
		for( int j=0; j<4; j++ ){
			if( len[j]<0 ){
				continue;
			}
			if( minlen > len[j] ){
				minj = j;
				minlen = len[j];
			}
		}
		if( minj>-1 ){
			Xeidx_P = i;
			Xeofs_P = minlen;
			switch( minj ){
				case 0:	Xetype_P = CETYPE_EGTOP0;	break;
				case 1:	Xetype_P = CETYPE_EGRIGHT0;	break;
				case 2:	Xetype_P = CETYPE_EGBOTTOM0;	break;
				case 3:	Xetype_P = CETYPE_EGLEFT0;	break;
			}
			break;
		}
		break; // ここまで来ることはないはず
	}
	if( Xetype_P == CETYPE_UNDEF ){ // 曲線全体が紙の外
		goto end;
	}
	if( Xetype_P == CETYPE_EGTOP1 ) { // 始点が紙の範囲内 -> 紙の端との交点 X*e2d, X*e

		// 延長する方向
		dx2d_Pe = c->Xx2d[c->Xeidx] - c->Xx2d[c->Xeidx-1];
		dy2d_Pe = c->Xy2d[c->Xeidx] - c->Xy2d[c->Xeidx-1];
		dx3d_Pe = c->Xx[c->Xeidx] - c->Xx[c->Xeidx-1];
		dy3d_Pe = c->Xy[c->Xeidx] - c->Xy[c->Xeidx-1];
		dz3d_Pe = c->Xz[c->Xeidx] - c->Xz[c->Xeidx-1];
		normalize_v2( &dx2d_Pe, &dy2d_Pe );
		normalize_v3( &dx3d_Pe, &dy3d_Pe, &dz3d_Pe );

		// 延長する長さ
		// ( Xx + edx*len, Xy + edy*len )
		// Xx + edx*len = psx  ->  len = (psx - Xx) / edx
		// Xy + edy*len = psy  ->  len = (psy - Xy) / edy
		double len[4];
		len[0] = (psy - c->Xy2d[c->Xeidx]) / dy2d_Pe;
		len[1] = (pex - c->Xx2d[c->Xeidx]) / dx2d_Pe;
		len[2] = (pey - c->Xy2d[c->Xeidx]) / dy2d_Pe;
		len[3] = (psx - c->Xx2d[c->Xeidx]) / dx2d_Pe;

		// 長さが正かつ最小のものを選択
		int minj=-1;
		double minlen=1000;
		for( int j=0; j<4; j++ ){
			if( len[j] > 0 && minlen > len[j] ){
				minj = j;
				minlen = len[j];
			}
		}
		if( minj>-1 ){
			Xeidx_P = c->Xeidx;
			Xeofs_P = minlen;
			switch( minj ){
				case 0:	Xetype_P = CETYPE_EGTOP1;	break;
				case 1:	Xetype_P = CETYPE_EGRIGHT1;	break;
				case 2:	Xetype_P = CETYPE_EGBOTTOM1;	break;
				case 3:	Xetype_P = CETYPE_EGLEFT1;	break;
			}
		}
	}
	// 終点
	if( c->rl==0 ){
		c->Xxe2d = c->Xx2d[Xeidx_P] + Xeofs_P * dx2d_Pe;
		c->Xye2d = c->Xy2d[Xeidx_P] + Xeofs_P * dy2d_Pe;
		c->Xxe = c->Xx[Xeidx_P] + Xeofs_P * dx3d_Pe;
		c->Xye = c->Xy[Xeidx_P] + Xeofs_P * dy3d_Pe;
		c->Xze = c->Xz[Xeidx_P] + Xeofs_P * dz3d_Pe;
	}
#if 0
	else if((Xetype_P==CETYPE_EGTOP0 || Xetype_P==CETYPE_EGTOP1) && c->Xye2d0==psy
		|| (Xetype_P==CETYPE_EGBOTTOM0 || Xetype_P==CETYPE_EGBOTTOM1) && c->Xye2d0==pey
		|| (Xetype_P==CETYPE_EGRIGHT0 || Xetype_P==CETYPE_EGRIGHT1) && c->Xxe2d0==pex
		|| (Xetype_P==CETYPE_EGLEFT0 || Xetype_P==CETYPE_EGLEFT1) && c->Xxe2d0==psx)	
	{
		c->Xxe2d = c->Xxe2d0;
		c->Xye2d = c->Xye2d0;
		//c->Xxe = c->Xye = c->Xze = 0;
		c->Xxe = c->Xx[Xeidx_P] + Xeofs_P * dx3d_Pe;
		c->Xye = c->Xy[Xeidx_P] + Xeofs_P * dy3d_Pe;
		c->Xze = c->Xz[Xeidx_P] + Xeofs_P * dz3d_Pe;
	} else {
		c->Xxe2d = c->Xxs2d0;
		c->Xye2d = c->Xys2d0;
		//c->Xxe = c->Xye = c->Xze = 0;
		c->Xxe = c->Xx[Xsidx_P] + Xsofs_P * dx3d_Ps;
		c->Xye = c->Xy[Xsidx_P] + Xsofs_P * dy3d_Ps;
		c->Xze = c->Xz[Xsidx_P] + Xsofs_P * dz3d_Ps;
	}
#endif
	if( c->rl!=0 ){

		if((Xstype_P==CETYPE_EGTOP0 || Xstype_P==CETYPE_EGTOP1) && c->Xys2d0==psy
			|| (Xstype_P==CETYPE_EGBOTTOM0 || Xstype_P==CETYPE_EGBOTTOM1) && c->Xys2d0==pey
			|| (Xstype_P==CETYPE_EGRIGHT0 || Xstype_P==CETYPE_EGRIGHT1) && c->Xxs2d0==pex
			|| (Xstype_P==CETYPE_EGLEFT0 || Xstype_P==CETYPE_EGLEFT1) && c->Xxs2d0==psx)
		{
			c->Xxs2d = c->Xxs2d0;
			c->Xys2d = c->Xys2d0;
			c->Xxs = c->Xys = c->Xzs = 0;
			//c->Xxs = c->Xx[Xsidx_P] + Xsofs_P * dx3d_Ps;
			//c->Xys = c->Xy[Xsidx_P] + Xsofs_P * dy3d_Ps;
			//c->Xzs = c->Xz[Xsidx_P] + Xsofs_P * dz3d_Ps;
		} else {
			c->Xxs2d = c->Xxe2d0;
			c->Xys2d = c->Xye2d0;
			c->Xxs = c->Xys = c->Xzs = 0;
			//c->Xxs = c->Xx[Xeidx_P] + Xeofs_P * dx3d_Pe;
			//c->Xys = c->Xy[Xeidx_P] + Xeofs_P * dy3d_Pe;
			//c->Xzs = c->Xz[Xeidx_P] + Xeofs_P * dz3d_Pe;
		}
		if((Xetype_P==CETYPE_EGTOP0 || Xetype_P==CETYPE_EGTOP1) && c->Xye2d0==psy
			|| (Xetype_P==CETYPE_EGBOTTOM0 || Xetype_P==CETYPE_EGBOTTOM1) && c->Xye2d0==pey
			|| (Xetype_P==CETYPE_EGRIGHT0 || Xetype_P==CETYPE_EGRIGHT1) && c->Xxe2d0==pex
			|| (Xetype_P==CETYPE_EGLEFT0 || Xetype_P==CETYPE_EGLEFT1) && c->Xxe2d0==psx)	
		{
			c->Xxe2d = c->Xxe2d0;
			c->Xye2d = c->Xye2d0;
			c->Xxe = c->Xye = c->Xze = 0;
			//c->Xxe = c->Xx[Xeidx_P] + Xeofs_P * dx3d_Pe;
			//c->Xye = c->Xy[Xeidx_P] + Xeofs_P * dy3d_Pe;
			//c->Xze = c->Xz[Xeidx_P] + Xeofs_P * dz3d_Pe;
		} else {
			c->Xxe2d = c->Xxs2d0;
			c->Xye2d = c->Xys2d0;
			c->Xxe = c->Xye = c->Xze = 0;
			//c->Xxe = c->Xx[Xsidx_P] + Xsofs_P * dx3d_Ps;
			//c->Xye = c->Xy[Xsidx_P] + Xsofs_P * dy3d_Ps;
			//c->Xze = c->Xz[Xsidx_P] + Xsofs_P * dz3d_Ps;
		}

	}
	//
	// 曲線との交点
	//
	for( int k=0; k<dccnt; k++ ){
		if( dcurve[k].ctype!=CTYPE_TRIM ){ // 0:trim
			continue;
		}

		int is0[MAX_SPCNT], is1[MAX_SPCNT], iscnt;
		double isl0[MAX_SPCNT], isl1[MAX_SPCNT];

		iscnt=0;
		is0[iscnt]=-1;	is1[iscnt]=Xsidx_P-1;	isl0[iscnt]=isl1[iscnt]=0.0;	iscnt++;

		// (Xsidx_P-1, Xsidx_P) との交差判定
		for( int j=0; j<dcurve[k].cvcnt-1; j++ ){
			double ix,iy,l0,l1, dx, dy, dd;
			int ret = intersectionOfLine( c->Xxs2d, c->Xys2d, c->Xx2d[Xsidx_P], c->Xy2d[Xsidx_P],
				dcurve[k].cvx[j], dcurve[k].cvy[j], dcurve[k].cvx[j+1], dcurve[k].cvy[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			dx = c->Xx2d[Xsidx_P] - c->Xxs2d;
			dy = c->Xy2d[Xsidx_P] - c->Xys2d;
			dd = sqrt( dx*dx + dy*dy );
			is0[iscnt] = Xsidx_P-1;
			is1[iscnt] = Xsidx_P;
			isl0[iscnt] = l0;
			isl1[iscnt] = dd-l0;
			iscnt++;
			break;
		}
		// (i, i+1) との交差判定
		for( int i=Xsidx_P; i<Xeidx_P; i++ ){
			for( int j=0; j<dcurve[k].cvcnt-1; j++ ){
				double ix,iy,l0,l1, dx, dy, dd;
				int ret = intersectionOfLine( c->Xx2d[i], c->Xy2d[i], c->Xx2d[i+1], c->Xy2d[i+1],
					dcurve[k].cvx[j], dcurve[k].cvy[j], dcurve[k].cvx[j+1], dcurve[k].cvy[j+1], &ix, &iy, &l0, &l1 );
				if( ret<0 ){
					continue;
				}
				dx = c->Xx2d[i+1] - c->Xx2d[i];
				dy = c->Xy2d[i+1] - c->Xy2d[i];
				dd = sqrt( dx*dx + dy*dy );
				is0[iscnt] = i;
				is1[iscnt] = i+1;
				isl0[iscnt] = l0;
				isl1[iscnt] = dd-l0;
				iscnt++;
				break;
			}
		}
		// (Xeidx_P, Xeidx_P+1) との交差判定
		for( int j=0; j<dcurve[k].cvcnt-1; j++ ){
			double ix,iy,l0,l1, dx, dy, dd;
			int ret = intersectionOfLine( c->Xx2d[Xeidx_P], c->Xy2d[Xeidx_P], c->Xxe2d, c->Xye2d,
				dcurve[k].cvx[j], dcurve[k].cvy[j], dcurve[k].cvx[j+1], dcurve[k].cvy[j+1], &ix, &iy, &l0, &l1 );
			if( ret<0 ){
				continue;
			}
			dx = c->Xxe2d - c->Xx2d[Xeidx_P];
			dy = c->Xye2d - c->Xy2d[Xeidx_P];
			dd = sqrt( dx*dx + dy*dy );
			is0[iscnt] = Xeidx_P;
			is1[iscnt] = Xeidx_P+1;
			isl0[iscnt] = l0;
			isl1[iscnt] = dd-l0;
			iscnt++;
			break;
		}
		is0[iscnt]=Xeidx_P+1;	is1[iscnt]=-1;	isl0[iscnt]=isl1[iscnt]=0.0;	iscnt++;

		// 曲線との交差で区切られた区間is*のうち、最も長い区間を採用
		int maxi=-1, maxrcnt=0;
		for( int i=0; i<iscnt-1; i++ ){
			int rcnt;
			if( is0[i]<0 ){
				rcnt = is0[i+1]-(c->Xsidx-2); // -2 ?
			} else {
				rcnt = is0[i+1]-is0[i];
			}
			if( maxrcnt < rcnt ){
				maxi = i;
				maxrcnt = rcnt;
			}
		}
		if( maxi>-1 ){
			double d=0;
			if( Xsidx_C < is1[maxi] || ( Xsidx_C == is1[maxi] && ( Xsofs_C < 0.0 || isl1[maxi] < Xsofs_C ))){
				Xsidx_C = is1[maxi];
				Xsofs_C = isl1[maxi];
				dx2d_Cs = c->Xx2d[Xsidx_C-1] - c->Xx2d[Xsidx_C];
				dy2d_Cs = c->Xy2d[Xsidx_C-1] - c->Xy2d[Xsidx_C];
				normalize_v2( &dx2d_Cs, &dy2d_Cs );
				dx3d_Cs = c->Xx[Xsidx_C-1] - c->Xx[Xsidx_C];
				dy3d_Cs = c->Xy[Xsidx_C-1] - c->Xy[Xsidx_C];
				dz3d_Cs = c->Xz[Xsidx_C-1] - c->Xz[Xsidx_C];
				normalize_v3( &dx3d_Cs, &dy3d_Cs, &dz3d_Cs );
			}
			if( Xeidx_C > is0[maxi+1] || ( Xeidx_C == is0[maxi+1] && ( Xeofs_C < 0.0 || isl0[maxi+1] < Xeofs_C ))){
				Xeidx_C = is0[maxi+1];
				Xeofs_C = isl0[maxi+1];
				dx2d_Ce = c->Xx2d[Xeidx_C+1] - c->Xx2d[Xeidx_C];
				dy2d_Ce = c->Xy2d[Xeidx_C+1] - c->Xy2d[Xeidx_C];
				normalize_v2( &dx2d_Ce, &dy2d_Ce );
				dx3d_Ce = c->Xx[Xeidx_C+1] - c->Xx[Xeidx_C];
				dy3d_Ce = c->Xy[Xeidx_C+1] - c->Xy[Xeidx_C];
				dz3d_Ce = c->Xz[Xeidx_C+1] - c->Xz[Xeidx_C];
				normalize_v3( &dx3d_Ce, &dy3d_Ce, &dz3d_Ce );
			}
		}
	} // k

	// 始点
	if( Xsidx_P > Xsidx_C || ( Xsidx_P == Xsidx_C && ( Xsofs_C < 0.0 || Xsofs_P < Xsofs_C ))){
#if 0 // 更新済み
		c->Xxs2d = c->Xx2d[c->Xsidx] + Xsofs_P * dx2d_Ps;
		c->Xys2d = c->Xy2d[c->Xsidx] + Xsofs_P * dy2d_Ps;
		c->Xxs = c->Xx[c->Xsidx] + Xsofs_P * dx3d_Ps;
		c->Xys = c->Xy[c->Xsidx] + Xsofs_P * dy3d_Ps;
		c->Xzs = c->Xz[c->Xsidx] + Xsofs_P * dz3d_Ps;
#endif
		c->Xsidx = Xsidx_P;
		c->Xstype = Xstype_P;
	} else {
		c->Xsidx = Xsidx_C;
		c->Xxs2d = c->Xx2d[c->Xsidx] + Xsofs_C * dx2d_Cs;
		c->Xys2d = c->Xy2d[c->Xsidx] + Xsofs_C * dy2d_Cs;
		c->Xxs = c->Xx[c->Xsidx] + Xsofs_C * dx3d_Cs;
		c->Xys = c->Xy[c->Xsidx] + Xsofs_C * dy3d_Cs;
		c->Xzs = c->Xz[c->Xsidx] + Xsofs_C * dz3d_Cs;
		c->Xstype = CETYPE_TRIM;
	}

	// 終点
	if( Xeidx_P < Xeidx_C || ( Xeidx_P == Xeidx_C && ( Xeofs_C < 0.0 || Xeofs_P < Xeofs_C ))){
#if 0 // 更新済み
		c->Xxe2d = c->Xx2d[c->Xeidx] + Xeofs_P * dx2d_Pe;
		c->Xye2d = c->Xy2d[c->Xeidx] + Xeofs_P * dy2d_Pe;
		c->Xxe = c->Xx[c->Xeidx] + Xeofs_P * dx3d_Pe;
		c->Xye = c->Xy[c->Xeidx] + Xeofs_P * dy3d_Pe;
		c->Xze = c->Xz[c->Xeidx] + Xeofs_P * dz3d_Pe;
#endif
		c->Xeidx = Xeidx_P;
		c->Xetype = Xetype_P;
	} else {
		c->Xeidx = Xeidx_C;
		c->Xxe2d = c->Xx2d[c->Xeidx] + Xeofs_C * dx2d_Ce;
		c->Xye2d = c->Xy2d[c->Xeidx] + Xeofs_C * dy2d_Ce;
		c->Xxe = c->Xx[c->Xeidx] + Xeofs_C * dx3d_Ce;
		c->Xye = c->Xy[c->Xeidx] + Xeofs_C * dy3d_Ce;
		c->Xze = c->Xz[c->Xeidx] + Xeofs_C * dz3d_Ce;
		c->Xetype = CETYPE_TRIM;
	}

end:
	if(  Xstype_P==CETYPE_UNDEF || Xetype_P==CETYPE_UNDEF ){ // 曲線全体が紙の外

	}
	return 0;
}

void crease::initEnds()
{
	Xsidx = CEMGN;
	Xeidx = Xcnt-CEMGN;
	setRLen( 20.0 );
	Xxs2d = Xxs2d0 = Xx2d[Xsidx]; 
	Xys2d = Xys2d0 = Xy2d[Xsidx];
	Xxs = Xx[Xsidx];
	Xys = Xy[Xsidx];
	Xzs = Xz[Xsidx];
	Xxe2d = Xxe2d0 = Xx2d[Xeidx];
	Xye2d = Xye2d0 = Xy2d[Xeidx];
	Xxe = Xx[Xeidx];
	Xye = Xy[Xeidx];
	Xze = Xz[Xeidx];
}
