#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "papermodel.h"
#include "Bspline.h"

int papermodel::CP2DC()
{
	int maxcvcnt, nv2[] = {0,0,0,0, 1, 2, 3,3,3,3}, ret=0;
	double *t1=NULL;
	Bspline bs;

	maxcvcnt = 100;
#if 1
	for( int i=0; i<dccnt; i++ ){
		if( ( dcurve[i].ctype==CTYPE_FOLD || dcurve[i].ctype==CTYPE_FOLD_LEFT || dcurve[i].ctype==CTYPE_FOLD_RIGHT )
			&& maxcvcnt < dcurve[i].cvcnt ){
				maxcvcnt = dcurve[i].cvcnt;
		}
	}
#endif
	t1 = new double[maxcvcnt+10];

	for( int i=1; i<crcnt; i++ ){
		crease *c=&(crs[i]);
		int cidx = c->org_idx;

		// FOLD以外と対応付いていることはあり得ないけど一応チェック
		if( 0<=cidx && cidx < dccnt
			&& !( dcurve[cidx].ctype==CTYPE_FOLD || dcurve[cidx].ctype==CTYPE_FOLD_LEFT || dcurve[cidx].ctype==CTYPE_FOLD_RIGHT ) )
		{
			continue;
		}

		// cv0を優先する場合は更新しない
		if( c->flg_org==0 ) // 元データ： 0:org* 1:CP*
		{
			continue;
		}

		// 対応するc0がない → 末尾に加える
		if( cidx<0 || dccnt <= cidx ){
			cidx = dccnt; dccnt++;
			dcurve[cidx].cvcnt = 0;
		}

		int ccnt = dcurve[cidx].cvcnt;
		if( ccnt==0 ){
			ccnt = dcurve[cidx].cvcnt = maxcvcnt;
			if( c->rl<0 ){
				dcurve[cidx].ctype = CTYPE_FOLD_LEFT;
			} else if( c->rl>0 ){
				dcurve[cidx].ctype = CTYPE_FOLD_RIGHT;
			} else {
				dcurve[cidx].ctype = CTYPE_FOLD;
			}
		}

		ret = bs.calcT( ccnt-2, 3, 0, t1);
		ret = bs.calcBSpline( c->CPx, c->CPy, NULL, CCNT, nv2, 10, 4/*degree*/,
			ccnt-2, t1, &(dcurve[cidx].cvx[1]), &(dcurve[cidx].cvy[1]), NULL );
		dcurve[cidx].cvx[0] = dcurve[cidx].cvx[1] + (dcurve[cidx].cvx[1]-dcurve[cidx].cvx[2]);
		dcurve[cidx].cvy[0] = dcurve[cidx].cvy[1] + (dcurve[cidx].cvy[1]-dcurve[cidx].cvy[2]);
		dcurve[cidx].cvx[ccnt-1] = dcurve[cidx].cvx[ccnt-2] + (dcurve[cidx].cvx[ccnt-2]-dcurve[cidx].cvx[ccnt-3]);
		dcurve[cidx].cvy[ccnt-1] = dcurve[cidx].cvy[ccnt-2] + (dcurve[cidx].cvy[ccnt-2]-dcurve[cidx].cvy[ccnt-3]);
	}
end:
	delete[] t1;
	return ret;
}

#if 0
int papermodel::DC2FC( )
{
	int ret=0;
	int c0cnt=0;


	cvertex_type c0type[MAXDCX0];
	int c0idx[MAXDCX0], c0cidx[MAXDCX0], c0sort[MAXDCX0];
	double c0x[MAXDCX0], c0y[MAXDCX0], c0len[MAXDCX0];
	for( int i=0; i<MAXDCX0; i++ ){
		c0type[i] = CVTYPE_UNDEF;
		c0idx[i] = c0cidx[i] = -1;
		c0sort[i] = cvcnt0;
		c0x[i] = c0y[i] = c0len[i] = 0.0;
	}

	//
	// ruling－曲線交点
	//
	for( int kk=0; kk<2; kk++ ){

		double *rx=NULL, *ry=NULL;
		switch( kk ){
				case 0:	rx=rlx_cp;	ry=rly_cp;	break;	// 左側 ruling－曲線交点
				case 1:	rx=rrx_cp;	ry=rry_cp;	break;	// 右側 ruling－曲線交点
		}
		for( int i=0; i<Xcnt; i++ ){
			int intersect=0;
			double maxlen=pw+ph;
			for( int j=0; j<cvcnt0-1; j++ ){
				double ix,iy,l0,l1;
				int ret = intersectionOfLine( Xx2d[i], Xy2d[i], Xx2d[i]+rx[i]*maxlen, Xy2d[i]+ry[i]*maxlen,
					cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix, &iy, &l0, &l1 );
				if( ret<0 ){
					continue;
				}
				// rulingと曲線が複数回交差する場合は最も短いものを使う
				if( maxlen > l0 ){
					intersect = 1;
					c0cidx[c0cnt] = j;
					c0len[c0cnt] = maxlen = l0;
					c0x[c0cnt] = ix;
					c0y[c0cnt] = iy;
				}
			}
			if( intersect ){
				switch( kk ){
						case 0:	c0type[c0cnt]=CVTYPE_RUL_LEFT;	break;	// rl
						case 1:	c0type[c0cnt]=CVTYPE_RUL_RIGHT;	break;	// rr
				}
				c0idx[c0cnt] = i;
				c0cnt++;
			}
		}
	} // kk

	//
	// 紙端－曲線交点
	//
	for( int j=0; j<cvcnt0-1; j++ ){
		int ret[4];
		double ix[4],iy[4],l0[4],l1[4];
		ret[0] = intersectionOfLine( psx, psy, pex, psy, cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix[0], &iy[0], &l0[0], &l1[0] );
		ret[1] = intersectionOfLine( pex, psy, pex, pey, cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix[1], &iy[1], &l0[1], &l1[1] );
		ret[2] = intersectionOfLine( pex, pey, psx, pey, cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix[2], &iy[2], &l0[2], &l1[2] );
		ret[3] = intersectionOfLine( psx, pey, psx, psy, cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix[3], &iy[3], &l0[3], &l1[3] );
		for( int i=0; i<4; i++ ){
			if( ret[i]<0 || l1[i]==0.0 ){
				continue;
			}
			c0type[c0cnt] = CVTYPE_PAPER_EDGE; // 紙端
			c0idx[c0cnt] = i; // 0:上, 1:右, 2:下, 3:左
			c0cidx[c0cnt] = j;
			c0len[c0cnt] = l0[i];
			c0x[c0cnt] = ix[i];
			c0y[c0cnt] = iy[i];
			c0cnt++;
		}
		//break;
	}

	if( ctype0==CTYPE_TRIM ){
		//
		// 折り線(primary)－曲線交点
		//
		int is[MAXDCX0], iscidx[MAXDCX0], iscnt=0, psidx=-1, peidx=-1, pscidx=-1, pecidx=-1;
		double isx[MAXDCX0], isy[MAXDCX0], isl0[MAXDCX0];
		is[0]=0;	iscnt=1;
		for( int i=CEMGN; i<Xcnt-CEMGN-1; i++ ){
			//for( int i=Xsidx; i<Xeidx; i++ ){
			for( int j=0; j<cvcnt0-1; j++ ){
				double ix,iy,l0,l1;
				int ret = intersectionOfLine( Xx2d[i], Xy2d[i], Xx2d[i+1], Xy2d[i+1],
					cvx0[j], cvy0[j], cvx0[j+1], cvy0[j+1], &ix, &iy, &l0, &l1 );
				if( ret<0 ){
					continue;
				}
				is[iscnt] = i;
				iscidx[iscnt] = j;
				isl0[iscnt] = l0;
				isx[iscnt] = ix;
				isy[iscnt] = iy;
				iscnt++;
				break;
			}
		}
		is[iscnt]=Xcnt-1;	iscnt++;

		// 曲線との交差で区切られた区間is*のうち、最も長い区間を採用
		int maxi=-1, maxrcnt=0;// cs=-1, ce=-1;
		for( int i=0; i<iscnt-1; i++ ){
			int rcnt = is[i+1]-is[i];
			if( maxrcnt < rcnt ){
				maxi = i;
				maxrcnt = rcnt;
			}
		}
		if( is[maxi]>0){	//cs = is[i]+1;
			c0type[c0cnt] = CVTYPE_FCURVE_START;
			c0idx[c0cnt] = psidx = is[maxi];
			c0cidx[c0cnt] = pscidx = iscidx[maxi]; // j
			c0len[c0cnt] = isl0[maxi]; // l0;
			c0x[c0cnt] = isx[maxi];
			c0y[c0cnt] = isy[maxi];
			c0cnt++;
		}
		if( is[maxi+1]>-1 && is[maxi+1]<Xcnt-1 ){	//ce = is[i+1];
			c0type[c0cnt] = CVTYPE_FCURVE_END;
			c0idx[c0cnt] = peidx = is[maxi+1];
			c0cidx[c0cnt] = pecidx = iscidx[maxi+1]; // j
			c0len[c0cnt] = isl0[maxi+1]; // l0;
			c0x[c0cnt] = isx[maxi+1];
			c0y[c0cnt] = isy[maxi+1];
			c0cnt++;
		}

	} // if( ctype0==CTYPE_TRIM ){

	// 曲線更新, データコピー
	int ccnt=0;
	double cx[MAXDCX0], cy[MAXDCX0];
	for( int i=0; i<c0cnt; i++ ){
		// idx小さい方から
		int minj = -1, mincidx = cvcnt0, jj;
		double d, dx, dy;
		for( int j=0; j<c0cnt; j++ ){
			if( ccnt <= c0sort[j] && mincidx > c0cidx[j] ){
				minj = j;
				mincidx = c0cidx[j];
			}
		}
		if( minj > -1 ){
			c0sort[minj] = ccnt;
			jj = c0idx[minj];
			switch( c0type[minj] ){
				case CVTYPE_RUL_LEFT:
					cx[ccnt] = Xx2d[jj] + rlx_cp[jj]*c0len[minj]; // cl0[i]
					cy[ccnt] = Xy2d[jj] + rly_cp[jj]*c0len[minj];
					ccnt++;
					break;
				case CVTYPE_RUL_RIGHT:
					cx[ccnt] = Xx2d[jj] + rrx_cp[jj]*c0len[minj];
					cy[ccnt] = Xy2d[jj] + rry_cp[jj]*c0len[minj];
					ccnt++;
					break;
				case CVTYPE_FCURVE_START:
				case CVTYPE_FCURVE_END:
					dx = Xx2d[jj+1]-Xx2d[jj];
					dy = Xy2d[jj+1]-Xy2d[jj];
					d = sqrt(dx*dx+dy*dy);
					cx[ccnt] = Xx2d[jj] + dx/d*c0len[minj];
					cy[ccnt] = Xy2d[jj] + dy/d*c0len[minj];
					ccnt++;
					break;
				case CVTYPE_PAPER_EDGE:
					switch( c0idx[minj] ){
				case 0:	// 上
					cx[ccnt] = psx + c0len[minj];
					cy[ccnt] = psy;
					break;
				case 1:	// 右
					cx[ccnt] = pex;
					cy[ccnt] = psy + c0len[minj];
					break;
				case 2:	// 下
					cx[ccnt] = pex - c0len[minj];
					cy[ccnt] = pey;
					break;
				case 3:	// 左
					cx[ccnt] = psx;
					cy[ccnt] = pey - c0len[minj];
					break;
					}
					ccnt++;
					break;
			}
			//c0cidx[minj] = ccnt_;
		}
	}

	// データコピー
	int cvcnt1_=0;
	for( int i=0; i<ccnt; i++ ){
		int minj=-1;
		for( int j=0; j<c0cnt; j++ ){
			if( c0sort[j] == i ){
				minj = j;
				break;
			}
		}
		if( minj<0 ){
			continue;
		}
		cvtype[cvcnt1_] = c0type[minj];
		cvidx[cvcnt1_] = c0idx[minj];
		cvlen[cvcnt1_] = c0len[minj];
		// 座標変換
		cvx1[cvcnt1_] = c0x[minj];
		cvy1[cvcnt1_] = c0y[minj];
		cvcnt1_++;
	}
	if( cvcnt1_ != ccnt ){
		printf("error\n");
	}
	cvcnt1=cvcnt1_;
	ctype1=ctype0[k];
	//ccnt1++;

#if 1	// 曲線の向きを紙に対して右回りに
	int i;
	for( i=0; i<cvcnt1-1; i++ ){
		if( cvtype[i]==CVTYPE_FCURVE_START ){
			double vx0,vy0,vx1,vy1;
			vx0 = Xx2d[Xsidx] - Xx2d[Xsidx+1];
			vy0 = Xy2d[Xsidx] - Xy2d[Xsidx+1];
			vx1 = cvx1[i+1] - cvx1[i];
			vy1 = cvy1[i+1] - cvy1[i];
			if( vx0*vy1-vx1*vy0 < 0 ){ // 左側
				inverse_c(k);
			}
			break;
		}
		if( cvtype[k][i]==CVTYPE_FCURVE_END ){
			double vx0,vy0,vx1,vy1;
			vx0 = Xx2d[Xeidx] - Xx2d[Xeidx-1];
			vy0 = Xy2d[Xeidx] - Xy2d[Xeidx-1];
			vx1 = cvx1[i+1] - cvx1[i];
			vy1 = cvy1[i+1] - cvy1[i];
			if( vx0*vy1-vx1*vy0 < 0 ){
				inverse_c();
			}
			break;
		}
	}
	if( i == cvcnt1[k]-1 ){ // breakしていない＝折り線と交わらない
		for( i=0; i<cvcnt1-1; i++ ){
			if( cvtype[i]==CVTYPE_RUL_LEFT || cvtype[i]==CVTYPE_RUL_RIGHT ){
				double vx0,vy0,vx1,vy1;
				vx0 = cvx1[i] - Xx2d[cvidx[i]];
				vy0 = cvy1[i] - Xy2d[cvidx[i]];
				vx1 = cvx1[i+1] - cvx1[i];
				vy1 = cvy1[i+1] - cvy1[i];
				if( vx0*vy1-vx1*vy0 < 0 ){
					inverse_c();
				}
				break;
			}
		}
	}
#endif

#if 0
	{
		FILE *fp = fopen( "output/curve.csv", "w" );
		if( fp ){
			for( int i=0; i<ccnt1; i++ ){
				for( int j=0; j<cvcnt1[i]; j++ ){
					fprintf( fp, "%f,%f\n", cvx1[i][j], cvy1[i][j] );
				}
				fprintf( fp, "\n" );
			}
			for( int i=Xsidx; i<Xeidx+1; i++ ){
				fprintf( fp, "%f,%f\n", Xx2d[i], Xy2d[i] );
			}
			fclose(fp); fp=NULL;
		}
	}
#endif

end:
	return ret;
}
#endif


// ---------------------------------------------------------------
// INPUT:
// 曲線情報（ruling交差位置，紙に対して右回りに）
// ctype1 = -1:undefined, 0:trim, 1: fold
// cvtype=0:rl, 1:rr, 2:曲線始点, 3:曲線終点, 4:紙端
// cvidx= 0:rl, 1:rr -> index, 2:曲線始点, 3:曲線終点, 4:紙端
//int ccnt1, ctype1[MAXCN], cvcnt1[MAXCN], cvtype[MAXCN][MAXDCX0], cvidx[MAXCN][MAX_SPCNT];
//double cvlen[MAXCN][MAX_SPCNT], cvx1[MAXCN][MAX_SPCNT], cvy1[MAXCN][MAX_SPCNT];
//
// OUTPUT:
//int crcnt;
//crease crs[MAX_CRS_CNT];
// ---------------------------------------------------------------
int papermodel::FC2Crs()
{
	int rl=0, ret = 0;	

	//
	// 左右判定, 順番 -> int cli[MAXCN], clicnt=0, cri[MAXCN], cricnt=0;
	//
	int cli[MAXCN], clicnt=0, cri[MAXCN], cricnt=0;
	int idx[MAXCN]={-1}, icnt=0, flg[MAXCN]={1}, id0=-1;
	double avelen[MAXCN]={0.0};
	for( int i=0; i<fccnt; i++ ){
		flg[i]=1;

		int avecnt=0;
		for( int j=0; j<fcurve[i].cvcnt; j++ ){
			if( fcurve[i].cvtype[j]==CVTYPE_RUL_LEFT ){
				avelen[i] += fcurve[i].cvlen[j]; avecnt++;
			} else if( fcurve[i].cvtype[j]==CVTYPE_RUL_RIGHT ){
				avelen[i] -= fcurve[i].cvlen[j]; avecnt++;
			}
		}
		avelen[i] /= (double)avecnt;
	}
	for( int j=0; j<fccnt; j++ ){
		int mini=-1;
		double minlen=pw+ph;
		for( int i=0; i<fccnt; i++ ){
			if( flg[i] && minlen > avelen[i] ){
				mini = i;
				minlen = avelen[i];
			}
		}
		if( mini>-1 ){
			idx[icnt] = mini; icnt++;
			flg[mini] = 0;
		} else {
			break;
		}
	}
	for( id0=0; id0<icnt; id0++ ){
		if( avelen[idx[id0]] > 0 ){
			break;
		}
	}
	for( int j=id0; j<icnt; j++ ){
		cli[clicnt]=idx[j]; clicnt++;
	}
	for( int j=id0-1; j>=0; j-- ){
		cri[cricnt]=idx[j]; cricnt++;
	}
#if 0
	printf(" cli = ");
	for( int i=0; i<clicnt; i++ ){
		printf("%d, ", cli[i]);
	}
	printf("\n cri = ");
	for( int i=0; i<cricnt; i++ ){
		printf("%d, ", cri[i]);
	}
	printf("\n");
#endif

	//
	// 並べ替え
	//
	int crli[MAXCN], crri[MAXCN];
	for( int i=0; i<clicnt; i++ ){
		crli[i]=0;
		for( int j=1; j<crcnt; j++ ){
			if( cli[i] == crs[j].org_idx ){
				crli[i] = j;
				break;
			}
		}
	}
	for( int i=0; i<cricnt; i++ ){
		crri[i]=0;
		for( int j=1; j<crcnt; j++ ){
			if( cri[i] == crs[j].org_idx ){
				crri[i] = j;
				break;
			}
		}
	}

	for( int i=0; i<clicnt; i++ ){
		if( crli[i]>0 ){
			lcrs[i+1] =  &(crs[ crli[i] ]);
		} else {
			lcrs[i+1] =  &(crs[ crcnt ]); crcnt++; lcrcnt++;
		}
		curveintersect *fc = &(fcurve[cli[i]]);
		lcrs[i+1]->org_idx = cli[i];
		//lcrs[i+1]->org_cnt = fc->cvcnt;
		//lcrs[i+1]->org_x = fc->cvx;
		//lcrs[i+1]->org_y = fc->cvy;
		lcrs[i+1]->org_cnt = fc->dc->cvcnt;
		lcrs[i+1]->org_x = fc->dc->cvx;
		lcrs[i+1]->org_y = fc->dc->cvy;
		lcrs[i+1]->rl = -1; // left

		int jj0=-1, jj1=-1;
		for( int j=0; j<fc->cvcnt; j++){
			if( fc->cvtype[j]==CVTYPE_PAPER_TOP || fc->cvtype[j]==CVTYPE_PAPER_BOTTOM
				|| fc->cvtype[j]==CVTYPE_PAPER_LEFT || fc->cvtype[j]==CVTYPE_PAPER_RIGHT )
			{
				jj0=j;
				lcrs[i+1]->Xxs2d0 = fc->cvx[j];
				lcrs[i+1]->Xys2d0 = fc->cvy[j];
				break; // j
			}
		}
		for( int j=fc->cvcnt-1; j>jj0; j--){
			if( fc->cvtype[j]==CVTYPE_PAPER_TOP || fc->cvtype[j]==CVTYPE_PAPER_BOTTOM
				|| fc->cvtype[j]==CVTYPE_PAPER_LEFT || fc->cvtype[j]==CVTYPE_PAPER_RIGHT )
			{
				jj1=j;
				lcrs[i+1]->Xxe2d0 = fc->cvx[j];
				lcrs[i+1]->Xye2d0 = fc->cvy[j];
				break; // j
			}
		}

		fc->ctype=CTYPE_FOLD_LEFT;
		fc->cridx=i+1; // 対応する lcrs[], rcrs[] の index
	}
	for( int i=0; i<cricnt; i++ ){
		if( crri[i]>0 ){
			rcrs[i+1] =  &(crs[ crri[i] ]);
		} else {
			rcrs[i+1] =  &(crs[ crcnt ]); crcnt++; rcrcnt++;
		}
		curveintersect *fc = &(fcurve[cri[i]]);
		rcrs[i+1]->org_idx = cri[i];
		//rcrs[i+1]->org_cnt = fc->cvcnt;
		//rcrs[i+1]->org_x = fc->cvx;
		//rcrs[i+1]->org_y = fc->cvy;
		rcrs[i+1]->org_cnt = fc->dc->cvcnt;
		rcrs[i+1]->org_x = fc->dc->cvx;
		rcrs[i+1]->org_y = fc->dc->cvy;
		rcrs[i+1]->rl = 1; // right

		int jj0=-1, jj1=-1;
		for( int j=0; j<fc->cvcnt; j++){
			if( fc->cvtype[j]==CVTYPE_PAPER_TOP || fc->cvtype[j]==CVTYPE_PAPER_BOTTOM
				|| fc->cvtype[j]==CVTYPE_PAPER_LEFT || fc->cvtype[j]==CVTYPE_PAPER_RIGHT)
			{
				jj0=j;
				rcrs[i+1]->Xxs2d0 = fc->cvx[j];
				rcrs[i+1]->Xys2d0 = fc->cvy[j];
				break; // j
			}
		}
		for( int j=fc->cvcnt-1; j>jj0; j--){
			if( fc->cvtype[j]==CVTYPE_PAPER_TOP || fc->cvtype[j]==CVTYPE_PAPER_BOTTOM
				|| fc->cvtype[j]==CVTYPE_PAPER_LEFT || fc->cvtype[j]==CVTYPE_PAPER_RIGHT)
			{
				jj1=j;
				rcrs[i+1]->Xxe2d0 = fc->cvx[j];
				rcrs[i+1]->Xye2d0 = fc->cvy[j];
				break; // j
			}
		}
		fc->ctype=CTYPE_FOLD_RIGHT;
		fc->cridx=i+1; // 対応する lcrs[], rcrs[] の index
	}

	// ↑でインクリメント済
	//crcnt = clicnt+cricnt+1;
	//lcrcnt = clicnt+1;
	//rcrcnt = cricnt+1;
	//for( int j=0; j<crcnt; j++ ){
	//	crs[j].init_rlflg();
	//	crs[j].init_rrflg();
	//}

end:
	return ret;
}


int papermodel::addCrease( int lcrs_idx, int rcrs_idx, double slen, double dlen )
{
	// int crcnt, lcrcnt, rcrcnt;
	// crease crs[MAX_CRS_CNT];
	// crease *lcrs[MAX_CRS_CNT], *rcrs[MAX_CRS_CNT];
	// が求まっていること

	int ret=0;
	if( lcrs_idx <= 0 ){

	} else {
		for( int j=lcrs_idx; j<lcrcnt; j++ ){
			crease *c0 = lcrs[j-1];
			crease *c1 = lcrs[j];

			if( slen>0.0 ){
				double len=slen;
				c0->init_rlflg();
				c1->initX();
				for( int i=c0->Xsidx; i<=c0->Xeidx; i++, len+=dlen ){
					if( len<0 ){ len=0.0; }
					c0->rllen[i] = len;
					c0->rlflg[i] = RETYPE_CURVE;
				}
			} else if( c1->flg_org ){
				c0->init_rlflg();
				c1->initX();
				ret = setRulLen_SplineCP( c0->Xx2d, c0->Xy2d, c0->rlx_cp, c0->rly_cp, c0->Xsidx, c0->Xeidx,
					c1->CPx, c1->CPy, c0->rllen, c0->rlflg );
			} else { // if( c1->flg_org )
				// cx1[],cy1[] 平滑化 -> c0->rllen
				c0->init_rlflg();
				c1->initX();
				ret = setRulLen_curve( &(c0->Xx2d[c0->Xsidx]), &(c0->Xy2d[c0->Xsidx]),
					&(c0->rlx_cp[c0->Xsidx]), &(c0->rly_cp[c0->Xsidx]),
					c0->Xeidx-c0->Xsidx+1,
					c1->org_x, c1->org_y, c1->org_cnt, c1->CPx, c1->CPy,
					&(c0->rllen[c0->Xsidx]), &(c0->rlflg[c0->Xsidx]) );
				c1->flg_org = 1;
			}
			c1->setL2R( c0 );
			c1->calcTN2d();
			c1->calcTNB( c1->Xsidx, c1->Xeidx );
			c1->calcLeft();
		} // j
	}

	if( rcrs_idx <= 0 ){

	} else {
		for( int j=rcrs_idx; j<rcrcnt; j++ ){
			crease *c0 = rcrs[j-1];
			crease *c1 = rcrs[j];

			if( slen>0.0 ){
				double len=slen;
				c0->init_rrflg();
				c1->initX();
				for( int i=c0->Xsidx; i<=c0->Xeidx; i++, len+=dlen ){
					if( len<0 ){ len=0.0; }
					c0->rrlen[i] = len;
					c0->rrflg[i] = RETYPE_CURVE;
				}
			} else if( c1->flg_org ){
				c0->init_rrflg();
				c1->initX();
				ret = setRulLen_SplineCP( c0->Xx2d, c0->Xy2d, c0->rrx_cp, c0->rry_cp, c0->Xsidx, c0->Xeidx,
					c1->CPx, c1->CPy, c0->rrlen, c0->rrflg );
			} else {
				// cx1[],cy1[] 平滑化 -> c0->rllen
				c0->init_rrflg();
				c1->initX();
				ret = setRulLen_curve( &(c0->Xx2d[c0->Xsidx]), &(c0->Xy2d[c0->Xsidx]),
					&(c0->rrx_cp[c0->Xsidx]), &(c0->rry_cp[c0->Xsidx]),
					c0->Xeidx-c0->Xsidx+1,
					c1->org_x, c1->org_y, c1->org_cnt, c1->CPx, c1->CPy, 
					&(c0->rrlen[c0->Xsidx]), &(c0->rrflg[c0->Xsidx]) );
				c1->flg_org = 1;
			}

			c1->setR2L( c0 );
			c1->calcTN2d();
			c1->calcTNB( c1->Xsidx, c1->Xeidx );
			c1->calcRight();
		} // j
	}
end:
	return ret;
}
