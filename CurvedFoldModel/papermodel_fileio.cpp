#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "papermodel.h"
#include "util.h"

int papermodel::dumpsvg0( char *fname )
{
	int i, j, li, ret=0, idxk[MAX_CPCNT], idxt[MAX_CPCNT], idxa[MAX_CPCNT];
	double sc=3.55;

	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"%f %f %f %f\" >\n",
			(double)psx*sc-18, (double)psy*sc-18, (double)(pex-psx)*sc+36, (double)(pey-psy)*sc+36);
		// #000: boundary
		// #f00: yamaori

		// checker pattern
		for( i=0; i<16; i++ ){
			for( j=0; j<16; j++ ){
				if( (i+j)%2 ){
					continue;
				}
				fprintf( fp, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(214,214,214)\" />\n",
					(double)i*200.0/16.0*sc, (double)j*200.0/16.0*sc, 200.0/16.0*sc, 200.0/16.0*sc );
			}
		}

		// crease curve
		for( int k=0; k<crcnt; k++ )
		{
			crease *c = &(crs[k]);
			fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				c->Xxs2d*sc, c->Xys2d*sc, c->Xx2d[c->Xsidx]*sc, c->Xy2d[c->Xsidx]*sc );
			for( i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					c->Xx2d[li]*sc, c->Xy2d[li]*sc, c->Xx2d[ i]*sc, c->Xy2d[ i]*sc );
				li=i;
			}
			fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				c->Xx2d[c->Xeidx]*sc, c->Xy2d[c->Xeidx]*sc, c->Xxe2d*sc, c->Xye2d*sc );
		}

		// trim
		for( int j=0; j<dccnt; j++ ){
			if( dcurve[j].ctype != CTYPE_TRIM ){
				continue;
			}
			for( int i=1; i<dcurve[j].cvcnt; i++ ){
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					dcurve[j].cvx[j]*sc, dcurve[j].cvy[i-1]*sc, dcurve[j].cvx[i]*sc, dcurve[j].cvy[i]*sc );
			}
		}

		// boundary
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc, 0*sc, 200*sc, 0*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc, 0*sc, 200*sc, 200*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc, 200*sc, 0*sc, 200*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc, 200*sc, 0*sc, 0*sc );

#if 0		// 印刷で端が切れるのを防ぐ
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc-18, 0*sc-18, 200*sc+18, 0*sc-18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc+18, 0*sc-18, 200*sc+18, 200*sc+18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc+18, 200*sc+18, 0*sc-18, 200*sc+18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc-18, 200*sc+18, 0*sc-18, 0*sc-18 );
#endif
		fprintf( fp, "</svg>\n" );
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

int papermodel::dumpsvg1( char *fname )
{
	int ret=0;// idxk[MAX_CPCNT], idxt[MAX_CPCNT], idxa[MAX_CPCNT];
	double sc=3.55;
	char strcol[5];

	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"%f %f %f %f\" >\n",
			(double)psx*sc-18, (double)psy*sc-18, (double)(pex-psx)*sc+36, (double)(pey-psy)*sc+36);

		// #000: boundary
		// #f00: yamaori
		// #00f: taniori
		// #ff0: papermodel

		// primary crease
		crease *c = &(crs[0]);
		fprintf( fp, "<line stroke=\"#f00\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			c->Xxs2d*sc, c->Xys2d*sc, c->Xx2d[c->Xsidx]*sc, c->Xy2d[c->Xsidx]*sc );
		for( int i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
			fprintf( fp, "<line stroke=\"#f00\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				c->Xx2d[li]*sc, c->Xy2d[li]*sc, c->Xx2d[ i]*sc, c->Xy2d[ i]*sc );
			li=i;
		}
		fprintf( fp, "<line stroke=\"#f00\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			c->Xx2d[c->Xeidx]*sc, c->Xy2d[c->Xeidx]*sc, c->Xxe2d*sc, c->Xye2d*sc );

		// crease left
		for( int j=1; j<lcrcnt; j++ ){
			crease *c = lcrs[j];
			if( j%2 ){
				sprintf( strcol, "#00f" );
			} else {
				sprintf( strcol, "#f00" );
			}
			fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				strcol, c->Xxs2d*sc, c->Xys2d*sc, c->Xx2d[c->Xsidx]*sc, c->Xy2d[c->Xsidx]*sc );
			for( int i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
				fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					strcol, c->Xx2d[li]*sc, c->Xy2d[li]*sc, c->Xx2d[ i]*sc, c->Xy2d[ i]*sc );
				li=i;
			}
			fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				strcol, c->Xx2d[c->Xeidx]*sc, c->Xy2d[c->Xeidx]*sc, c->Xxe2d*sc, c->Xye2d*sc );
		} // j

		// crease right
		for( int j=1; j<rcrcnt; j++ ){
			crease *c = rcrs[j];
			if( j%2 ){
				sprintf( strcol, "#00f" );
			} else {
				sprintf( strcol, "#f00" );
			}
			fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				strcol, c->Xxs2d*sc, c->Xys2d*sc, c->Xx2d[c->Xsidx]*sc, c->Xy2d[c->Xsidx]*sc );
			for( int i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
				fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					strcol, c->Xx2d[li]*sc, c->Xy2d[li]*sc, c->Xx2d[ i]*sc, c->Xy2d[ i]*sc );
				li=i;
			}
			fprintf( fp, "<line stroke=\"%s\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				strcol, c->Xx2d[c->Xeidx]*sc, c->Xy2d[c->Xeidx]*sc, c->Xxe2d*sc, c->Xye2d*sc );
		} // j

		// rulings left
		for( int j=0; j<lcrcnt; j++ ){
			crease *c = lcrs[j];
			for( int i=c->Xsidx; i<=c->Xeidx; i++ ){
				if( c->rlflg[i]==RETYPE_UNDEF || c->rllen[i]==0.0 ){
					continue;
				}
				fprintf( fp, "<line stroke=\"#ff0\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					c->Xx2d[i]*sc, c->Xy2d[i]*sc,
					(c->Xx2d[i] + c->rlx_cp[i] * c->rllen[i]) *sc,
					(c->Xy2d[i] + c->rly_cp[i] * c->rllen[i]) *sc );
			}
		}

		// rulings right
		for( int j=0; j<rcrcnt; j++ ){
			crease *c = rcrs[j];
			for( int i=c->Xsidx; i<=c->Xeidx; i++ ){
				if( c->rrflg[i]==RETYPE_UNDEF || c->rrlen[i]==0.0 ){
					continue;
				}
				fprintf( fp, "<line stroke=\"#ff0\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					c->Xx2d[i]*sc, c->Xy2d[i]*sc,
					(c->Xx2d[i] + c->rrx_cp[i] * c->rrlen[i]) *sc,
					(c->Xy2d[i] + c->rry_cp[i] * c->rrlen[i]) *sc );
			}
		}

		{	// edge left
			crease *c = lcrs[lcrcnt-1];
			for( int i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
				if( c->rlflg[i]==RETYPE_UNDEF || c->rllen[i]==0.0 || c->rlflg[li]==RETYPE_UNDEF || c->rllen[li]==0.0 ){
					li=i;
					continue;
				}
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					(c->Xx2d[li] + c->rlx_cp[li] * c->rllen[li]) *sc,
					(c->Xy2d[li] + c->rly_cp[li] * c->rllen[li]) *sc,
					(c->Xx2d[ i] + c->rlx_cp[ i] * c->rllen[ i]) *sc,
					(c->Xy2d[ i] + c->rly_cp[ i] * c->rllen[ i]) *sc );
				li=i;
			}
		}
		{	// edge right
			crease *c = rcrs[rcrcnt-1];
			for( int i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
				if( c->rrflg[i]==RETYPE_UNDEF || c->rrlen[i]==0.0 || c->rrflg[li]==RETYPE_UNDEF || c->rrlen[li]==0.0 ){
					li=i;
					continue;
				}
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					(c->Xx2d[li] + c->rrx_cp[li] * c->rrlen[li]) *sc,
					(c->Xy2d[li] + c->rry_cp[li] * c->rrlen[li]) *sc,
					(c->Xx2d[ i] + c->rrx_cp[ i] * c->rrlen[ i]) *sc,
					(c->Xy2d[ i] + c->rry_cp[ i] * c->rrlen[ i]) *sc );
				li=i;
			}
		}

		// trim
		for( int j=0; j<dccnt; j++ ){
			if( dcurve[j].ctype != CTYPE_TRIM ){
				continue;
			}
			for( int i=1; i<dcurve[j].cvcnt; i++ ){
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					dcurve[j].cvx[j]*sc, dcurve[j].cvy[i-1]*sc, dcurve[j].cvx[i]*sc, dcurve[j].cvy[i]*sc );
			}
		}

		// boundary
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc, 0*sc, 200*sc, 0*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc, 0*sc, 200*sc, 200*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc, 200*sc, 0*sc, 200*sc );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc, 200*sc, 0*sc, 0*sc );

#if 0	// 印刷で端が切れるのを防ぐ
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc-18, 0*sc-18, 200*sc+18, 0*sc-18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc+18, 0*sc-18, 200*sc+18, 200*sc+18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			200*sc+18, 200*sc+18, 0*sc-18, 200*sc+18 );
		fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
			0*sc-18, 200*sc+18, 0*sc-18, 0*sc-18 );
#endif
		fprintf( fp, "</svg>\n" );
		fclose(fp); fp=NULL;
	}

end:
	return ret;
}

// 折り線部分の3D形状
// Obj
//	int o_vcnt, o_fcnt;
//	double o_vx[MAX_SPCNT*(MAXCN+2)*2], o_vy[MAX_SPCNT*(MAXCN+2)*2], o_vz[MAX_SPCNT*(MAXCN+2)*2];
//	int o_fi[MAX_SPCNT*(MAXCN+1)][4];
int papermodel::dumpObj( char *fname, double scl=1.0 )
{
	int ret=0;
	double minx=100000, maxx=-100000, miny=100000, maxy=-100000, minz=100000, maxz=-100000;

	FILE *fp = fopen( fname, "w" );
	if( fp==NULL ){
		ret = -1;
		goto end;
	}
	for( int i=0; i<o_vcnt; i++ ){
		double xx=o_vx[i]*scl, yy=o_vy[i]*scl, zz=o_vz[i]*scl;
		fprintf( fp, "v %f %f %f\n", xx, yy, zz);
		minx = minx < xx ? minx : xx;
		maxx = maxx > xx ? maxx : xx;
		miny = miny < yy ? miny : yy;
		maxy = maxy > yy ? maxy : yy;
		minz = minz < zz ? minz : zz;
		maxz = maxz > zz ? maxz : zz;
	}
	for( int i=0; i<o_fcnt; i++ ){
		switch( o_fvcnt[i] ){
			case 3:	fprintf( fp, "f %d %d %d\n", o_fi[i][0]+1, o_fi[i][1]+1, o_fi[i][2]+1 );	break;
			case 4:	fprintf( fp, "f %d %d %d %d\n", o_fi[i][0]+1, o_fi[i][1]+1, o_fi[i][2]+1, o_fi[i][3]+1 );	break;
		}
	}

	printf("x = [ %f, %f ] %f, y = [ %f, %f ] %f, z = [ %f, %f ] %f\n",
		minx,maxx,maxx-minx, miny,maxy,maxy-miny, minz,maxz,maxz-minz );
end:
	if(fp) fclose(fp); fp=NULL;
	return ret;
}

// 頂点周りの角度
int papermodel::check180( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");
	if( !fp ){
		ret = -1; goto end;
	}

	ret = crs[0].check180( fp );
	for( int j=1; j<lcrcnt; j++ )
	{
		fprintf( fp, "lcrs[%d]\n", j );
		ret = lcrs[j]->check180( fp );
	}
	for( int j=1; j<rcrcnt; j++ )
	{
		fprintf( fp, "rcrs[%d]\n", j );
		ret = rcrs[j]->check180( fp );
	}
end:
	if( fp ){ fclose(fp); }
	return ret;
}

// quad の平面度
int papermodel::checkquatplane( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");
	if( !fp ){
		ret = -1; goto end;
	}
	ret = crs[0].checkquadplane( fp );
	for( int j=1; j<lcrcnt; j++ )
	{
		fprintf( fp, "lcrs[%d]\n", j );
		ret = lcrs[j]->checkquadplane( fp );
	}
	for( int j=1; j<rcrcnt; j++ )
	{
		fprintf( fp, "rcrs[%d]\n", j );
		ret = rcrs[j]->checkquadplane( fp );
	}
end:
	if( fp ){ fclose(fp); }
	return ret;
}

int papermodel::dumpcv0( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");
	if( fp ){
		fprintf(fp, "%d	// dccnt\n", dccnt);
		for( int i=0; i<dccnt; i++ ){
			fprintf(fp, "%d %d // ctype[%d], cvcnt[%d]\n", dcurve[i].ctype, dcurve[i].cvcnt, i, i );
			for( int j=0; j<dcurve[i].cvcnt; j++ ){
				fprintf(fp, "%f %f\n", dcurve[i].cvx[j], dcurve[i].cvy[j]);
			}
		}
		fclose(fp); fp=NULL;
	}
	return ret;
}

int papermodel::loadcv0( char *fname )
{
	int ret=0, tmp_ctype=0, tmp_cvcnt=0;
	double tmp_cvx, tmp_cvy;
	char buf[1024];
	FILE *fp=fopen(fname, "r");
	if( fp ){
		fgets(buf, 1024, fp);	sscanf( buf, "%d", &dccnt );
		for( int i=0; i<dccnt; i++ ){
			fgets(buf, 1024, fp);	sscanf( buf, "%d %d", &tmp_ctype, &tmp_cvcnt );
			switch( tmp_ctype ){
				case 1:	dcurve[i].ctype = CTYPE_TRIM;	break;
				case 2:	dcurve[i].ctype = CTYPE_FOLD;	break;
				case 3:	dcurve[i].ctype = CTYPE_FOLD_LEFT;	break;
				case 4:	dcurve[i].ctype = CTYPE_FOLD_RIGHT;	break;
				default:	dcurve[i].ctype = CTYPE_UNDEF;	break;
			}
			dcurve[i].cvcnt = tmp_cvcnt;
			for( int j=0; j<dcurve[i].cvcnt; j++ ){
				fgets( buf, 1024, fp );
				sscanf( buf, "%lf %lf", &tmp_cvx, &tmp_cvy );
				dcurve[i].cvx[j] = tmp_cvx;
				dcurve[i].cvy[j] = tmp_cvy;
			}
		}
	} else {
		dccnt = 0;
	}
	return ret;
}

#if 0
int papermodel::dumpcv1( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");
	if( fp ){
		//int ccnt1, ctype1[MAXCN], cvcnt1[MAXCN], cvtype[MAXCN][MAX_SPCNT], cvidx[MAXCN][MAX_SPCNT];
		//double cvlen[MAXCN][MAX_SPCNT], cvx1[MAXCN][MAX_SPCNT], cvy1[MAXCN][MAX_SPCNT];
		//double cvrllen[MAXCN][MAX_SPCNT], cvrrlen[MAXCN][MAX_SPCNT];	// 曲線で切った後のruling長さ
		fprintf(fp, "%d	// ccnt1\n", ccnt1);
		for( int i=0; i<ccnt1; i++ ){
			fprintf(fp, "%d // ctype1[%d]\n", ctype1[i], i);
			fprintf(fp, "%d // cvcnt1[%d]\n", cvcnt1[i], i);
			for( int j=0; j<cvcnt1[i]; j++ ){
				fprintf(fp, "%f %f %f %d %d // cvx1, cvy1, cvlen, cvtype, cvidx\n",
					cvx1[i][j], cvy1[i][j], cvlen[i][j], cvtype[i][j], cvidx[i][j]);
			}
		}
		fclose(fp); fp=NULL;
	}
	return ret;
}

int papermodel::loadcv1( char *fname )
{
	int ret=0, tmp_cvtype=0;
	char buf[1024];
	FILE *fp=fopen(fname, "r");
	if( fp ){
		fgets(buf, 1024, fp);	sscanf( buf, "%d", &ccnt1 );
		for( int i=0; i<ccnt1; i++ ){
			fgets(buf, 1024, fp);	sscanf( buf, "%d", &(ctype1[i]) );
			fgets(buf, 1024, fp);	sscanf( buf, "%d", &(cvcnt1[i]) );
			for( int j=0; j<cvcnt1[i]; j++ ){
				fgets(buf, 1024, fp);
				sscanf(buf, "%lf %lf %lf %d %d",
					&(cvx1[i][j]), &(cvy1[i][j]), &(cvlen[i][j]), &(tmp_cvtype), &(cvidx[i][j]));
				// 交差対称 0:rl, 1:rr, 2:曲線始点, 3:曲線終点, 4:紙端
				switch( tmp_cvtype){
					case 1: cvtype[i][j] = CVTYPE_RUL_LEFT; break;
					case 2: cvtype[i][j] = CVTYPE_RUL_RIGHT; break;
					case 3: cvtype[i][j] = CVTYPE_FCURVE_START; break;
					case 4: cvtype[i][j] = CVTYPE_FCURVE_END; break;
					case 5: cvtype[i][j] = CVTYPE_PAPER_EDGE; break;
					default: cvtype[i][j] = CVTYPE_UNDEF; break;
				}
			}
		}
		fclose(fp); fp=NULL;
	}
	return ret;
}
#endif

int papermodel::dumpBsplineCP( char *fname )
{
	int ret=0;
	FILE *fp=fopen(fname, "w");
	if( fp ){
		fprintf( fp, "%d %d // lcrcnt rcrcnt\n", lcrcnt-1, rcrcnt-1 );
		for( int i=1; i<lcrcnt; i++ ){
			crease *c = lcrs[i];
			fprintf(fp, "%d %d %d // c0/CP rl org_c0_idx\n", c->flg_org, c->rl, c->org_idx);
			for( int j=0; j<CCNT; j++ ){
				fprintf(fp, "%f %f\n", c->CPx[j], c->CPy[j]);
			}
		}
		for( int i=1; i<rcrcnt; i++ ){
			crease *c = rcrs[i];
			fprintf(fp, "%d %d %d // c0/CP rl org_c0_idx\n", c->flg_org, c->rl, c->org_idx);
			for( int j=0; j<CCNT; j++ ){
				fprintf(fp, "%f %f\n", c->CPx[j], c->CPy[j]);
			}
		}
		fclose(fp); fp=NULL;
	}
	return ret;
}

int papermodel::loadBsplineCP( char *fname )
{
	int cnt=1, ret=0;
	char buf[1024];
	FILE *fp=fopen(fname, "r");
	if( fp ){
		//fprintf( fp, "%d %d// lcrcnt rcrcnt\n", lcrcnt-1, rcrcnt-1 );
		fgets(buf,1024,fp); sscanf( buf, "%d %d", &lcrcnt, &rcrcnt );
		crcnt = lcrcnt+rcrcnt+cnt; lcrcnt+=cnt; rcrcnt+=cnt;
		for( int i=1; i<lcrcnt; i++ ){
			crease *c = lcrs[i] = &(crs[cnt]);
			//fprintf(fp, "%d %d %d // c0/CP rl org_c0_idx\n", c->flg_org, c->rl, c->org_idx);
			fgets(buf,1024,fp); sscanf( buf, "%d %d %d", &(c->flg_org), &(c->rl), &(c->org_idx) );
			c->org_cnt = fcurve[c->org_idx].cvcnt;
			c->org_x = fcurve[c->org_idx].cvx;
			c->org_y = fcurve[c->org_idx].cvy;
			for( int j=0; j<CCNT; j++ ){
				//fprintf(fp, "%f %f\n", c->CPx[j], c->CPy[j]);
				fgets(buf,1024,fp); sscanf( buf, "%lf %lf", &(c->CPx[j]), &(c->CPy[j]) );
			}
			cnt++;
		}
		for( int i=1; i<rcrcnt; i++ ){
			crease *c = rcrs[i] = &(crs[cnt]);
			//fprintf(fp, "%d %d %d // c0/CP rl org_c0_idx\n", c->flg_org, c->rl, c->org_idx);
			fgets(buf,1024,fp); sscanf( buf, "%d %d %d", &(c->flg_org), &(c->rl), &(c->org_idx) );
			c->org_cnt = fcurve[c->org_idx].cvcnt;
			c->org_x = fcurve[c->org_idx].cvx;
			c->org_y = fcurve[c->org_idx].cvy;
			for( int j=0; j<CCNT; j++ ){
				//fprintf(fp, "%f %f\n", c->CPx[j], c->CPy[j]);
				fgets(buf,1024,fp); sscanf( buf, "%lf %lf", &(c->CPx[j]), &(c->CPy[j]) );
			}
			cnt++;
		}
		fclose(fp); fp=NULL;
	}
	return ret;
}

