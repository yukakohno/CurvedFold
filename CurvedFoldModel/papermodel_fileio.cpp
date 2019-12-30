#include <stdio.h>
#include <memory.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "papermodel.h"
#include "util.h"

int papermodel::dumpsvg00( char *fname, int divnum )
{
	int i, j, li, ret=0, idxk[MAX_CPCNT], idxt[MAX_CPCNT], idxa[MAX_CPCNT];
	double sc=3.55;

	FILE *fp=fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"%f %f %f %f\" >\n",
			(double)-pex*sc*1.3, (double)-pey*sc*1.5, (double)pex*2*sc*1.3, (double)pey*2*sc*1.5);
		// #000: boundary
		// #f00: yamaori

#if 1	// checker pattern
		for( i=-20; i<20; i++ ){
			for( j=-22; j<22; j++ ){
				if( (i+j)%2==0 ){
					continue;
				}
				fprintf( fp, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(214,214,214)\" />\n",
					(double)i*200.0/16.0*sc, (double)j*200.0/16.0*sc, 200.0/16.0*sc, 200.0/16.0*sc );
			}
		}
#endif
		for( int d=0; d<divnum; d++ ){
			double dangle = (double)d*2.0*M_PI/(double)divnum -15.0/180.0*M_PI;
			double cosd = cos(dangle);
			double sind = sin(dangle);
			double x0, y0, x1, y1;
#if 0		// checker pattern
			for( i=0; i<16; i++ ){
				for( j=0; j<16; j++ ){
					if( (i+j)%2 ){
						continue;
					}
					x0 = (double)i*200.0/16.0*sc;
					y0 = (double)j*200.0/16.0*sc;
					x1 = 200.0/16.0*sc;
					y1 = 200.0/16.0*sc;
					fprintf( fp, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" style=\"fill:rgb(214,214,214)\" />\n",
						//x0, y0, x1, y1 );
						cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );
				}
			}
#endif
			// crease curve
			for( int k=0; k<2; k++ )
			{
				crease *c = &(crs[k]);
				x0 = -c->Xxs2d*sc;
				y0 = c->Xys2d*sc;
				x1 = -c->Xx2d[c->Xsidx]*sc;
				y1 = c->Xy2d[c->Xsidx]*sc;
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					//x0, y0, x1, y1 );
					cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );
				for( i=c->Xsidx+1, li=c->Xsidx; i<=c->Xeidx; i++ ){
					x0 = -c->Xx2d[li]*sc;
					y0 = c->Xy2d[li]*sc;
					x1 = -c->Xx2d[ i]*sc;
					y1 = c->Xy2d[ i]*sc;
					fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
						//x0, y0, x1, y1 );
						cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );
					li=i;
				}
				x0 = -c->Xx2d[c->Xeidx]*sc;
				y0 = c->Xy2d[c->Xeidx]*sc;
				x1 = -c->Xxe2d*sc;
				y1 = c->Xye2d*sc;
				fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
					//x0, y0, x1, y1 );
					cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );
			}
			x0 = -crs[0].Xxe2d*sc;
			y0 = crs[0].Xye2d*sc;
			x1 = -crs[1].Xxe2d*sc;
			y1 = crs[1].Xye2d*sc;
			fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				//x0, y0, x1, y1 );
				cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );

			x0 = -crs[0].Xxe2d*sc;
			y0 = crs[0].Xye2d*sc;
			x1 = -pex*sc;
			y1 = pey*sc;
			fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				//x0, y0, x1, y1 );
				cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );

			x0 = -pex*sc;
			y0 = pey*sc;
			x1 = -crs[2].Xxe2d*sc;
			y1 = crs[2].Xye2d*sc;
			fprintf( fp, "<line stroke=\"#000\" opacity=\"1\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"2\"/>\n",
				//x0, y0, x1, y1 );
				cosd*x0+sind*y0, -sind*x0+cosd*y0, cosd*x1+sind*y1, -sind*x1+cosd*y1 );
		}
#if 0	// trim
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
#endif
#if 0		// àÛç¸Ç≈í[Ç™êÿÇÍÇÈÇÃÇñhÇÆ
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

#if 0		// àÛç¸Ç≈í[Ç™êÿÇÍÇÈÇÃÇñhÇÆ
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

#if 0	// àÛç¸Ç≈í[Ç™êÿÇÍÇÈÇÃÇñhÇÆ
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

// ê‹ÇËê¸ïîï™ÇÃ3Då`èÛ
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

// í∏ì_é¸ÇËÇÃäpìx
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

// quad ÇÃïΩñ ìx
int papermodel::checkquadplane( char *fname )
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

int papermodel::checkGap( char *fname, int divnum )
{
	int ret=0, spcnt, p3cnt0, p3cnt1;
	double space = 10;
	double stx0[MAX_STP_CNT], sty0[MAX_STP_CNT], stz0[MAX_STP_CNT];
	double stx1[MAX_STP_CNT], sty1[MAX_STP_CNT], stz1[MAX_STP_CNT];
	double dif[MAX_STP_CNT];

	if( lcrcnt<=1 && rcrcnt<=1 ){
		return -1;
	}
#if 1
	double errdata[ MAX_STP_CNT*7 ];	memset( errdata, 0, sizeof(double)*MAX_STP_CNT*7 );
	spcnt = checkGap( errdata, MAX_STP_CNT, 7, divnum );

	FILE *fp = fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "i, dif, stx0, sty0, stz0, stx1, sty1, stz1\n" );
		for( int i=0; i<spcnt; i++ ){
			double *ed = &(errdata[i*7]);
			fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f\n",
				i, ed[0], ed[1], ed[2], ed[3], ed[4], ed[5], ed[6] );
		}
		fclose( fp );
	}
#else
	p3cnt0 = lcrs[1]->getSamplePoints3D( space, stx0, sty0, stz0 );
	p3cnt1 = rcrs[1]->getSamplePoints3D( space, stx1, sty1, stz1 );
	spcnt = p3cnt0 < p3cnt1 ? p3cnt0 : p3cnt1;

	// rotation & mid points
	double axis[3]={0.0,0.0,1.0}, quat[4], rmat[16];
	axis_ang_quat( axis, 2.0*M_PI/(double)divnum, quat );
	quat_mat( quat, rmat );

	for( int i=0; i<spcnt; i++ )
	{
		double tmpx, tmpy, tmpz;
		tmpx = rmat[0]*stx1[i] + rmat[4]*sty1[i] + rmat[8]*stz1[i] + rmat[12];
		tmpy = rmat[1]*stx1[i] + rmat[5]*sty1[i] + rmat[9]*stz1[i] + rmat[13];
		tmpz = rmat[2]*stx1[i] + rmat[6]*sty1[i] + rmat[10]*stz1[i] + rmat[14];
		stx1[i] = tmpx;
		sty1[i] = tmpy;
		stz1[i] = tmpz;
		tmpx = stx1[i]-stx0[i];
		tmpy = sty1[i]-sty0[i];
		tmpz = stz1[i]-stz0[i];
		dif[i] = sqrt( tmpx*tmpx + tmpy*tmpy + tmpz*tmpz );
	}

	FILE *fp = fopen( fname, "w" );
	if( fp ){
		fprintf( fp, "i, dif, stx0, sty0, stz0, stx1, sty1, stz1\n" );
		for( int i=0; i<spcnt; i++ ){
			fprintf( fp, "%d,%f,%f,%f,%f,%f,%f,%f\n",
				i, dif[i], stx0[i], sty0[i], stz0[i], stx1[i], sty1[i], stz1[i] );
		}
		fclose( fp );
	}
#endif
	return ret;
}

