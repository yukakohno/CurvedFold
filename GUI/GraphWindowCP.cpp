#include <FL/fl_draw.H>
#define _USE_MATH_DEFINES
#include <math.h>
#include <opencv2/opencv.hpp>

#include "GraphWindowCP.h"
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

GraphWindowCP::GraphWindowCP(int X, int Y, int W, int H, const char *L) : Fl_Double_Window(X, Y, W, H, L)
{
	box(FL_FLAT_BOX);
	imgbuf = new unsigned char[w() * h() * 3];
	jpg = NULL;
	init();
	draw();
}

GraphWindowCP::GraphWindowCP(int X, int Y, int W, int H) : Fl_Double_Window(X, Y, W, H)
{
	box(FL_FLAT_BOX);
	imgbuf = new unsigned char[w() * h() * 3];
	jpg = NULL;
	init();
	draw();
}

void GraphWindowCP::init()
{
	ppm = NULL;

	ox = w()*0.1;
	oy = h()*0.1;
	wscl = (double)TEXWIDTH / (double)PPWIDTH;
	hscl = (double)TEXHEIGHT / (double)PPHEIGHT;

	ofs_x = ofs_y = ofs_a = ofs_sina = 0.0;
	ofs_cosa = 1.0;
	ofs_psx = ofs_psy = ofs_pex = ofs_pey = 0;

	ccnt = 0;

	if( jpg ){ delete jpg; }
	jpg = new Fl_JPEG_Image( TEXFNAME );

	disp_modTGT = -2;
}

void GraphWindowCP::draw()
{
	fl_color(FL_WHITE);
	fl_rectf(0,0,w(),h());

	if( ppm==NULL ){
		return;
	}

	int pw=(int)(ppm->pw*wscl), ph=(int)(ppm->ph*hscl);
	int psx=(int)(ppm->psx*wscl), psy=(int)(ppm->psy*hscl), pex=(int)(ppm->pex*wscl), pey=(int)(ppm->pey*hscl);
	int cnt = ppm->crs[0].Xcnt;
	double len0 = 5;
	double *m2 = ppm->crs[0].m2;
	int plcnt=ppm->plcnt, *plvcnt=ppm->plvcnt;
	double *plx=ppm->plx_cp, *ply=ppm->ply_cp;

	double tmpx0,tmpy0,tmpx1,tmpy1;
	double mr[9],mt1[9],mtofs[9],xx=0.0,yy=0.0,ca=0.0,sa=1.0;

	fl_color(230, 230, 230);
	fl_rectf( (int)(ox+psx+ofs_psx), (int)(oy+psy+ofs_psy), (int)(pw+ofs_pex-ofs_psx), (int)(ph+ofs_pey-ofs_psy) );

	if( jpg && disp_PTN ){
		jpg->draw( (int)(ox+psx+ofs_psx), (int)(oy+psy+ofs_psy), (int)(pw+ofs_pex-ofs_psx), (int)(ph+ofs_pey-ofs_psy) );
	}

	if( m2 ){
		xx=m2[6]*wscl;
		yy=m2[7]*hscl;
		ca=m2[0];
		sa=m2[1];
	}

	unit_m33(mt);	mt[2]=ox+xx;	mt[5]=oy+yy;
	unit_m33(mr);	mr[0] = wscl;	mr[4] = hscl;
	unit_m33(mt1);	mt1[2] = -xx/wscl;	mt1[5] = -yy/hscl;
	mult_m33_n33( mt, mr, 1);	// dst=1 : n44[16] = n44[16] * m44[16], n44にm44を後ろからかける 
	mult_m33_n33( mt, mt1, 1);

	unit_m33(mtofs);	mtofs[2]=ox+xx+ofs_x;	mtofs[5]=oy+yy+ofs_y;
	unit_m33(mr);
	mr[0] = ofs_cosa*wscl;	mr[1] = -ofs_sina*wscl;
	mr[3] = ofs_sina*hscl;	mr[4] = ofs_cosa*hscl;
	mult_m33_n33( mtofs, mr, 1);	// dst=1 : n44[16] = n44[16] * m44[16], n44にm44を後ろからかける 
	mult_m33_n33( mtofs, mt1, 1);

	// addcurve
	fl_color(0, 0, 0);
	fl_line_style( FL_SOLID, 2 );
	for( int i=1; i<ccnt; i++ ){
		fl_line( cx[i-1], cy[i-1], cx[i], cy[i] );
	}

	for( int i=0; i<ppm->dccnt; i++ ){
		curvedraw *dc = &(ppm->dcurve[i]);
		switch( dc->ctype ){
			case CTYPE_TRIM: fl_color(255, 127, 127); break;
			case CTYPE_FOLD:
			case CTYPE_FOLD_LEFT:
			case CTYPE_FOLD_RIGHT:
				fl_color(127, 255, 127);
				break;
			default: fl_color(127, 127, 127); break;
		}
		fl_line_style( FL_SOLID, 2 );
		for( int j=0; j<dc->cvcnt-1; j++ ){
			fl_line((int)(mt[0]*dc->cvx[j] + mt[1]*dc->cvy[j] +mt[2]),
				(int)(mt[3]*dc->cvx[j  ] + mt[4]*dc->cvy[j  ] +mt[5]),
				(int)(mt[0]*dc->cvx[j+1] + mt[1]*dc->cvy[j+1] +mt[2]),
				(int)(mt[3]*dc->cvx[j+1] + mt[4]*dc->cvy[j+1] +mt[5]) );
		}
	}

	fl_color(0, 127, 0);
	fl_line_style( FL_SOLID, 2 );
	for( int i=0; i<ppm->fccnt; i++ ){
		curveintersect *fc = &(ppm->fcurve[i]);
		for( int j=0; j<fc->cvcnt-1; j++ ){
			fl_line((int)(mt[0]*fc->cvx[j] + mt[1]*fc->cvy[j] +mt[2]),
				(int)(mt[3]*fc->cvx[j  ] + mt[4]*fc->cvy[j  ] +mt[5]),
				(int)(mt[0]*fc->cvx[j+1] + mt[1]*fc->cvy[j+1] +mt[2]),
				(int)(mt[3]*fc->cvx[j+1] + mt[4]*fc->cvy[j+1] +mt[5]) );
		}
	}
	fl_color(127, 0, 0);
	for( int i=0; i<ppm->tccnt; i++ ){
		curveintersect *tc = &(ppm->tcurve[i]);
		for( int j=0; j<tc->cvcnt-1; j++ ){
			fl_line((int)(mt[0]*tc->cvx[j] + mt[1]*tc->cvy[j] +mt[2]),
				(int)(mt[3]*tc->cvx[j  ] + mt[4]*tc->cvy[j  ] +mt[5]),
				(int)(mt[0]*tc->cvx[j+1] + mt[1]*tc->cvy[j+1] +mt[2]),
				(int)(mt[3]*tc->cvx[j+1] + mt[4]*tc->cvy[j+1] +mt[5]) );
		}
	}

	fl_line_style( FL_SOLID, 2 );
	// polygons
	if( disp_PLY && plcnt>0 && plvcnt && plx && ply ){
		fl_color(0, 0, 0);
		for( int i=0; i<plcnt; i++ ){
			if( disp_PRI && ppm->pl_cridx[i]>0 ){
				continue;
			}
			if( plvcnt[i]==4 ){
				fl_loop((int)(mt[0] * plx[i*4  ] + mt[1] * ply[i*4  ] +mt[2]),
					//fl_polygon((int)(mt[0] * plx[i*4  ] + mt[1] * ply[i*4  ] +mt[2]),
					(int)(mt[3] * plx[i*4  ] + mt[4] * ply[i*4  ] +mt[5]),
					(int)(mt[0] * plx[i*4+1] + mt[1] * ply[i*4+1] +mt[2]),
					(int)(mt[3] * plx[i*4+1] + mt[4] * ply[i*4+1] +mt[5]),
					(int)(mt[0] * plx[i*4+2] + mt[1] * ply[i*4+2] +mt[2]),
					(int)(mt[3] * plx[i*4+2] + mt[4] * ply[i*4+2] +mt[5]),
					(int)(mt[0] * plx[i*4+3] + mt[1] * ply[i*4+3] +mt[2]),
					(int)(mt[3] * plx[i*4+3] + mt[4] * ply[i*4+3] +mt[5]));
			} else if( plvcnt[i]==3 ){
				fl_loop((int)(mt[0] * plx[i*4  ] + mt[1] * ply[i*4  ] +mt[2]),
					//fl_polygon((int)(mt[0] * plx[i*4  ] + mt[1] * ply[i*4  ] +mt[2]),
					(int)(mt[3] * plx[i*4  ] + mt[4] * ply[i*4  ] +mt[5]),
					(int)(mt[0] * plx[i*4+1] + mt[1] * ply[i*4+1] +mt[2]),
					(int)(mt[3] * plx[i*4+1] + mt[4] * ply[i*4+1] +mt[5]),
					(int)(mt[0] * plx[i*4+2] + mt[1] * ply[i*4+2] +mt[2]),
					(int)(mt[3] * plx[i*4+2] + mt[4] * ply[i*4+2] +mt[5]));
			}
		}
	}

	// 折り線
	if( disp_X ){
		int x0,y0,x1,y1;
		fl_color(0, 0, 0);
		fl_line_style( FL_SOLID, 2 );

		crease *c = &(ppm->crs[0]);
		double *cpx0 = c->Xx2d, *cpy0 = c->Xy2d;
		x0=c->Xxs2d;	y0=c->Xys2d;	x1=cpx0[c->Xsidx];	y1=cpy0[c->Xsidx];
		fl_line((int)(mtofs[0]*x0+mtofs[1]*y0+mtofs[2]), (int)(mtofs[3]*x0+mtofs[4]*y0+mtofs[5]), (int)(mtofs[0]*x1+mtofs[1]*y1+mtofs[2]), (int)(mtofs[3]*x1+mtofs[4]*y1+mtofs[5]));
		for( int i=c->Xsidx; i<c->Xeidx; i++ ){
			x0=cpx0[i];	y0=cpy0[i];	x1=cpx0[i+1];	y1=cpy0[i+1];
			fl_line((int)(mtofs[0]*x0+mtofs[1]*y0+mtofs[2]), (int)(mtofs[3]*x0+mtofs[4]*y0+mtofs[5]), (int)(mtofs[0]*x1+mtofs[1]*y1+mtofs[2]), (int)(mtofs[3]*x1+mtofs[4]*y1+mtofs[5]));
		}
		x0=cpx0[c->Xeidx];	y0=cpy0[c->Xeidx];	x1=c->Xxe2d;	y1=c->Xye2d;
		fl_line((int)(mtofs[0]*x0+mtofs[1]*y0+mtofs[2]), (int)(mtofs[3]*x0+mtofs[4]*y0+mtofs[5]), (int)(mtofs[0]*x1+mtofs[1]*y1+mtofs[2]), (int)(mtofs[3]*x1+mtofs[4]*y1+mtofs[5]));

		if( !disp_PRI ){
			// crease left
			for( int j=1; j<ppm->lcrcnt; j++ ){
				crease *c = ppm->lcrs[j];
				double *cpx0 = c->Xx2d, *cpy0 = c->Xy2d;
				x0=c->Xxs2d;	y0=c->Xys2d;	x1=cpx0[c->Xsidx];	y1=cpy0[c->Xsidx];
				fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
				for( int i=c->Xsidx; i<c->Xeidx; i++ ){
					x0=cpx0[i];	y0=cpy0[i];	x1=cpx0[i+1];	y1=cpy0[i+1];
					fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
				}
				x0=cpx0[c->Xeidx];	y0=cpy0[c->Xeidx];	x1=c->Xxe2d;	y1=c->Xye2d;
				fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
			} // j
			// crease right
			for( int j=1; j<ppm->rcrcnt; j++ ){
				crease *c = ppm->rcrs[j];
				double *cpx0 = c->Xx2d, *cpy0 = c->Xy2d;
				x0=c->Xxs2d;	y0=c->Xys2d;	x1=cpx0[c->Xsidx];	y1=cpy0[c->Xsidx];
				fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
				for( int i=c->Xsidx; i<c->Xeidx; i++ ){
					x0=cpx0[i];	y0=cpy0[i];	x1=cpx0[i+1];	y1=cpy0[i+1];
					fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
				}
				x0=cpx0[c->Xeidx];	y0=cpy0[c->Xeidx];	x1=c->Xxe2d;	y1=c->Xye2d;
				fl_line((int)(mt[0]*x0+mt[1]*y0+mt[2]), (int)(mt[3]*x0+mt[4]*y0+mt[5]), (int)(mt[0]*x1+mt[1]*y1+mt[2]), (int)(mt[3]*x1+mt[4]*y1+mt[5]));
			} // j
		}
	}

	fl_line_style( FL_SOLID, 1 );

	// crease ruling
	if(	disp_R ){
		for( int j=0; j<ppm->lcrcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = ppm->lcrs[j];
			if( j%2==0 ){
				//fl_color(150, 100, 230); // Lavender: B57EDC
				fl_color(82, 32, 118); // Dark Violet: 522076
			} else {
				fl_color(230, 0, 100); // Raspberry: E30B5D
			}
			//fl_color(255, 255, 0); // YELLOW
			fl_line_style( FL_SOLID, 2 );
			for( int i=c->Xsidx; i<c->Xeidx+1; i++ ){
#ifdef EVERYOTHER
				if( i%2 ) continue;
#endif
				if( c->rllen[i]==0.0 ) continue;
				tmpx0 = c->Xx2d[i];
				tmpy0 = c->Xy2d[i];
				tmpx1 = c->Xx2d[i] + c->rlx_cp[i]*c->rllen[i];
				tmpy1 = c->Xy2d[i] + c->rly_cp[i]*c->rllen[i];
				fl_line((int)(mt[0]*tmpx0 + mt[1]*tmpy0 +mt[2]), (int)(mt[3]*tmpx0 + mt[4]*tmpy0 +mt[5]),
					(int)(mt[0]*tmpx1 + mt[1]*tmpy1 +mt[2]), (int)(mt[3]*tmpx1 + mt[4]*tmpy1 +mt[5]) );
			} // i
		} // j
		for( int j=0; j<ppm->rcrcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = ppm->rcrs[j];
			if( j%2==0 ){
				fl_color(230, 0, 100); // Raspberry: E30B5D
			} else {
				//fl_color(150, 100, 230); // Lavender: B57EDC
				fl_color(82, 32, 118); // Dark Violet: 522076
			}
			//fl_color(255, 255, 0); // YELLOW
			fl_line_style( FL_SOLID, 2 );
			for( int i=c->Xsidx; i<c->Xeidx+1; i++ ){
#ifdef EVERYOTHER
				if( i%2 ) continue;
#endif
				if( c->rrlen[i]==0.0 ) continue;
				tmpx0 = c->Xx2d[i];
				tmpy0 = c->Xy2d[i];
				tmpx1 = c->Xx2d[i] + c->rrx_cp[i]*c->rrlen[i];
				tmpy1 = c->Xy2d[i] + c->rry_cp[i]*c->rrlen[i];
				fl_line((int)(mt[0]*tmpx0 + mt[1]*tmpy0 +mt[2]), (int)(mt[3]*tmpx0 + mt[4]*tmpy0 +mt[5]),
					(int)(mt[0]*tmpx1 + mt[1]*tmpy1 +mt[2]), (int)(mt[3]*tmpx1 + mt[4]*tmpy1 +mt[5]) );
			} // i
		} // j
	}
	fl_line_style( FL_SOLID, 1 );

	// Target points
	if( disp_TGT ){
		for( int i=0; i<ppm->tgcnt; i++ )
		{
			int r,g,b;
			if ( i == disp_modTGT ) {
				r = 255; g = 255; b = 0;
			} else if( ppm->tgap[i] >= 0 ){
				r = (int)(ppm->tgap[i]*0.05*256); r = r < 255 ? r : 255;
				g = 255 - r;
				b = 0;
			} else {
				b = (int)(-ppm->tgap[i]*0.05*255); b = b < 255 ? b : 255;
				g = 255 - b;
				r = 0;
			}
			fl_color(r, g, b);
			fl_line_style( FL_SOLID, 2 );
			fl_circle( (int)(mt[0]*ppm->ogx_cp[i]+mt[1]*ppm->ogy_cp[i]+mt[2]),
				(int)(mt[3]*ppm->ogx_cp[i]+mt[4]*ppm->ogy_cp[i]+mt[5]), 3*wscl);

			//glTranslatef( ppm->ogx[i], ppm->ogy[i], ppm->ogz[i] );
			//glutSolidSphere( 1.5, 8, 8 );
		}
	}

	if (disp_stitch) {
		fl_color(0, 127, 0);
		fl_line_style(FL_SOLID, 3);
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < ppm->stpcnt; i++) {
				fl_circle((int)(mt[0] * ppm->stx_cp[j][i] + mt[1] * ppm->sty_cp[j][i] + mt[2]),
					(int)(mt[3] * ppm->stx_cp[j][i] + mt[4] * ppm->sty_cp[j][i] + mt[5]), 2);
			}
		}
	}

	// TN
	if(	disp_TNB ){
		fl_color(255, 0, 0);
		for( int j=0; j<ppm->crcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = &(ppm->crs[j]);
			for( int i=c->Xsidx; i<c->Xeidx+1; i++ ){
				tmpx0 = c->Xx2d[i];
				tmpy0 = c->Xy2d[i];
				tmpx1 = c->Xx2d[i] + c->Tx2d[i]*len0;
				tmpy1 = c->Xy2d[i] + c->Ty2d[i]*len0;
				fl_line((int)(mt[0]*tmpx0 +mt[1]*tmpy0 +mt[2]),
					(int)(mt[3]*tmpx0 +mt[4]*tmpy0 +mt[5]),
					(int)(mt[0]*tmpx1 +mt[1]*tmpy1 +mt[2]),
					(int)(mt[3]*tmpx1 +mt[4]*tmpy1 +mt[5]) );
			} // i
		} // j

		fl_color(0, 255, 0);
		for( int j=0; j<ppm->crcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = &(ppm->crs[j]);
			for( int i=c->Xsidx; i<c->Xeidx+1; i++ ){
				tmpx0 = c->Xx2d[i];
				tmpy0 = c->Xy2d[i];
				tmpx1 = c->Xx2d[i] + c->Nx2d[i]*len0;
				tmpy1 = c->Xy2d[i] + c->Ny2d[i]*len0;
				fl_line((int)(mt[0]*tmpx0 +mt[1]*tmpy0 +mt[2]),
					(int)(mt[3]*tmpx0 +mt[4]*tmpy0 +mt[5]),
					(int)(mt[0]*tmpx1 +mt[1]*tmpy1 +mt[2]),
					(int)(mt[3]*tmpx1 +mt[4]*tmpy1 +mt[5]) );
			} // i
		} // j
	}

	if( disp_CP && !disp_PRI ){
		// control points of crease curve
		for( int j=1; j<ppm->crcnt; j++ ){
			crease *c = &(ppm->crs[j]);
			for( int i=0; i<CCNT; i++ ){
				fl_color(127, 127, 127);
				fl_line_style( FL_SOLID, 3 );
				fl_circle( (int)(mt[0]*c->CPx[i]+mt[1]*c->CPy[i]+mt[2]),
					(int)(mt[3]*c->CPx[i]+mt[4]*c->CPy[i]+mt[5]), 4 );
			}
			fl_line_style( FL_SOLID, 1 );
			fl_color(127, 127, 127);
			for( int i=0; i<CCNT-1; i++ ){
				fl_line((int)(mt[0]*c->CPx[i]+mt[1]*c->CPy[i]+mt[2]),
					(int)(mt[3]*c->CPx[i]+mt[4]*c->CPy[i]+mt[5]),
					(int)(mt[0]*c->CPx[i+1]+mt[1]*c->CPy[i+1]+mt[2]),
					(int)(mt[3]*c->CPx[i+1]+mt[4]*c->CPy[i+1]+mt[5]) );
			}
		}
	}

	// control points
	if( disp_CP && ppos>-1 && pprm>-1 ){
		int Pcnt = ppm->crs[0].Pcnt;
		int Xcnt = ppm->crs[0].Xcnt;
		double *cpx = ppm->crs[0].Xx2d;
		double *cpy = ppm->crs[0].Xy2d;
		switch(pprm){
			case P_CV2D:	fl_color(255, 0, 255);	break;
			case P_CV3D:	fl_color(0, 255, 0);	break;
			case P_TRSN:	fl_color(0, 0, 255);	break;
			case P_FLDA:	fl_color(255, 0, 0);	break;
			case P_RULL:	fl_color(82, 32, 118);	break;	// Dark Violet: 522076
			case P_RULR:	fl_color(230, 0, 100);	break;	// Raspberry: E30B5D
			case P_TRSN1:	fl_color(0, 0, 255);	break;
			case P_FLDA1:	fl_color(255, 0, 0);	break;
			case P_RULL1:	fl_color(82, 32, 118);	break;	// Dark Violet: 522076
			case P_RULR1:	fl_color(230, 0, 100);	break;	// Raspberry: E30B5D
		}
		int ii = (int) ( (double)ppos/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		fl_line_style( FL_SOLID, 5 );
		fl_circle( (int)(mt[0]*cpx[ii]+mt[1]*cpy[ii]+mt[2]),
			(int)(mt[3]*cpx[ii]+mt[4]*cpy[ii]+mt[5]),
			3*wscl);

		int prev, next;
		if( ppos==0 ){ prev=0; } else { prev=ppos-1; }
		if( ppos==Pcnt-1 ){ next=Pcnt-1; } else { next=ppos+1; }
		int i0 = (int) ( (double)prev/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		int i1 = (int) ( (double)next/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		fl_line_style( FL_SOLID, 3 );
		for( int i=i0+1; i<i1; i++ ){
			fl_line((int)(mt[0]*cpx[i-1]+mt[1]*cpy[i-1]+mt[2]),
				(int)(mt[3]*cpx[i-1]+mt[4]*cpy[i-1]+mt[5]),
				(int)(mt[0]*cpx[i  ]+mt[1]*cpy[i  ]+mt[2]),
				(int)(mt[3]*cpx[i  ]+mt[4]*cpy[i  ]+mt[5]) );
		}

		fl_color(127, 127, 127);
		fl_line_style( FL_SOLID, 3 );
		//for( int jj=0; jj<CCNT; jj++ ){
		for( int jj=0; jj<Pcnt; jj++ ){
			if( jj==ppos ){
				continue;
			}
			ii = (int) ( (double)jj/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
			fl_circle( (int)(mt[0]*cpx[ii]+mt[1]*cpy[ii]+mt[2]),
				(int)(mt[3]*cpx[ii]+mt[4]*cpy[ii]+mt[5]), 4 );
		}
	}

	// coordinates
	if( disp_axis ){
		fl_line_style( FL_SOLID, 3 );
		fl_color(255, 0, 0);	fl_line( ox, oy, ox+100*wscl, oy);
		fl_color(0, 255, 0);	fl_line( ox, oy, ox, oy+100*hscl);
	}

	imgbuf = fl_read_image(imgbuf, 0, 0, w(), h());
}

int GraphWindowCP::handle(int event)
{
	int key;

	if( ppm==NULL ){
		return Fl_Double_Window::handle(event);
	}
	int pw = ppm->pw;
	int ph = ppm->ph;
	int psx = ppm->psx;
	int psy = ppm->psy;
	int pex = ppm->pex;
	int pey = ppm->pey;
	double *m2 = ppm->crs[0].m2;
	double *m3 = ppm->crs[0].m3;

	switch(event){
	case FL_FOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
		key = Fl::event_key();
		if(key == 'r'){
			ofs_a += 0.05;
			ofs_cosa = cos( ofs_a );
			ofs_sina = sin( ofs_a );
			redraw();
			return 1;
		} else if(key == 'l'){
			ofs_a -= 0.05;
			ofs_cosa = cos( ofs_a );
			ofs_sina = sin( ofs_a );
			redraw();
			return 1;
		} else if(key == 65363){
			ofs_x += 5;
			redraw();
			return 1;
		} else if(key == 65361){
			ofs_x -= 5;
			redraw();
			return 1;
		} else if(key == 65362){
			ofs_y -= 5;
			redraw();
			return 1;
		} else if(key == 65364){
			ofs_y += 5;
			redraw();
			return 1;
		}
		break;
	case FL_KEYUP:
		key = Fl::event_key();
#if 1		// delete target point
		if (key == 'd' && disp_modTGT > -1) {
			if (disp_modTGT < ppm->tgcnt - 1) {
				int size = ppm->tgcnt - disp_modTGT;
				memcpy(&(ppm->ogx_cp[disp_modTGT]), &(ppm->ogx_cp[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->ogy_cp[disp_modTGT]), &(ppm->ogy_cp[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->ogx[disp_modTGT]), &(ppm->ogx[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->ogy[disp_modTGT]), &(ppm->ogy[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->ogz[disp_modTGT]), &(ppm->ogz[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->tgx[disp_modTGT]), &(ppm->tgx[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->tgy[disp_modTGT]), &(ppm->tgy[disp_modTGT + 1]), sizeof(double) * size);
				memcpy(&(ppm->tgz[disp_modTGT]), &(ppm->tgz[disp_modTGT + 1]), sizeof(double) * size);
			}
			ppm->tgcnt--;
			if (disp_modTGT+1 > ppm->tgcnt) {
				disp_modTGT = ppm->tgcnt - 1;
			}
			redraw();
			((ControlPanel*)cwin)->gwin->redraw();
		}
#endif
		break;
	case FL_PUSH:
		if( Fl::event_button() == FL_RIGHT_MOUSE || Fl::event_button() == FL_LEFT_MOUSE ){
			push_x = Fl::event_x();
			push_y = Fl::event_y();
			if( Fl::event_button() == FL_LEFT_MOUSE ){
				if (ppm!=NULL && -2 < disp_modTGT && disp_modTGT < ppm->tgcnt ) {
					// choose nearest target point
					int thres = 8;
					int minidx = -1;
					double mindist = 100000;
					for (int ii = 0; ii < ppm->tgcnt; ii++ ) {
						double dx = ox + ppm->ogx_cp[ii] * wscl - (double)push_x;
						double dy = oy + ppm->ogy_cp[ii] * hscl - (double)push_y;
						double dist2 = dx * dx + dy * dy;
						if (dist2 < mindist) {
							minidx = ii;
							mindist = dist2;
						}
					}
					if (minidx > -1 && mindist < thres*thres) {
						((ControlPanel*)cwin)->gwin->disp_modTGT = disp_modTGT = minidx;
						flg_psize = 10;
					} // add new target point
					else if(ppm->tgcnt+1 < MAX_TGT_CNT){
						((ControlPanel*)cwin)->gwin->disp_modTGT = disp_modTGT = ppm->tgcnt;
						ppm->ogx_cp[disp_modTGT] = ((double)push_x - ox) / wscl;
						ppm->ogy_cp[disp_modTGT] = ((double)push_y - oy) / hscl;
						ppm->tgcnt++;
						ppm->getTgt2D3D();
						ppm->tgx[disp_modTGT] = ppm->ogx[disp_modTGT];
						ppm->tgy[disp_modTGT] = ppm->ogy[disp_modTGT];
						ppm->tgz[disp_modTGT] = ppm->ogz[disp_modTGT];
						flg_psize = 10;
					}
				}
				else {
					int thres = 8;
					if (abs((int)(ox + psx * wscl) - push_x) < thres) {
						flg_psize = 1;
					}
					else if (abs((int)(ox + pex * wscl) - push_x) < thres) {
						flg_psize = 2;
					}
					else if (abs((int)(oy + psy * hscl) - push_y) < thres) {
						flg_psize = 3;
					}
					else if (abs((int)(oy + pey * hscl) - push_y) < thres) {
						flg_psize = 4;
					}
					else {
						flg_psize = 0;
					}
				}
			}
			//printf( "FL_PUSH: x: %d, y: %d, flg:%d\n", push_x, push_y, flg_psize );
			return 1; // ドラッグイベントを有効にする
		}
		break;
	case FL_DRAG:
		if( Fl::event_button() == FL_RIGHT_MOUSE ){
			ofs_x = Fl::event_x() - push_x;
			ofs_y = Fl::event_y() - push_y;
			//printf( "FL_DRAG: dx: %.2f, dy:%.2f\n", ofs_x, ofs_y );
			redraw();
		}
		if( Fl::event_button() == FL_LEFT_MOUSE ){
			switch( flg_psize ){
	case 1:
		ofs_psx = Fl::event_x() - push_x;
		break;
	case 2:
		ofs_pex = Fl::event_x() - push_x;
		break;
	case 3:
		ofs_psy = Fl::event_y() - push_y;
		break;
	case 4:
		ofs_pey = Fl::event_y() - push_y;
		break;
	case 10:
		if ( -1 < disp_modTGT && disp_modTGT < ppm->tgcnt ) {
			ppm->ogx_cp[disp_modTGT] = ((double)Fl::event_x() - ox) / wscl;
			ppm->ogy_cp[disp_modTGT] = ((double)Fl::event_y() - oy) / hscl;
		}
		break;
	default:
		ofs_a = (Fl::event_y() - push_y) * M_PI/this->h();
		ofs_sina = sin(ofs_a);
		ofs_cosa = cos(ofs_a);
		break;
			}
			if( 1<=flg_psize && flg_psize<5 ){
				//printf( "FL_DRAG: psx: %d, psy: %d, pex: %d, pey: %d, da: %.2f\n",
				//	ofs_psx, ofs_psy, ofs_pex, ofs_pey, ofs_a*180.0/M_PI );
			}
			redraw();
		}
		break;
	case FL_RELEASE:
		if( Fl::event_button() == FL_RIGHT_MOUSE ){
			ofs_x = Fl::event_x() - push_x;
			ofs_y = Fl::event_y() - push_y;
			if( m2 ){
				m2[6] += ofs_x/wscl;
				m2[7] += ofs_y/hscl;
			}
			if( m3 ){
				m3[12] += ofs_x/wscl;
				m3[13] += ofs_y/hscl;
			}
			if( ppm ){
				ppm->set_postproc_type( PPTYPE_RTCURVE );
				ppm->crs[0].calcXA_CP( ppm->flg_interpolate, &ppm->rp );
				ppm->postproc();
				ppm->set_postproc_type( PPTYPE_UNDEF );
			}
			//printf( "FL_RELEASE: dx: %.2f, dy: %.2f\n", ofs_x, ofs_y );
			ofs_x = ofs_y = 0.0;
			redraw();
			((ControlPanel*)cwin)->gwin->redraw();
			((ControlPanel*)cwin)->gwin_gr->redraw();
		}
		if( Fl::event_button() == FL_LEFT_MOUSE ){
			switch( flg_psize ){
	case 1:
		ofs_psx = Fl::event_x() - push_x;
		if( ppm ){
			ppm->psx += (int)(ofs_psx/wscl);
			ppm->pw = ppm->pex - ppm->psx;
			ppm->set_postproc_type( PPTYPE_PPEDGE );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
		}
		//printf( "FL_RELEASE: psx: %d\n", ofs_psx );
		ofs_psx = 0;
		break;
	case 2:
		ofs_pex = Fl::event_x() - push_x;
		if( ppm ){
			ppm->pex += (int)(ofs_pex/wscl);
			ppm->pw = ppm->pex - ppm->psx;
			ppm->set_postproc_type( PPTYPE_PPEDGE );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
		}
		//printf( "FL_RELEASE: pex: %d\n", ofs_pex );
		ofs_pex = 0;
		break;
	case 3:
		ofs_psy = Fl::event_y() - push_y;
		if( ppm ){
			ppm->psy += (int)(ofs_psy/hscl);
			ppm->ph = ppm->pey - ppm->psy;
			ppm->set_postproc_type( PPTYPE_PPEDGE );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
		}
		//printf( "FL_RELEASE: psy: %d\n", ofs_psy );
		ofs_psy = 0;
		break;
	case 4:
		ofs_pey = Fl::event_y() - push_y;
		if( ppm ){
			ppm->pey += (int)(ofs_pey/hscl);
			ppm->ph = ppm->pey - ppm->psy;
			ppm->set_postproc_type( PPTYPE_PPEDGE );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
		}
		//printf( "FL_RELEASE: pey: %d\n", ofs_pey );
		ofs_pey = 0;
		break;
	case 10:
		if (-1 < disp_modTGT && disp_modTGT < ppm->tgcnt) {
			ppm->ogx_cp[disp_modTGT] = ((double)Fl::event_x() - ox) / wscl;
			ppm->ogy_cp[disp_modTGT] = ((double)Fl::event_y() - oy) / hscl;
			ppm->getTgt2D3D();
		}
		break;
	default:
		ofs_a = (Fl::event_y() - push_y) * M_PI/this->h();
		if( m2 ){
			double sina = sin(ofs_a), cosa = cos(ofs_a);
			double m00, m01, m10, m11;
			m00 = m2[0]*cosa - m2[1]*sina;
			m01 = m2[0]*sina + m2[1]*cosa;
			m10 = m2[3]*cosa - m2[4]*sina;
			m11 = m2[3]*sina + m2[4]*cosa;
			m2[0] = m00;
			m2[1] = m01;
			m2[3] = m10;
			m2[4] = m11;
		}
		if( m3 ){
			double sina = sin(ofs_a), cosa = cos(ofs_a);
			double m00, m01, m10, m11;
			m00 = m3[0]*cosa - m3[1]*sina;
			m01 = m3[0]*sina + m3[1]*cosa;
			m10 = m3[4]*cosa - m3[5]*sina;
			m11 = m3[4]*sina + m3[5]*cosa;
			m3[0] = m00;
			m3[1] = m01;
			m3[4] = m10;
			m3[5] = m11;
		}
		if( ppm ){
			ppm->set_postproc_type( PPTYPE_RTCURVE );
			ppm->crs[0].calcXA_CP( ppm->flg_interpolate, &ppm->rp );
			ppm->postproc();
			ppm->set_postproc_type( PPTYPE_UNDEF );
		}
		//printf( "FL_RELEASE: da: %.2f\n", ofs_a*180.0/M_PI );
		ofs_a = ofs_sina = 0.0;
		ofs_cosa = 1.0;
		break;
			}
			redraw();
			((ControlPanel*)cwin)->gwin->redraw();
			((ControlPanel*)cwin)->gwin_gr->redraw();
		}
		return 1;
		break;
	case FL_MOVE:
		break;
	default:
		break;
	}
	return Fl_Double_Window::handle(event);
}

void GraphWindowCP::exportImage(char* filename)
{
	cv::Mat img = cv::Mat(h(), w(), CV_8UC3, imgbuf, w() * 3);
	cv::imwrite(filename, img);
}

void GraphWindowCP::exportCroppedImage(char* filename)
{
	cv::Mat img = cv::Mat(h(), w(), CV_8UC3, imgbuf, w() * 3);
	int pw = (int)(ppm->pw * wscl), ph = (int)(ppm->ph * hscl);
	int psx = (int)(ppm->psx * wscl), psy = (int)(ppm->psy * hscl), pex = (int)(ppm->pex * wscl), pey = (int)(ppm->pey * hscl);
	cv::Rect roi(cv::Point((int)(ox + psx + ofs_psx), (int)(oy + psy + ofs_psy)), cv::Size((int)(pw + ofs_pex - ofs_psx), (int)(ph + ofs_pey - ofs_psy)));
	cv::Mat crpimg = img(roi);
	cv::imwrite(filename, crpimg);
}
