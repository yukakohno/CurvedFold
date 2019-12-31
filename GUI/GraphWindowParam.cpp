#include <FL/fl_draw.H>
#define _USE_MATH_DEFINES
#include <math.h>
#include "GraphWindowParam.h"
#include "../CurvedFoldModel/crease.h"

GraphWindowParam::GraphWindowParam(int X, int Y, int W, int H, const char *L) : Fl_Double_Window(X, Y, W, H, L)
{
	init();
}

GraphWindowParam::GraphWindowParam(int X, int Y, int W, int H) : Fl_Double_Window(X, Y, W, H)
{
	init();
}

void GraphWindowParam::init()
{
	ppm = NULL;
	cidx = 0;
}

#define MAXV 9
void GraphWindowParam::draw()
{
	int i,j, gflg[MAXV]={1,1,1,1,1,1,1,0,0};
	int gsx[MAXV], gsy[MAXV], gw[MAXV], gh[MAXV], goy[MAXV], lw[MAXV];
	double *gcp[MAXV], *gv[MAXV], go[MAXV], gc[MAXV];
	unsigned char col[MAXV][3]={ {255,0,255}, {0,200,0}, {0,0,255}/*{127,127,255}*/, {255,0,0},
	{255,0,0}, {82, 32, 118}/*{150,100,230}*/, {230,0,100},
	{255,0,0}, {150,100,230}};//, {150,100,230}, {230,0,100}, {230,0,100}};
	crease *c = NULL;
	int Pcnt = ppm->crs[0].Pcnt;
	int Xcnt = ppm->crs[0].Xcnt;

	fl_color(FL_WHITE);
	fl_rectf(0,0,w(),h());

	if( ppm==NULL ){
		return;
	}

	for( i=0; i<MAXV; i++ ){
		gsx[i] = 5;
		gw[i] = w()-10;
		gh[i] = h()/5-10;
	}

	gsy[0] = gsy[1] = 0;
	gsy[3] = gsy[0] + gh[0] +10;
	gsy[4] = gsy[3] + gh[3] +10;
	gsy[2] = gsy[4] + gh[4] +10;
	gsy[7] = gsy[8] = gsy[4] + gh[4] +10;
	gsy[5] = gsy[6] = gsy[8] + gh[8] +10;

	//goy[0] = goy[1] = gsy[0]+gh[0]/2;
	goy[0] = goy[1] = gsy[0]+gh[0];
	goy[3] = gsy[3]+gh[3]/2;
	goy[4] = gsy[4]+gh[4]/2;
	//goy[4] = gsy[4] + gh[4];
	goy[2] = gsy[2] + gh[2]/2;
	goy[5] = goy[6] = gsy[5]+gh[5]/2;
	goy[7] = goy[8] = gsy[7] + gh[7];

	for( i=0; i<MAXV; i++ ){ lw[i]=2; }
	for( i=0; i<MAXV; i++ ){ gcp[i]=NULL; }
	for( i=0; i<MAXV; i++ ){ gv[i]=NULL; }
	if( -ppm->lcrcnt<cidx && cidx<ppm->rcrcnt ){
		if( cidx==0 ){
			c = &(ppm->crs[0]);
			gcp[0] = c->Px2d;
			gcp[1] = c->Px;
			gcp[2] = c->Py;
			gcp[4] = c->Pa;
			gcp[5] = c->Pbl;
			gcp[6] = c->Pbr;
		} else if( cidx<0 ){
			c = ppm->lcrs[-cidx];
		} else if( 0<cidx ){
			c = ppm->rcrs[cidx];
		}
		gv[0] = c->k2d;		go[0] = goy[0];			gc[0] = -gh[0]/0.05;
		gv[1] = c->kv;		go[1] = goy[1];			gc[1] = -gh[1]/0.05;		lw[1]=1;
		gv[2] = c->tr;		go[2] = goy[2];			gc[2] = -gh[2]/0.05;
		gv[3] = c->da;		go[3] = goy[3];			gc[3] = -gh[3]/0.05;		lw[3]=1;
		gv[4] = c->alpha;	go[4] = goy[4];			gc[4] = -gh[4]/M_PI;		//lw[4]=1;
		gv[5] = c->betal;	go[5] = gsy[5]+gh[5];	gc[5] = -gh[5]/M_PI;
		gv[6] = c->betar;	go[6] = gsy[6]+gh[6];	gc[6] = -gh[6]/M_PI;
		//gv[7] = ppm->cverr;		go[7] = goy[7];		gc[7] = -gh[7];
		//gv[8] = ppm->cvminerr;	go[8] = goy[8];		gc[8] = -gh[8];			lw[8]=1;
	}

#if 0
	{
		crease *c0 = &(ppm->crs[0]);
		int sx=gsx[4];
		double o=go[4], c=gc[4], x, dxP = (double)gw[4]/(double)(Pcnt-1);
		//int k=20;
		for( int k=c0->FM_fidx_min; k<=c0->FM_fidx_max; k++ )
		{
			if( c0->flg_FM_Pa[k]>0 ){
				fl_line_style( FL_SOLID, 2 );
				fl_color( col[4][0], col[4][1], col[4][2] );
			} else {
				fl_line_style( FL_SOLID, 1 );
				fl_color( FL_GRAY );
			}
			for( j=0, x=sx; j<Pcnt-1; j++, x+=dxP ){
				fl_line( (int)x, (int)(o+ c0->FM_Pa[k][j]*c), (int)(x+dxP), (int)(o+ c0->FM_Pa[k][j+1]*c) );
			}
			for( j=0, x=sx; j<Pcnt; j++, x+=dxP ){
				fl_circle( (int)x, (int)(o+ c0->FM_Pa[k][j]*c), 3 );
			}
		}
	}
#endif
#if 0
	{
		crease *c0 = &(ppm->crs[0]);
		int sx=gsx[7];
		double o=go[7], c=-gh[7]/0.01, x;
		double dxP = (double)gw[7]/(double)(Pcnt-1);
		double dxX = (double)gw[7]/(double)(Xcnt-1);
#if 0
		for( int k=c0->FM_fidx_min; k<=c0->FM_fidx_max; k++ )
		{
			if( c0->flg_FM_Pt[k]>0 ){
				fl_line_style( FL_SOLID, 2 );
				fl_color( col[2][0], col[2][1], col[2][2] );
			} else {
				fl_line_style( FL_SOLID, 1 );
				fl_color( FL_GRAY );
			}
			for( j=0, x=sx; j<Pcnt-1; j++, x+=dxP ){
				double t0 = log(fabs( c0->FM_Pt[k][j]  *100) +1.0 )*0.01;
				double t1 = log(fabs( c0->FM_Pt[k][j+1]*100) +1.0 )*0.01;
				fl_line( (int)x, (int)(o+t0*c), (int)(x+dxP), (int)(o+t1*c) );
			}
			for( j=0, x=sx; j<Pcnt; j++, x+=dxP ){
				double t = log(fabs( c0->FM_Pt[k][j]  *100) +1.0 )*0.01;
				fl_circle( (int)x, (int)(o+ t*c), 3 );
			}
		}
#endif
		fl_color( col[2][0], col[2][1], col[2][2] );
		fl_line_style( FL_SOLID, 2 );
		for( j=0, x=sx; j<Xcnt-1; j++, x+=dxX ){
			double t0 = log(fabs( c0->tr[j]  *100) +1.0 )*0.01;
			double t1 = log(fabs( c0->tr[j+1]*100) +1.0 )*0.01;
			fl_line( (int)x, (int)(o+t0*c), (int)(x+dxX), (int)(o+t1*c) );
		}
		fl_line_style( FL_SOLID, 3 );
		for( j=0, x=sx; j<Pcnt; j++, x+=dxP ){
			double t = log(fabs( c0->Py[j]  *100) +1.0 )*0.01;
			fl_circle( (int)x, (int)(o+ t*c), 4 );
		}
	}
#endif

	for( i=0; i<MAXV; i++ )
	{
		int sx, sy, w, h, oy;
		double *cp, *v, o, c, x, dxX, dxP;

		if( gflg[i]==0 ){
			continue;
		}

		dxX = (double)gw[i]/(double)(Xcnt-1);
		dxP = (double)gw[i]/(double)(Pcnt-1);
		sx = gsx[i];	sy= gsy[i];	w = gw[i];	h = gh[i];	oy = goy[i];
		cp = gcp[i];	v = gv[i];	o = go[i];	c = gc[i];

		fl_color( FL_GRAY );
		fl_line_style( FL_SOLID, 1 );
		fl_rect( sx, sy, w, h );
		fl_line( sx, oy, sx+w, oy );

		if( v ){
			fl_color( col[i][0], col[i][1], col[i][2] );
			//fl_line_style( FL_SOLID, lw[i] );
			fl_line_style( FL_SOLID, 2 );
			if( lw[i]>1 ){
				for( j=0, x=sx; j<Xcnt-1; j++, x+=dxX ){
					fl_line( (int)x, (int)(o+v[j]*c), (int)(x+dxX), (int)(o+v[j+1]*c) );
				}
			} else {
				for( j=0, x=sx; j<Xcnt-1; j+=2, x+=2*dxX ){
					fl_line( (int)x, (int)(o+v[j]*c), (int)(x+dxX), (int)(o+v[j+1]*c) );
				}
			}
		}

		if( cp ){
			fl_line_style( FL_SOLID, 3 );
			for( j=0, x=sx; j<Pcnt; j++, x+=dxP ){
				fl_circle( (int)x, (int)(o+cp[j]*c), 4 );
			}
		}
	}
}

int GraphWindowParam::handle(int event)
{
	return Fl_Double_Window::handle(event);
}

