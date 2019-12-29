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
}

#define MAXV 7
void GraphWindowParam::draw()
{
	int i,j, gflg[MAXV]={1,1,1,1,1,1,1};
	int gsx[MAXV], gsy[MAXV], gw[MAXV], gh[MAXV], goy[MAXV], lw[MAXV];
	double *gcp[MAXV], *gv[MAXV], go[MAXV], gc[MAXV];
	unsigned char col[MAXV][3]={ {255,0,255}, {0,200,0}, {127,127,255}, {255,0,0},
	{255,0,0}, {82, 32, 118}/*{150,100,230}*/, {230,0,100}};
	crease *c = NULL;
	int Pcnt = ppm->crs[0].Pcnt;
	int Xcnt = ppm->crs[0].Xcnt;

	fl_color(FL_WHITE);
	fl_rectf(0,0,w(),h());

	// TEST
	//const char text[] = "string";
	//fl_draw(text, 10,10);

	if( ppm==NULL ){
		return;
	}

	for( i=0; i<MAXV; i++ ){
		gsx[i] = 5;
		gw[i] = w()-10;
		gh[i] = h()/5-10;
	}

	gsy[0] = gsy[1] = 0;
	gsy[2] = gsy[3] = gsy[0] + gh[0] +10;
	gsy[4] = gsy[2] + gh[2] +10;
	gsy[5] = gsy[6] = gsy[4] + gh[4] +10;

	goy[0] = goy[1] = gsy[0]+gh[0]/2;
	goy[2] = goy[3] = gsy[2]+gh[2]/2;
	goy[4] = gsy[4]+gh[4]/2;
	goy[5] = goy[6] = gsy[5]+gh[5]/2;

	for( i=0; i<MAXV; i++ ){ lw[i]=2; }
	for( i=0; i<MAXV; i++ ){ gcp[i]=NULL; }
	for( i=0; i<MAXV; i++ ){ gv[i]=NULL; }
	c = &(ppm->crs[0]);
	gcp[0] = c->Px2d;
	gcp[1] = c->Px;
	gcp[2] = c->Py;
	gcp[4] = c->Pa;
	gv[0] = c->k2d;		go[0] = goy[0];			gc[0] = -gh[0]/0.1;
	gv[1] = c->kv;		go[1] = goy[1];			gc[1] = -gh[1]/0.1;		lw[1]=1;
	gv[2] = c->tr;		go[2] = goy[2];			gc[2] = -gh[2]/0.1;
	gv[3] = c->da;		go[3] = goy[3];			gc[3] = -gh[3]/0.1;		lw[3]=1;
	gv[4] = c->alpha;	go[4] = goy[4];			gc[4] = -gh[4]/M_PI;
	gv[5] = c->betal;	go[5] = gsy[5]+gh[5];	gc[5] = -gh[5]/M_PI;
	gv[6] = c->betar;	go[6] = gsy[6]+gh[6];	gc[6] = -gh[6]/M_PI;

#if 0	// labels -> ‚È‚º‚©•\Ž¦‚³‚ê‚È‚¢
	char label[10][64] = { "0.1", "-0.1", "90", "0", "90", "-90", "", "", "", "" };
	fl_draw( label[0], strlen(label[0]), 5, gsy[0]);
	fl_draw( label[1], strlen(label[1]), 5, gsy[0]+gh[0]-10);
	fl_draw( label[2], strlen(label[2]), 5, gsy[3]);
	fl_draw( label[3], strlen(label[3]), 5, gsy[3]+gh[3]-10);
	fl_draw( label[4], strlen(label[4]), 5, gsy[4]);
	fl_draw( label[5], strlen(label[5]), 5, gsy[4]+gh[4]-10);
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

