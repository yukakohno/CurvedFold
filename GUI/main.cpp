
#if 1
#pragma comment(lib,"opengl32.lib")
//#pragma comment(lib,"glu32.lib")
//#pragma comment(lib,"glut32.lib")
#pragma comment(lib,"freeglut.lib")
#pragma comment(lib,"fltk.lib")
#pragma comment(lib,"fltkforms.lib")
#pragma comment(lib,"fltkgl.lib")
#pragma comment(lib,"fltkimages.lib")
#pragma comment(lib,"fltkjpeg.lib")
#pragma comment(lib,"fltkpng.lib")
#pragma comment(lib,"fltkzlib.lib")
#endif

#ifdef _DEBUG
#pragma comment(lib,"opencv_world420d.lib")
#else
#pragma comment(lib,"opencv_world420.lib")
#endif

#include "GraphWindow3DCF.h"
#include "GraphWindowParam.h"
#include "ControlPanel.h"

#define GCF_W 480
#define GCF_H 480
#define CP_W 200
#define LBL_W 50
#define GCP_W 320

GraphWindow3DCF *gwin=NULL;
ControlPanel *cwin=NULL;
GraphWindowCP *gwin_cp=NULL;
GraphWindowParam *gwin_gr=NULL;

int main(int argc, char* argv[])
{
	Fl_Window *window = new Fl_Window( 5, 50, GCF_W + GCP_W + LBL_W + GCP_W + CP_W, GCF_H);

	gwin = new GraphWindow3DCF(0, 0, GCF_W, GCF_H);
	gwin->end();

	gwin_cp = new GraphWindowCP(GCF_W, 0, GCP_W, GCF_H);
	gwin_cp->end();

	unsigned char col[9][3]={ {255,0,255}, {0,200,0}, {127,127,255}, {255,0,0},
	{255,0,0}, {82, 32, 118}/*{150,100,230}*/, {230,0,100},
	{255,0,0}, {150,100,230}};
	int gr1_hgt = GCF_H/5, bx_y = gr1_hgt/2;
	Fl_Box *bx_bg= new Fl_Box(GCF_W+GCP_W, 0, LBL_W, GCF_H,"");	bx_bg->box( FL_FLAT_BOX );	bx_bg->color( FL_WHITE );
	Fl_Box *bx_lbl0= new Fl_Box(GCF_W+GCP_W, bx_y-40, LBL_W, 40, "curv\n2D");	bx_lbl0->labelcolor( fl_rgb_color(col[0][0],col[0][1],col[0][2]));
	Fl_Box *bx_lbl1= new Fl_Box(GCF_W+GCP_W, bx_y   , LBL_W, 40, "curv\n3D");	bx_lbl1->labelcolor( fl_rgb_color(col[1][0],col[1][1],col[1][2]));	bx_y+=gr1_hgt;
	Fl_Box *bx_lbl2= new Fl_Box(GCF_W+GCP_W, bx_y-40, LBL_W, 20, "torsion");	bx_lbl2->labelcolor( fl_rgb_color(col[2][0],col[2][1],col[2][2]));
	Fl_Box *bx_lbl3= new Fl_Box(GCF_W+GCP_W, bx_y-20, LBL_W, 60, "diff.\nfold\nangle");	bx_lbl3->labelcolor( fl_rgb_color(col[3][0],col[3][1],col[3][2]));	bx_y+=gr1_hgt;
	Fl_Box *bx_lbl4= new Fl_Box(GCF_W+GCP_W, bx_y-20, LBL_W, 40, "fold\nangle");	bx_lbl4->labelcolor( fl_rgb_color(col[4][0],col[4][1],col[4][2]));	bx_y+=gr1_hgt;
	Fl_Box *bx_lbl5= new Fl_Box(GCF_W+GCP_W, bx_y-60, LBL_W, 60, "ruling\nangle\n(Left)");	bx_lbl5->labelcolor( fl_rgb_color(col[5][0],col[5][1],col[5][2]));
	Fl_Box *bx_lbl6= new Fl_Box(GCF_W+GCP_W, bx_y   , LBL_W, 60, "ruling\nangle\n(Right)");	bx_lbl6->labelcolor( fl_rgb_color(col[6][0],col[6][1],col[6][2]));	bx_y+=gr1_hgt;
	Fl_Box *bx_lbl7= new Fl_Box(GCF_W+GCP_W, bx_y-10, LBL_W, 20, "error");	bx_lbl7->labelcolor( fl_rgb_color(col[7][0],col[7][1],col[7][2]));

	gwin_gr = new GraphWindowParam(GCF_W+GCP_W+LBL_W, 0, GCP_W, GCF_H );
	gwin_gr->end();

	cwin = new ControlPanel(GCF_W+GCP_W+LBL_W+GCP_W, 0, CP_W, GCF_H, "aaa", gwin);
	cwin->end();

	window->add(gwin);
	//window->add(bx_bg);
	window->add(gwin_cp);
	window->add(gwin_gr);
	window->add(cwin);
	window->end();
	window->show();

	gwin_cp->ppm = &(cwin->ppm);
	gwin_cp->disp_axis = gwin->disp_axis;
	gwin_cp->disp_X = gwin->disp_X;
	gwin_cp->disp_TNB = gwin->disp_TNB;
	gwin_cp->disp_R = gwin->disp_R;
	gwin_cp->disp_PLY = gwin->disp_PLY;
	gwin_cp->disp_PTN = gwin->disp_PTN;
	gwin_cp->disp_CP = gwin->disp_CP;
	gwin_cp->ppos = gwin->ppos;
	gwin_cp->pprm = gwin->pprm;
	gwin_gr->ppm = &(cwin->ppm);
	cwin->gwin = gwin;
	cwin->gwin_cp = gwin_cp;
	cwin->gwin_gr = gwin_gr;
	gwin_cp->cwin = (void*)cwin;

	int ret = cwin->ppm.crs[0].load("input/P.txt");
	ret = cwin->ppm.crs[0].loadm2m3("input/m2m3.txt");
	cwin->rb_fix[CMODE_A]->setonly();
	cwin->ppm.set_postproc_type( PPTYPE_OPEN );
	cwin->refresh(1);
	cwin->ppm.set_postproc_type( PPTYPE_UNDEF );

	crease *c0=&(cwin->ppm.crs[0]);
	cwin->vs_ppos->bounds( 0, c0->Pcnt-1 );
	Fl::run();

	return 0;
}

