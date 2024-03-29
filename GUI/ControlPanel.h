
/*************************************************************
	ControlPanel.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Window.h>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Group.h>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Round_Button.h>
#include <FL/Fl_Value_Slider.h>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Tabs.H>
#include "GraphWindow3DCF.h"
#include "GraphWindowCP.h"
#include "GraphWindowParam.h"
#include "../CurvedFoldModel/papermodel.h"
#include "common.h"

#define RECDIR "input/"

//#define TEXFNAME "texture/no_tex.jpg"
#define TEXFNAME "texture/grid_bw3.jpg"
//#define TEXFNAME "texture/cp.jpg"
#define TEXWIDTH  256	// texture width
#define TEXHEIGHT 256	// texture height

#define MAX_HISTORY 100

//#define EVERYOTHER	// ruling を１本おきに表示

#ifndef CPANEL
#define CPANEL

enum DISP { D_X, D_TNB, D_R, D_RLEN, D_PLY, D_PTN, D_CP, D_AX2, D_ONE,
			D_LSMT, D_OFST, D_TGT };

class ControlPanel : public Fl_Window
{
public:
	GraphWindow3DCF *gwin;
	GraphWindowCP *gwin_cp;
	GraphWindowParam *gwin_gr;

	papermodel ppm;
	char filepath[256];

	// for animation
	double kv_org[MAX_SPCNT], k2d_org[MAX_SPCNT], tr_org[MAX_SPCNT], alpha_org[MAX_SPCNT];
	double kv_dif0[MAX_SPCNT], k2d_dif0[MAX_SPCNT], tr_dif0[MAX_SPCNT], alpha_dif0[MAX_SPCNT];
	double kv_dif1[MAX_SPCNT], k2d_dif1[MAX_SPCNT], tr_dif1[MAX_SPCNT], alpha_dif1[MAX_SPCNT];
	int an_fcnt;

	double kmin,kmax,kstep,kstepk, tmin,tmax,tstep,tstepk, amin,amax,astep,astepk;

	ControlPanel(int X, int Y, int W, int H, GraphWindow3DCF *_gwin );
	ControlPanel(int X, int Y, int W, int H, const char *L, GraphWindow3DCF *_gwin );
	~ControlPanel();
	void initPanel();
	void createPanel();
	void refresh(int init);
	//void refresh2();
	void setorg();

private:
	int handle(int event);
	int value_grpfix();
	int value_grpparam();

	// ------------------------- fILE_IO -------------------------------------------
public:
	Fl_Button *btn_load;
	Fl_File_Chooser *fc;
	Fl_Button *btn_loadtex;
	Fl_Button *btn_loadtpt;
	Fl_File_Chooser* fc_tpt;
	char fname_tpt[1024];
	Fl_Button* btn_loadtptmask;
	Fl_File_Chooser* fc_tptmask;
	Fl_Button *btn_loadrul;
	Fl_File_Chooser* fc_rul;
	Fl_Button* btn_loadall;
	Fl_File_Chooser* fc_all;

	Fl_Button *btn_savelog;
	Fl_Button *btn_savescreen;
	Fl_Button *btn_savetpt;
	Fl_Button *btn_saverul;
	Fl_Button* btn_saveall;

private:
	static void cb_btn_load( Fl_Widget *wgt, void *idx);
	static void cb_btn_savelog( Fl_Widget *wgt, void *idx);
	static void cb_btn_savescreen( Fl_Widget *wgt, void *idx);
	static void cb_btn_loadrul( Fl_Widget *wgt, void *idx);
	static void cb_btn_saverul( Fl_Widget *wgt, void *idx);
	static void cb_btn_loadtpt( Fl_Widget *wgt, void *idx);
	static void cb_btn_loadtptmask(Fl_Widget* wgt, void* idx);
	static void cb_btn_savetpt( Fl_Widget *wgt, void *idx);
	static void cb_btn_loadall(Fl_Widget* wgt, void* idx);
	static void cb_btn_saveall(Fl_Widget* wgt, void* idx);

	// ------------------------- EVALUATE -------------------------------------------
public:
	Fl_Button *btn_eval_rulingcross;

private:
	static void cb_btn_eval_rulingcross( Fl_Widget *wgt, void *idx);

	// ------------------------- DISPLAY -------------------------------------------
public:
	Fl_Button *btn_resetRot;
	Fl_Button *btn_resetScale;
	Fl_Button *btn_resetTrans;
	Fl_Check_Button *cb_dispaxis;
	Fl_Check_Button *cb_disp[20]; /* enum DISP */
	Fl_Check_Button *cb_CurveEnd;
private:
	static void cb_cb_dispaxis(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispX(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispTNB(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispR(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispRlen(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispAx2(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispPLY(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispPTN(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispCP(Fl_Widget *wgt, void *idx);
	static void cb_cb_CurveEnd(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispLSMT(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispOFST(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispTGT(Fl_Widget *wgt, void *idx);

	// ------------------------- CONFIGURATION -------------------------------------------
public:
	Fl_Value_Slider *vs_rulres;

private:
	static void cb_vs_rulres(Fl_Widget *wgt, void* idx);

	// ------------------------- DESIGN -------------------------------------------
public:
	// MODE
	Fl_Group *grp_fix;
	Fl_Round_Button *rb_fix[N_C_MODE];

	// PARAMETER
	Fl_Group *grp_param0, *grp_param1;
	Fl_Round_Button *rb_param[N_P_IDX];
	Fl_Value_Slider *vs_ppos;
	Fl_Value_Slider *vs_pval;
	Fl_Button *btn_R2TA;
	Fl_Button *btn_RPar;

	Fl_Value_Slider *vs_fmot;
	Fl_Button *btn_apply;

private:
	static void cb_rb_gwin(Fl_Widget *wgt, void *idx);
	static void cb_rb_fix0(Fl_Widget *wgt, void* idx);
	static void cb_rb_param0(Fl_Widget *wgt, void* idx);
	static void cb_vs_ppos(Fl_Widget *wgt, void *idx);
	static void cb_vs_pval(Fl_Widget *wgt, void *idx);
	static void cb_btn_R2TA(Fl_Widget *wgt, void *idx);
	static void cb_btn_RPar(Fl_Widget *wgt, void *idx);
	static void cb_btn_apply(Fl_Widget *wgt, void *idx);

	// ------------------------- RULING2CURVE --------------------------------------
public:
	Fl_Value_Slider *vs_xmang0;
	Fl_Value_Slider *vs_xmang1;
	Fl_Button *btn_switchRuling;
	Fl_Button *btn_R2TA0;
	Fl_Button *btn_optmat;
	Fl_Button *btn_start;
	Fl_Button *btn_stop;
	Fl_Check_Button *cb_optmat;
	Fl_Check_Button* cb_optrul;

private:
	int fcnt;
	int acnt;
	bool acnt_inc;
	bool flg_idle_active;
	static void cb_btn_R2TA0(Fl_Widget *wgt, void *idx);
	static void cb_btn_switchRuling(Fl_Widget* wgt, void* idx);
	static void cb_vs_xmang(Fl_Widget *wgt, void *idx);
	static void idle(void *idx);
	static void cb_btn_optmat(Fl_Widget *wgt, void *idx);
	static void cb_btn_start(Fl_Widget *wgt, void *idx);
	static void cb_btn_stop(Fl_Widget *wgt, void *idx);

	// ------------------------- OPTIMIZATION --------------------------------------
public:
	Fl_Button* btn_optfold2;
	Fl_Button* btn_randrul;
	Fl_Button* btn_randrul2;

private:
	static void cb_btn_optfold2(Fl_Widget* wgt, void* idx);
	static void cb_btn_randrul(Fl_Widget* wgt, void* idx);
	static void cb_btn_randrul2(Fl_Widget* wgt, void* idx);

	// ------------------------- HISTORY --------------------------------------
public:
	Fl_Check_Button* cb_history;
	Fl_Value_Slider* vs_history;

private:
	static void cb_vs_history(Fl_Widget* wgt, void* idx);

	int hist_head, hist_size;
	int hist_mode[MAX_HISTORY], hist_cpcnt[MAX_HISTORY];
	double histPx2d[MAX_HISTORY][MAX_CPCNT], histPy2d[MAX_HISTORY][MAX_CPCNT];
	double histPx[MAX_HISTORY][MAX_CPCNT], histPy[MAX_HISTORY][MAX_CPCNT], histPz[MAX_HISTORY][MAX_CPCNT];
	double histPa[MAX_HISTORY][MAX_CPCNT];
	double histPbl[MAX_HISTORY][MAX_CPCNT], histPbr[MAX_HISTORY][MAX_CPCNT];
	double hist_m3[MAX_HISTORY][16];
	int push_hist(int cpcnt, int mode, double* Px2d, double* Py2d,
		double* Px, double* Py, double* Pz,
		double* Pa, double* Pbl, double* Pbr, double* _m3);
	int access_hist(int idx, int* cpcnt, int* mode, double* Px2d, double* Py2d,
		double* Px, double* Py, double* Pz,
		double* Pa, double* Pbl, double* Pbr, double* _m3);

	// ------------------------- RECTIFY -------------------------------------------
public:
	// 補正あり／なし
	Fl_Check_Button *cb_rectifyA;
	Fl_Check_Button *cb_rectifyT;
	Fl_Check_Button *cb_rectifyR;

private:
	static void cb_cb_rectifyA(Fl_Widget *wgt, void *idx);
	static void cb_cb_rectifyT(Fl_Widget *wgt, void *idx);
	static void cb_cb_rectifyR(Fl_Widget *wgt, void *idx);

	// ------------------------- FOLD_MOTION -------------------------------------------
public:

private:
	static void cb_vs_fmot(Fl_Widget *wgt, void *idx);
};

#endif	// CPANEL
