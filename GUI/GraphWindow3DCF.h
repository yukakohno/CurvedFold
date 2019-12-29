
/*************************************************************
	GraphWindow3DCF.h
*************************************************************/

#include "GraphWindow3D.h"
#include "../CurvedFoldModel/papermodel.h"
#include "common.h"

#ifndef GRAPHWINDOWCF
#define GRAPHWINDOWCF

class GraphWindow3DCF : public GraphWindow3D
{
protected:
	void draw();
	int handle(int);

	void SetRTS2();
	void draw3DCurveFold();

public:
	int disp_axis;
	int disp_X;
	int disp_TNB;
	int disp_R;
	int disp_PLY;
	int disp_PTN;
	int disp_CP, ppos, pprm;

	int flg_addcurve; // 0:-, 1:trim, 2:fold

	papermodel *ppm;
	int disp_maxrlen;
	int sts_alt2, sts_sft2, sts_ctrl2;
	int push_button2;
	int push_x2, push_y2;

	CQuaternion quat2;
	float prevquat2[4], currquat2[4], dispquat2[4];
	float prevtrans2[4], currtrans2[4], disptrans2[4];
	float mObject[16];

	int vcnt;
	double *vx,*vy,*vz;

public:
	GraphWindow3DCF(int X, int Y, int W, int H, const char *L);
	GraphWindow3DCF(int X, int Y, int W, int H);
	~GraphWindow3DCF();

	void init();
	void initTexture();
	void initObject();
};

#endif	// GRAPHWINDOWCF
