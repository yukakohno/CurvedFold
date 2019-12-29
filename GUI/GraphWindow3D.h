
/*************************************************************
	GraphWindow3D.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_GL_Window.h>
#include "Quaternion.h"

#ifndef GRAPHWINDOW3D
#define GRAPHWINDOW3D

class GraphWindow3D : public Fl_Gl_Window
{
protected:
	void draw(); // to be overload
	int handle(int);
	void init();

	void SetLightMat();
	void SetRTS();
	void drawAxis();

	int n_plane;	// near plane
	int f_plane;	// far plane
	int w_w, w_h;	// window size

	int sts_alt, sts_sft, sts_ctrl;	// key status of Alt,Shift,Ctrl
	int push_button;	// button clicked
	int push_x, push_y;	// position clicked (used in dragging)

	CQuaternion quat;
	float prevquat[4], currquat[4], dispquat[4]; // before mouse drag, during mouse drag, for display
	float mCamera[16]; // camera matrix

public:
	int projection;
	float prevscale, currscale, dispscale;			// before mouse drag, during mouse drag, for display
	float prevtrans[4], currtrans[4], disptrans[4];	// before mouse drag, during mouse drag, for display

	void resetRotation(int up = 0); // up = 0:y-up, 1:z-up, 2:x-up
	void resetScale();
	void resetTranslation();

	GraphWindow3D(int X, int Y, int W, int H, const char *L);
	GraphWindow3D(int X, int Y, int W, int H);
};

#endif
