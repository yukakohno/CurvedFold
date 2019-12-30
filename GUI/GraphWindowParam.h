
/*************************************************************
	GraphWindowParam.h
	papermodel �̊e��p�����[�^�l���O���t�\��
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Double_Window.h>
#include "../CurvedFoldModel/papermodel.h"

#ifndef GRAPHWINDOWPARAM
#define GRAPHWINDOWPARAM

class GraphWindowParam : public Fl_Double_Window
{
protected:
	void draw();
	int handle(int);
	void init();

public:
	GraphWindowParam(int X, int Y, int W, int H, const char *L);
	GraphWindowParam(int X, int Y, int W, int H);

	papermodel *ppm;
	int cidx;
};

#endif
