#define _USE_MATH_DEFINES
#include <math.h>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_vs_history(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	int hidx = (int)This->vs_history->value();
	if ( hidx<0 || MAX_HISTORY <= hidx ) {
		return;
	}
	int cpcnt, mode;
#if 0
	double tmpPx2d[MAX_CPCNT], tmpPy2d[MAX_CPCNT];
	double tmpPx[MAX_CPCNT], tmpPy[MAX_CPCNT], tmpPz[MAX_CPCNT];
	double tmpPa[MAX_CPCNT], tmpPbl[MAX_CPCNT], tmpPbr[MAX_CPCNT];
	int ret = access_hist(hidx, &cpcnt, &mode, tmpPx2d, tmpPy2d, tmpPx, tmpPy, tmpPz, tmpPa, tmpPbl, tmpPbr);
#else
	int ret = This->access_hist(hidx, &cpcnt, &mode, c->Px2d, c->Py2d, c->Px, c->Py, c->Pz, c->Pa, c->Pbl, c->Pbr, c->m3);
#endif
	if (ret < 0) {
		return;
	}

	// apply
	switch (mode) {
	case CMODE_A:	// A: 3D‹Èü‚ÆÜ‚èŠp“x‚©‚çAÜ‚èü‚ð‹‚ß‚é
		c->calcXA_CP(ppm->flg_interpolate, &ppm->rp);
		break;
	case CMODE_B:	// B: Ü‚èü‚ÆÜ‚èŠp“x‚©‚çA3D‹Èü‚ð‹‚ß‚é
		c->calcCPA_X(ppm->flg_interpolate, &ppm->rp);
		break;
	case CMODE_C:	// C: Ü‚èü‚Æ3D‹Èü‚©‚çAÜ‚èŠp“x‚ð‹‚ß‚é
		c->calcCPX_A(ppm->flg_interpolate, &ppm->rp);
		break;
	case CMODE_R:	// ruling angle -> torsion, alpha -> 3D‹ÈüAÜ‚èŠp“x
		c->calcR_TA(ppm->flg_interpolate, &ppm->rp, -1, -1, 2 * M_PI);
		break;
	}
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);
	//if (This->cb_optmat->value() && ppm->tgcnt > 3) {
	//	ppm->optMat(CMODE_R);
	//}

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

int ControlPanel::push_hist(int cpcnt, int mode, double* Px2d, double* Py2d,
	double* Px, double* Py, double* Pz, double* Pa, double* Pbl, double* Pbr, double* m3)
{
	bool flg_history = (bool)this->cb_history->value();
	if ( !flg_history ) {
		return -1;
	}
#if 0
	int head, size;
	int hist_mode[MAX_HISTORY], hist_cpcnt[MAX_HISTORY];
	double histPx2d[MAX_HISTORY][MAX_CPCNT], histPy2d[MAX_HISTORY][MAX_CPCNT];
	double histPx[MAX_HISTORY][MAX_CPCNT], histPy[MAX_HISTORY][MAX_CPCNT], histPz[MAX_HISTORY][MAX_CPCNT];
	double histPa[MAX_HISTORY][MAX_CPCNT];
	double histPbl[MAX_HISTORY][MAX_CPCNT], histPbr[MAX_HISTORY][MAX_CPCNT];
#endif
	if (hist_size < MAX_HISTORY ) {
		hist_mode[hist_size] = mode;
		hist_cpcnt[hist_size] = cpcnt;
		memcpy(histPx2d[hist_size], Px2d, sizeof(double)*cpcnt);
		memcpy(histPy2d[hist_size], Py2d, sizeof(double)* cpcnt);
		memcpy(histPx[hist_size], Px, sizeof(double)* cpcnt);
		memcpy(histPy[hist_size], Py, sizeof(double)* cpcnt);
		memcpy(histPz[hist_size], Pz, sizeof(double)* cpcnt);
		memcpy(histPa[hist_size], Pa, sizeof(double)* cpcnt);
		memcpy(histPbl[hist_size], Pbl, sizeof(double)* cpcnt);
		memcpy(histPbr[hist_size], Pbr, sizeof(double)* cpcnt);
		memcpy(hist_m3[hist_size], m3, sizeof(double) * 16);
		printf("push_hist: hist[%d] ", hist_size);
		hist_size++;
		hist_head = hist_size;
		if (hist_head == MAX_HISTORY) {
			hist_head = 0;
		}
	} else {
		hist_mode[hist_head] = mode;
		hist_cpcnt[hist_head] = cpcnt;
		memcpy(histPx2d[hist_head], Px2d, sizeof(double) * cpcnt);
		memcpy(histPy2d[hist_head], Py2d, sizeof(double) * cpcnt);
		memcpy(histPx[hist_head], Px, sizeof(double) * cpcnt);
		memcpy(histPy[hist_head], Py, sizeof(double) * cpcnt);
		memcpy(histPz[hist_head], Pz, sizeof(double) * cpcnt);
		memcpy(histPa[hist_head], Pa, sizeof(double) * cpcnt);
		memcpy(histPbl[hist_head], Pbl, sizeof(double) * cpcnt);
		memcpy(histPbr[hist_head], Pbr, sizeof(double) * cpcnt);
		memcpy(hist_m3[hist_head], m3, sizeof(double) * 16);
		printf("push_hist: hist[%d] ", hist_head);
		hist_head++;
		if (hist_head == MAX_HISTORY) {
			hist_head = 0;
		}
	}
	printf(" head=%d, size=%d\n", hist_head, hist_size);

	return 0;
}

int ControlPanel::access_hist(int idx, int* _cpcnt, int* _mode, double* Px2d, double* Py2d,
	double* Px, double* Py, double* Pz, double* Pa, double* Pbl, double* Pbr, double* m3)
{
	int idx0 = (hist_head-1) -idx;
	if (idx0<0) {
		idx0 = MAX_HISTORY + idx0;
	}
	if (idx0 >= hist_size) { return -1; }
	int cpcnt;
	*_mode = hist_mode[idx0];
	*_cpcnt = cpcnt = hist_cpcnt[idx0];
	memcpy(Px2d, histPx2d[idx0], sizeof(double) * cpcnt);
	memcpy(Py2d, histPy2d[idx0], sizeof(double) * cpcnt);
	memcpy(Px, histPx[idx0], sizeof(double) * cpcnt);
	memcpy(Py, histPy[idx0], sizeof(double) * cpcnt);
	memcpy(Pz, histPz[idx0], sizeof(double) * cpcnt);
	memcpy(Pa, histPa[idx0], sizeof(double) * cpcnt);
	memcpy(Pbl, histPbl[idx0], sizeof(double) * cpcnt);
	memcpy(Pbr, histPbr[idx0], sizeof(double) * cpcnt);
	memcpy(m3, hist_m3[idx0], sizeof(double) * 16);
	printf("access_hist: hist[%d], head=%d, size=%d\n", idx0, hist_head, hist_size);

	return 0;
}
