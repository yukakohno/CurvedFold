#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <windows.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

// ------------------------- PARAM LIST --------------------------------------

void ControlPanel::cb_btn_makePrmList(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	This->pcnt = This->pidx = 0;
	int lastpcnt = 0;

	//
	// load from file
	//
	std::ifstream fin("input/paramlist.csv");
	if (fin) {
		std::string str_buf;
		std::string str_conma_buf;
		while (getline(fin, str_buf)) {

			if (str_buf[0] == '\n') {
				continue;
			}

			std::istringstream i_stream(str_buf);

			for (int j = 0; j < c->Pcnt; j++) {
				if (getline(i_stream, str_conma_buf, ',')) {
					This->prm[This->pcnt][j] = atof(str_conma_buf.c_str());
				}
			}
			for (int j = 0; j < c->Pcnt; j++) {
				if (getline(i_stream, str_conma_buf, ',')) {
					This->prm[This->pcnt][c->Pcnt + j] = atof(str_conma_buf.c_str());
				}
			}
			This->pcnt++;
		}
		fin.close();
		printf("parameters loaded from input/paramlist.csv. pcnt = %d\n", This->pcnt);
		return;
	}

	//
	// reset paramlist.csv
	//
	std::ofstream fout("./output/paramlist.csv");
	fout.close();

	//
	// create random params, check, file output
	//
	int pmax = 20;
	while (This->pcnt < pmax )
	{
		//int i0[10], j0[10];
		for (int i = 0; i < c->Pcnt; i++) {
			This->prm[This->pcnt][i] = ((double)rand() / (double)RAND_MAX - 0.5) * 0.1;// [-0.05, 0.05]
			This->prm[This->pcnt][c->Pcnt + i] = ((double)rand() / (double)RAND_MAX) * 0.5 * M_PI; // [0, PI/2]
		}

#if 1	// check parameter value
		int j = 0;
		for (j = 0; j < c->Pcnt; j++) {
			int jj = c->Pcnt + j;
			if (This->prm[This->pcnt][j] <= This->tmin || This->tmax <= This->prm[This->pcnt][j]
				|| This->prm[This->pcnt][jj] <= This->amin || This->amax <= This->prm[This->pcnt][jj])
			{
				break;
			}
		}
		if (j < c->Pcnt) {
			continue;
		}
#endif
		// check ruling crossing
		for (int j = 0; j < c->Pcnt; j++) {
			int jj = c->Pcnt + j;
			c->Py[j] = This->prm[This->pcnt][j]; // torsion
			c->Pa[j] = This->prm[This->pcnt][jj];// fold angle
		}
		c->calcCPA_X(1, &ppm->rp);
		c->Xsidx = CEMGN;
		c->Xeidx = c->Xcnt - CEMGN - 1;
		c->setEnds(ppm->psx, ppm->psy, ppm->pex, ppm->pey, 0, NULL);
		c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// Ž†’[‚Ü‚Å‚Ìruling’·‚³
		//c->setRLen( 30.0 );

		int res_rul0 = c->checkRulingAngle(); // 0:OK, -1: left ruling NG, 1: right ruling NG
		if (res_rul0 != 0) {
			continue;
		}
		int res_rul1 = c->checkRulingCross();
		if (res_rul1) {
			continue;
		}
		int res_rul2 = c->checkRulCreaseCross();
		if (res_rul2) {
			continue;
		}

		std::ofstream fout("./output/paramlist.csv", std::ios_base::app);
		for (int i = 0; i < c->Pcnt; i++) { fout << c->Py[i] << ","; }
		for (int i = 0; i < c->Pcnt; i++) { fout << c->Pa[i] << ","; }
		fout << std::endl;
		fout.close();

		std::cout << "This->pcnt = " << This->pcnt  << std::endl;
		for (int i = 0; i < c->Pcnt; i++) { std::cout << c->Py[i] << ","; }
		std::cout << std::endl;
		for (int i = 0; i < c->Pcnt; i++) { std::cout << c->Pa[i] << ","; }
		std::cout << std::endl;

		This->pcnt++;
	}
}

void ControlPanel::idleprm(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	Sleep(500);
	std::cout << This->pidx << "/" << This->pcnt << std::endl;

	if (This->pidx >= This->pcnt)
	{
		//Fl::remove_idle(ControlPanel::idle, This);
		//return;
		This->pidx = 0;
	}

	for (int i = 0; i < c->Pcnt; i++) {
		c->Py[i] = This->prm[This->pidx][i]; // torsion
		c->Pa[i] = This->prm[This->pidx][c->Pcnt + i];// fold angle
	}
	//std::cout << This->pidx << "/" << This->pcnt << std::endl;

	c->calcCPA_X(1, &ppm->rp);
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);

	if (This->pidx < This->listtgcnt) {
		ppm->tgcnt = This->tgcnt[This->pidx];
		memcpy(ppm->tgx, This->tgx[This->pidx], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->tgy, This->tgy[This->pidx], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->tgz, This->tgz[This->pidx], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->ogx_cp, This->ogx_cp[This->pidx], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->ogy_cp, This->ogy_cp[This->pidx], sizeof(double) * MAX_TGT_CNT);
		ppm->getTgt2D3D();
#if 0
		for (int i = 0; i < ppm->tgcnt; i++) {
			std::cout << i << ":" << ppm->tgx[i] << "," << ppm->tgy[i] << "," << ppm->tgz[i]
				<< "," << ppm->ogx[i] << "," << ppm->ogy[i] << "," << ppm->ogz[i]
				<< "," << ppm->ogx_cp[i] << "," << ppm->ogy_cp[i] << std::endl;
		}
		std::cout << std::endl;
#endif
		//This->btn_optfold2->do_callback();
	}

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	This->pidx++;
}

void ControlPanel::cb_btn_startprm(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (!This->flg_idle_activeprm) {
		This->flg_idle_activeprm = true;
		Fl::add_idle(ControlPanel::idleprm, This);
	}
}

void ControlPanel::cb_btn_stopprm(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (This->flg_idle_activeprm) {
		This->flg_idle_activeprm = false;
		Fl::remove_idle(ControlPanel::idleprm, This);
	}
}

// ------------------------- PARAM LIST2 --------------------------------------

void ControlPanel::cb_btn_makePrmList2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	This->pcnt2 = This->pidx2 = 0;

	//
	// load from file
	//
	std::ifstream fin("input/paramlist2.csv");
	if (fin) {
		std::string str_buf;
		std::string str_conma_buf;
		while (getline(fin, str_buf)) {

			if (str_buf[0] == '\n') {
				continue;
			}

			std::istringstream i_stream(str_buf);

			for (int j = 0; j < c->Pcnt; j++) {
				if (getline(i_stream, str_conma_buf, ',')) {
					This->prm2[This->pcnt2][j] = atof(str_conma_buf.c_str());
				}
			}
			for (int j = 0; j < c->Pcnt; j++) {
				if (getline(i_stream, str_conma_buf, ',')) {
					This->prm2[This->pcnt2][c->Pcnt + j] = atof(str_conma_buf.c_str());
				}
			}
			for (int j = 0; j < c->Pcnt; j++) {
				if (getline(i_stream, str_conma_buf, ',')) {
					This->prm2[This->pcnt2][c->Pcnt*2 + j] = atof(str_conma_buf.c_str());
				}
			}
			This->pcnt2++;
		}
		ppm->re_sidx = 0;
		ppm->re_eidx = c->Xcnt;
		fin.close();
		std::cout << "parameters loaded from input/paramlist.csv. pcnt2 = " << This->pcnt2 << std::endl;
		return;
	}

	//
	// reset paramlist.csv
	//
	std::ofstream fout("./output/paramlist2.csv");
	fout.close();

	//
	// create random params, check, file output
	//

	crease c0; memcpy(&c0, c, sizeof(crease));
	std::vector<double> tmpPbl;
	std::vector<double> tmpPbr;
	int tmppmax = 200;

	// random rulings left -> tmpPbl
	for (int i = 0; i < 1000000; i++)
	{
		for (int j = 0; j < c0.Pcnt; j++) {
			c0.Pbl[j] = ((double)rand() / (double)RAND_MAX) * M_PI; // [0, PI]
		}
		crease::interpolate_spline_RulAngle(c0.Pbl, c0.Pcnt, c0.betal, c0.Xcnt, c0.cotbl, c0.cosbl, c0.sinbl);
		for (int j = c0.Xsidx; j <= c0.Xeidx; j++) {
			c0.rlx_cp[j] = c0.cosbl[j] * c0.Tx2d[j] - c0.sinbl[j] * c0.Ty2d[j];
			c0.rly_cp[j] = c0.sinbl[j] * c0.Tx2d[j] + c0.cosbl[j] * c0.Ty2d[j];
			normalize_v2(&(c0.rlx_cp[j]), &(c0.rly_cp[j]));
			c0.rllen[j] = ppm->pw * ppm->ph;
		}
		c0.calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj
		if (c0.checkRulingCross(-1) > 0) { // -1:left, 1:right
			continue;
		}
		for (int k = 0; k < c0.Pcnt; k++) {
			tmpPbl.push_back(c0.Pbl[k]);
		}
		if (tmpPbl.size() >= tmppmax * c0.Pcnt) {
			break;
		}
	}

	// random rulings right -> tmpPbr
	for (int i = 0; i < 1000000; i++)
	{
		for (int j = 0; j < c0.Pcnt; j++) {
			c0.Pbr[j] = ((double)rand() / (double)RAND_MAX) * M_PI; // [0, PI]
		}
		crease::interpolate_spline_RulAngle(c0.Pbr, c0.Pcnt, c0.betar, c0.Xcnt, c0.cotbr, c0.cosbr, c0.sinbr);
		for (int j = c0.Xsidx; j <= c0.Xeidx; j++) {
			c0.rrx_cp[j] = c0.cosbr[j] * c0.Tx2d[j] + c0.sinbr[j] * c0.Ty2d[j];
			c0.rry_cp[j] = -c0.sinbr[j] * c0.Tx2d[j] + c0.cosbr[j] * c0.Ty2d[j];
			normalize_v2(&(c0.rrx_cp[j]), &(c0.rry_cp[j]));
			c0.rrlen[j] = ppm->pw * ppm->ph;
		}
		c0.calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj
		if (c0.checkRulingCross(1) > 0) { // -1:left, 1:right
			continue;
		}
		for (int k = 0; k < c0.Pcnt; k++) {
			tmpPbr.push_back(c0.Pbr[k]);
		}
		if (tmpPbr.size() >= tmppmax * c0.Pcnt) {
			break;
		}
	}

	int cntl = (tmpPbl.size() / c0.Pcnt);
	int cntr = (tmpPbr.size() / c0.Pcnt);
	int cnt = cntl < cntr ? cntl : cntr;
	int pmax2 = 100;

	std::cout << "random rulings: left = " << cntl << ", right =" << cntr << std::endl;

	// This->vPbl + This->vPbl + random folding angle -> This->tgt*
	for (int j = 0; j < cnt; j++) {

		for (int i = 0; i < c0.Pcnt; i++) {
			c0.Pbl[i] = tmpPbl[j * c0.Pcnt + i];
		}
		for (int i = 0; i < c0.Pcnt; i++) {
			c0.Pbr[i] = tmpPbr[j * c0.Pcnt + i];
		}
		bool flg = false;
		for (int i = 0; i < 1000000; i++)
		{
			double mang = ((double)rand() / (double)RAND_MAX - 0.5) * M_PI; // [-PI/2, PI/2]
			ppm->re_sidx = 0;
			ppm->re_eidx = c0.Xcnt;
			int ret = c0.calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, mang, 0);
			if (ret == 0) {
				flg = true;
				break;
			}
		}
		if (flg) {
			This->tgPcnt[This->pcnt2] = c0.Pcnt;
			memcpy(&(This->prm2[This->pcnt2][0]), c0.Pbl, sizeof(double) * MAX_CPCNT);
			memcpy(&(This->prm2[This->pcnt2][c0.Pcnt]), c0.Pbr, sizeof(double) * MAX_CPCNT);
			memcpy(&(This->prm2[This->pcnt2][c0.Pcnt*2]), c0.Pa, sizeof(double) * MAX_CPCNT);

			int ret = c0.calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, c0.Pa[c0.Pcnt / 2], 0);
			c0.Xsidx = CEMGN;
			c0.Xeidx = c0.Xcnt - CEMGN - 1;
			c0.setEnds(ppm->psx, ppm->psy, ppm->pex, ppm->pey, 0, NULL);
			c0.calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// Ž†’[‚Ü‚Å‚Ìruling’·‚³

			int res_rul0 = c0.checkRulingAngle(); // 0:OK, -1: left ruling NG, 1: right ruling NG
			if (res_rul0 != 0) {
				continue;
			}
			int res_rul1 = c0.checkRulingCross();
			if (res_rul1) {
				continue;
			}
			int res_rul2 = c0.checkRulCreaseCross();
			if (res_rul2) {
				continue;
			}

			std::ofstream fout("./output/paramlist2.csv", std::ios_base::app);
			for (int i = 0; i < c0.Pcnt; i++) { fout << c0.Pbl[i] << ","; }
			for (int i = 0; i < c0.Pcnt; i++) { fout << c0.Pbr[i] << ","; }
			for (int i = 0; i < c0.Pcnt; i++) { fout << c0.Pa[i] << ","; }
			fout << std::endl;
			fout.close();

			std::cout << "This->pcnt2 = " << This->pcnt2 << std::endl;
			for (int i = 0; i < c0.Pcnt; i++) { std::cout << c0.Pbl[i] << ","; }
			std::cout << std::endl;
			for (int i = 0; i < c0.Pcnt; i++) { std::cout << c0.Pbr[i] << ","; }
			std::cout << std::endl;
			for (int i = 0; i < c0.Pcnt; i++) { std::cout << c0.Pa[i] << ","; }
			std::cout << std::endl;

			This->pcnt2++;
			if (This->pcnt2 >= pmax2) {
				break;
			}
		}
	}
}

void ControlPanel::idleprm2(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	Sleep(500);
	std::cout << This->pidx2 << "/" << This->pcnt2 << std::endl;

	if (This->pidx2 >= This->pcnt2)
	{
		//Fl::remove_idle(ControlPanel::idle2, This);
		//return;
		This->pidx2 = 0;
	}

	memcpy(c->Pbl, &(This->prm2[This->pidx2][0]), sizeof(double) * c->Pcnt);
	memcpy(c->Pbr, &(This->prm2[This->pidx2][c->Pcnt]), sizeof(double) * c->Pcnt);
	memcpy(c->Pa, &(This->prm2[This->pidx2][c->Pcnt*2]), sizeof(double) * c->Pcnt);

	int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, c->Pa[c->Pcnt / 2], 0);
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);

	if(This->pidx2 < This->listtgcnt){
		ppm->tgcnt = This->tgcnt[This->pidx2];
		memcpy(ppm->tgx, This->tgx[This->pidx2], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->tgy, This->tgy[This->pidx2], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->tgz, This->tgz[This->pidx2], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->ogx_cp, This->ogx_cp[This->pidx2], sizeof(double) * MAX_TGT_CNT);
		memcpy(ppm->ogy_cp, This->ogy_cp[This->pidx2], sizeof(double) * MAX_TGT_CNT);
		ppm->getTgt2D3D();
#if 0
		for (int i = 0; i < ppm->tgcnt; i++) {
			std::cout << i << ":" << ppm->tgx[i] << "," << ppm->tgy[i] << "," << ppm->tgz[i]
				<< "," << ppm->ogx[i] << "," << ppm->ogy[i] << "," << ppm->ogz[i]
				<< "," << ppm->ogx_cp[i] << "," << ppm->ogy_cp[i] << std::endl;
		}
		std::cout << std::endl;
#endif
		//This->btn_optfold2->do_callback();
	}
	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	This->pidx2++;
}

void ControlPanel::cb_btn_startprm2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (!This->flg_idle_activeprm2) {
		This->flg_idle_activeprm2 = true;
		Fl::add_idle(ControlPanel::idleprm2, This);
	}
}

void ControlPanel::cb_btn_stopprm2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (This->flg_idle_activeprm2) {
		This->flg_idle_activeprm2 = false;
		Fl::remove_idle(ControlPanel::idleprm2, This);
	}
}

// ------------------------- SAMPLE TARGET POINTS --------------------------------------

void ControlPanel::cb_btn_listtgt(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	crease c0; memcpy(&c0, c, sizeof(crease));

	This->listtgcnt = 0;
	for (int j = 0; j < This->pcnt; j++) {
		for (int i = 0; i < c->Pcnt; i++) {
			c->Py[i] = This->prm[j][i]; // torsion
			c->Pa[i] = This->prm[j][c->Pcnt + i];// fold angle
		}
		c->calcCPA_X(1, &ppm->rp);
		ppm->set_postproc_type(PPTYPE_PRICURVE);
		ppm->postproc();
		ppm->set_postproc_type(PPTYPE_UNDEF);
		{
			std::stringstream ss; ss << "./output/P" << std::setw(3) << std::setfill('0') << This->listtgcnt << ".txt";
			std::cout << "export " << ss.str().c_str() << std::endl;
			c->dumpP((char*)ss.str().c_str(), CMODE_B);
		}

		This->tgPcnt[This->listtgcnt] = c->Pcnt;
		This->tgcnt[This->listtgcnt] = ppm->tgcnt;
		memcpy(This->tgx[This->listtgcnt], ppm->ogx, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->tgy[This->listtgcnt], ppm->ogy, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->tgz[This->listtgcnt], ppm->ogz, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->ogx_cp[This->listtgcnt], ppm->ogx_cp, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->ogy_cp[This->listtgcnt], ppm->ogy_cp, sizeof(double) * MAX_TGT_CNT);
		{
			int i = This->listtgcnt;
			std::stringstream ss; ss << "./output/target" << std::setw(3) << std::setfill('0') << i << ".txt";
			std::cout << "export " << ss.str() << std::endl;
			std::ofstream ofs(ss.str());
			for (int j = 0; j < This->tgcnt[i]; j++) {
				ofs << This->tgx[i][j] << "\t" << This->tgy[i][j] << "\t" << This->tgz[i][j]
					<< "\t" << This->ogx_cp[i][j] << "\t" << This->ogy_cp[i][j] << std::endl;
			}
			ofs.close();
		}

		This->listtgcnt++;
		std::cout << "listtgcnt=" << This->listtgcnt << std::endl;

		if (This->listtgcnt >= MAX_TGTLST) {
			break;
		}
	}

	memcpy(c, &c0, sizeof(crease));
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);
}

void ControlPanel::cb_btn_listtgt2(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);
	crease c0; memcpy(&c0, c, sizeof(crease));

	This->listtgcnt = 0;
	for (int j = 0; j < This->pcnt2; j++) {

		memcpy(c->Pbl, &(This->prm2[j][0]), sizeof(double) * c->Pcnt);
		memcpy(c->Pbr, &(This->prm2[j][c->Pcnt]), sizeof(double) * c->Pcnt);
		memcpy(c->Pa, &(This->prm2[j][c->Pcnt * 2]), sizeof(double) * c->Pcnt);

		int ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, c->Pa[c->Pcnt / 2], 0);
		ppm->set_postproc_type(PPTYPE_PRICURVE);
		ppm->postproc();
		ppm->set_postproc_type(PPTYPE_UNDEF);
		{
			std::stringstream ss; ss << "./output/rulings" << std::setw(3) << std::setfill('0') << This->listtgcnt << ".txt";
			std::ofstream ofs(ss.str());
			for (int j = 0; j < c->Pcnt; j++) { ofs << c->Pbl[j] << "\t"; }
			ofs << std::endl;
			for (int j = 0; j < c->Pcnt; j++) { ofs << c->Pbr[j] << "\t"; }
			ofs << std::endl;
			ofs.close();
		}

		This->tgPcnt[This->listtgcnt] = c->Pcnt;
		This->tgcnt[This->listtgcnt] = ppm->tgcnt;
		memcpy(This->tgx[This->listtgcnt], ppm->ogx, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->tgy[This->listtgcnt], ppm->ogy, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->tgz[This->listtgcnt], ppm->ogz, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->ogx_cp[This->listtgcnt], ppm->ogx_cp, sizeof(double) * MAX_TGT_CNT);
		memcpy(This->ogy_cp[This->listtgcnt], ppm->ogy_cp, sizeof(double) * MAX_TGT_CNT);
		{
			int i = This->listtgcnt;
			std::stringstream ss; ss << "./output/target" << std::setw(3) << std::setfill('0') << i << ".txt";
			std::cout << "export " << ss.str() << std::endl;
			std::ofstream ofs(ss.str());
			for (int j = 0; j < This->tgcnt[i]; j++) {
				ofs << This->tgx[i][j] << "\t" << This->tgy[i][j] << "\t" << This->tgz[i][j]
					<< "\t" << This->ogx_cp[i][j] << "\t" << This->ogy_cp[i][j] << std::endl;
			}
			ofs.close();
		}

		This->listtgcnt++;
		std::cout << "listtgcnt=" << This->listtgcnt << std::endl;

		if (This->listtgcnt >= MAX_TGTLST) {
			break;
		}
	}

	memcpy(c, &c0, sizeof(crease));
	ppm->set_postproc_type(PPTYPE_PRICURVE);
	ppm->postproc();
	ppm->set_postproc_type(PPTYPE_UNDEF);
}

void ControlPanel::cb_btn_listtgtmasks(Fl_Widget* wgt, void* idx)
{
#define TGTSIZE 9
	int tgtnum = 0;
	int tgt[TGTSIZE * TGTSIZE] = { 0 };

	std::cout << "input number of target points" << std::endl;
	std::cin >> tgtnum;

	int tgtcnt = 0;
	while (tgtcnt < tgtnum) {
		int idx = rand() % (TGTSIZE * TGTSIZE);
		if (tgt[idx] == 0) {
			tgt[idx] = 1;
			tgtcnt++;
		}
	}

	std::stringstream ss; ss << "./output/tmask" << std::setw(2) << std::setfill('0') << tgtnum << ".txt";
	std::cout << "export " << ss.str() << std::endl;
	std::ofstream ofs(ss.str());
	if(ofs){
		for (int j = 0; j < TGTSIZE; j++) {
			for (int i = 0; i < TGTSIZE; i++) {
				ofs << tgt[i + j * TGTSIZE] << " ";
			}
			ofs << std::endl;
		}
		ofs.close();
	}
}

#if 0
// ------------------------- RULINGS LIST --------------------------------------

#define IPP 10
int checkRulCross_recursive(papermodel* ppm, crease* c, int pidx, int rl, int* cnt, std::vector<double>& vPb)
{
	int idx_crossing = -1;

	if (pidx == c->Pcnt - 1) {

		// TODO: check starting ruling angle
		int i_start = 1;

		for (int i = i_start; i < 180; i += IPP) {
			if (rl < 0) {
				c->Pbl[pidx] = (double)i / 180.0 * M_PI;
				crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
				for (int j = c->Xsidx; j <= c->Xeidx; j++) {
					//c->betal[j] = atan(1.0 / c->cotbl[j]);
					//if (c->betal[j] < 0.0) {
					//	c->betal[j] += M_PI;
					//}
					c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
					c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
					normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
					c->rllen[j] = ppm->pw * ppm->ph;
				}
			}
			else {
				c->Pbr[pidx] = (double)i / 180.0 * M_PI;
				crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
				for (int j = c->Xsidx; j <= c->Xeidx; j++) {
					//c->betar[j] = atan(1.0 / c->cotbr[j]);
					//if (c->betar[j] < 0.0) {
					//	c->betar[j] += M_PI;
					//}
					c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
					c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
					normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
					c->rrlen[j] = ppm->pw * ppm->ph;
				}
			}
			c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj

			// return index of rulings crossing, -1 if no crossing
			idx_crossing = c->checkRulingCross(rl); // -1:left, 1:right
			if (idx_crossing > -1) {
#if 0
				std::cout << "idx_crossing=" << idx_crossing;
				if (rl < 0) {
					for (int k = 0; k < c->Pcnt; k++) { std::cout << ", " << c->Pbl[k]; }
				}
				else {
					for (int k = 0; k < c->Pcnt; k++) { std::cout << ", " << c->Pbr[k]; }
				}
				std::cout << std::endl;
#endif
				break;
			}

			(*cnt)++;
#if 1
			{
				//std::cout << "OK, cnt="<< (*cnt);
				if (rl < 0) {
					for (int k = 0; k < c->Pcnt; k++) {
						//std::cout << ", " << c->Pbl[k];
						vPb.push_back(c->Pbl[k]);
					}
				}
				else {
					for (int k = 0; k < c->Pcnt; k++) {
						//std::cout << ", " << c->Pbr[k];
						vPb.push_back(c->Pbr[k]);
					}
				}
				//std::cout << std::endl;
			}
#endif
		}
	}
	else {
		double* Pb = NULL;
		if (rl < 0) {
			Pb = c->Pbl;
		}
		else {
			Pb = c->Pbr;
		}

		// TODO: check starting ruling angle
		int i_start = 1;

		for (int i = i_start; i < 180; i += IPP) {
			Pb[pidx] = (double)i / 180.0 * M_PI;

			idx_crossing = checkRulCross_recursive(ppm, c, pidx + 1, rl, cnt, vPb);
			if (idx_crossing > -1) {
				continue;
			}
		}
	}
	return idx_crossing;
}

void ControlPanel::cb_btn_makeRulList(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* _c = &(ppm->crs[0]);
	crease c; memcpy(&c, _c, sizeof(crease));
	int idx_crossing = -1, lcnt = 0, rcnt = 0;

	if (This->rb_listleft->value() != 0) {
		This->vPbl.clear();
		idx_crossing = checkRulCross_recursive(ppm, &c, 0, -1, &lcnt, This->vPbl);
		std::cout << "vPbl->size = " << This->vPbl.size() << ", " << This->vPbl.size() / c.Pcnt << std::endl;
	}
	else {
		This->vPbr.clear();
		idx_crossing = checkRulCross_recursive(ppm, &c, 0, 1, &rcnt, This->vPbr);
		std::cout << "vPbr->size = " << This->vPbr.size() << ", " << This->vPbr.size() / c.Pcnt << std::endl;
	}
}

void ControlPanel::idlerul(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if (This->rb_listleft->value() != 0) {

		if (This->idlerul_idx >= This->vPbl.size() / c->Pcnt) {
			This->idlerul_idx = 0;
		}
#if 0
		std::cout << "idlerul_idx = " << This->idlerul_idx;
		for (int i = 0; i < c->Pcnt; i++) {
			std::cout << ", " << This->vPbl[This->idlerul_idx * c->Pcnt + i];
		}
		std::cout << std::endl;
#endif
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbl[i] = This->vPbl[This->idlerul_idx * c->Pcnt + i];
		}
		crease::interpolate_spline_RulAngle(c->Pbl, c->Pcnt, c->betal, c->Xcnt, c->cotbl, c->cosbl, c->sinbl);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rlx_cp[j] = c->cosbl[j] * c->Tx2d[j] - c->sinbl[j] * c->Ty2d[j];
			c->rly_cp[j] = c->sinbl[j] * c->Tx2d[j] + c->cosbl[j] * c->Ty2d[j];
			normalize_v2(&(c->rlx_cp[j]), &(c->rly_cp[j]));
			c->rllen[j] = ppm->pw * ppm->ph;
		}
	}
	else {
		if (This->idlerul_idx >= This->vPbr.size() / c->Pcnt) {
			This->idlerul_idx = 0;
		}
#if 0
		std::cout << "idlerul_idx = " << This->idlerul_idx;
		for (int i = 0; i < c->Pcnt; i++) {
			std::cout << ", " << This->vPbr[This->idlerul_idx * c->Pcnt + i];
		}
		std::cout << std::endl;
#endif
		for (int i = 0; i < c->Pcnt; i++) {
			c->Pbr[i] = This->vPbr[This->idlerul_idx * c->Pcnt + i];
		}
		crease::interpolate_spline_RulAngle(c->Pbr, c->Pcnt, c->betar, c->Xcnt, c->cotbr, c->cosbr, c->sinbr);
		for (int j = c->Xsidx; j <= c->Xeidx; j++) {
			c->rrx_cp[j] = c->cosbr[j] * c->Tx2d[j] + c->sinbr[j] * c->Ty2d[j];
			c->rry_cp[j] = -c->sinbr[j] * c->Tx2d[j] + c->cosbr[j] * c->Ty2d[j];
			normalize_v2(&(c->rrx_cp[j]), &(c->rry_cp[j]));
			c->rrlen[j] = ppm->pw * ppm->ph;
		}
	}
	c->calcRLenP(ppm->psx, ppm->psy, ppm->pex, ppm->pey);		// ruling’·‚³C³i˜güA‹Èü‚Ü‚Å‚Ì’·‚³—Dæj

	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
	Sleep(1);

	This->idlerul_idx++;
}
void ControlPanel::cb_btn_startrul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	//This->acnt = (int)This->vs_xmang1->value();
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if (!This->flg_idlerul_active) {
		if (This->rb_listleft->value() != 0 && This->vPbl.size() == 0) {
			return;
		}
		if (This->rb_listright->value() != 0 && This->vPbr.size() == 0) {
			return;
		}
		This->flg_idlerul_active = true;
		Fl::add_idle(ControlPanel::idlerul, This);
	}
}

void ControlPanel::cb_btn_stoprul(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	if (This->flg_idlerul_active) {
		This->flg_idlerul_active = false;
		Fl::remove_idle(ControlPanel::idlerul, This);
	}
}

#endif
