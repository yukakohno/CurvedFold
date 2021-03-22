#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

// 15 initial shapes
// 15 target shapes
// 10 taget point masks
std::string fname_init[15] = { "./input/param_synth/rulings00.txt",
	"./input/param_synth/rulings01.txt",
	"./input/param_synth/rulings04.txt",
	"./input/param_synth/rulings05.txt",
	"./input/param_synth/rulings06.txt",
	"./input/param_rulings/rulings000.txt",
	"./input/param_rulings/rulings001.txt",
	"./input/param_rulings/rulings002.txt",
	"./input/param_rulings/rulings003.txt",
	"./input/param_rulings/rulings004.txt",
	"./input/param_tr_ang/P001.txt",
	"./input/param_tr_ang/P002.txt",
	"./input/param_tr_ang/P003.txt",
	"./input/param_tr_ang/P004.txt",
	"./input/param_tr_ang/P005.txt" };
std::string fname_target[15] = { "./input/param_synth/target00.txt",
	"./input/param_synth/target01.txt",
	"./input/param_synth/target04.txt",
	"./input/param_synth/target05.txt",
	"./input/param_synth/target06.txt",
	"./input/param_rulings/target000.txt",
	"./input/param_rulings/target001.txt",
	"./input/param_rulings/target002.txt",
	"./input/param_rulings/target003.txt",
	"./input/param_rulings/target004.txt",
	"./input/param_tr_ang/target000.txt",
	"./input/param_tr_ang/target001.txt",
	"./input/param_tr_ang/target002.txt",
	"./input/param_tr_ang/target003.txt",
	"./input/param_tr_ang/target004.txt" };
std::string fname_tmask[10] = { "./input/tmasks/tmask81.txt",
	"./input/tmasks/tmask_left.txt",
	"./input/tmasks/tmask_convex.txt",
	"./input/tmasks/tmask_concave.txt",
	"./input/tmasks/tmask40_even.txt",
	"./input/tmasks/tmask10_even.txt",
	"./input/tmasks/tmask05_even.txt",
	"./input/tmasks/tmask40_random.txt",
	"./input/tmasks/tmask10_random.txt",
	"./input/tmasks/tmask05_random.txt" };

// results
#define INPUT_DATA_SIZE 15
#define TARGET_DATA_SIZE 15
#define TMASK_DATA_SIZE 10
#define TRIAL_SIZE (INPUT_DATA_SIZE*TARGET_DATA_SIZE*TMASK_DATA_SIZE)
double proc_time;
int batch_i = -1, batch_j = -1, batch_k = -1, batch_phase = -1;

#define FILENAME_RESULT "./output/result.csv"
#define FILENAME_PROCESS "./output/process.csv"
#define FILENAME_PARAM "./output/param.csv"

void ControlPanel::idle_batchproc(void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	if ((This->cnt_batchproc == 0 && This->phase_batchproc == 0) || batch_i == batch_j || batch_phase == 3) {
	} else{
		char fname[128];
		sprintf(fname, "./output/input%02d_target%02d_mask%02d_phase%d_model.bmp", batch_i, batch_j, batch_k, batch_phase);
		This->gwin->exportImage(fname);
		sprintf(fname, "./output/input%02d_target%02d_mask%02d_phase%d_crease.bmp", batch_i, batch_j, batch_k, batch_phase);
		This->gwin_cp->exportImage(fname);
		sprintf(fname, "./output/input%02d_target%02d_mask%02d_phase%d_creasecrop.bmp", batch_i, batch_j, batch_k, batch_phase);
		This->gwin_cp->exportCroppedImage(fname);
	}

	int i = batch_i = This->cnt_batchproc / (TARGET_DATA_SIZE * TMASK_DATA_SIZE);
	int j = batch_j = (This->cnt_batchproc / TMASK_DATA_SIZE) % TARGET_DATA_SIZE;
	int k = batch_k = This->cnt_batchproc % TMASK_DATA_SIZE;
	int phase = batch_phase = This->phase_batchproc;

	This->in_input_no->value(std::to_string(i).c_str());
	This->in_target_no->value(std::to_string(j).c_str());
	This->in_mask_no->value(std::to_string(k).c_str());

	if (This->cnt_batchproc >= TRIAL_SIZE) {
		std::cout << "process index = " << i << ", " << j << ", " << k << ", phase = " << phase << ", remove_idle()" << std::endl;
		Fl::remove_idle(ControlPanel::idle_batchproc, This);
		return;
	}

	std::cout << "process index = " << i << ", " << j << ", " << k << ", phase = " << phase << std::endl;

	if (This->phase_batchproc == 3) {
		This->phase_batchproc = 0;
		This->cnt_batchproc++;
	}
	else {
		This->phase_batchproc++;
	}

	if (i == j) {
		return;
	}

	switch (phase)
	{
	case 0:
	{
		//
		// load initial shape
		//
		ppm->tgcnt = 0;

		if (i < 10) { // load rulings
			c->load("input/P0.txt");
			c->loadm2m3("input/P0_m2m3.txt");
			c->loadPb((char*)fname_init[i].c_str());
			ppm->re_sidx = 0;
			ppm->re_eidx = c->Xcnt;
			int ret = 0;
			if (i < 5) { // Synthetic data
				ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, 30.0 / 180.0 * M_PI, 0);
			}
			else { // random data
				ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, This->prm2[i - 5][c->Pcnt * 2 + c->Pcnt / 2], 0);
			}
			if (ret != 0) {
				std::cout << "ERROR (" << i << ", " << j << ", " << k << ")" << std::endl;
			}
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
		}
		else { // load random data (torsion and folding angle)
			c->load((char*)fname_init[i].c_str());
			c->loadm2m3("input/P0_m2m3.txt");
			c->calcCPA_X(ppm->flg_interpolate, &ppm->rp);
			ppm->set_postproc_type(PPTYPE_PRICURVE);
			ppm->postproc();
			ppm->set_postproc_type(PPTYPE_UNDEF);
		}
	}
	break;
	case 1:
	{
		//
		// load target points & masks
		//
		ppm->loadTgtMask((char*)fname_target[j].c_str(), (char*)fname_tmask[k].c_str());
	}
	break;
	case 2:
	{
		//
		// optimize
		//
		clock_t start_clock, end_clock;
		start_clock = clock();

		This->btn_randrul2->do_callback();

		end_clock = clock();
		proc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
		//std::cout << "process time: " << proc_time << std::endl;
	}
	break;
	case 3:
	{
		//
		// evaluate
		//
		ppm->loadTgtMask((char*)fname_target[j].c_str(), (char*)fname_tmask[0].c_str());
		ppm->calcAvetgapMat(); // calculate ppm->avetgap;

		//
		// save result
		//
		{
			time_t t0 = time(NULL);
			struct tm* t1 = localtime(&t0);

			std::ofstream ofs(FILENAME_RESULT, std::ios_base::app);
			ofs << i << "," << j << "," << k << ","
				<< This->optlog_itr << "," << This->optlog_minerr[This->optlog_itr - 1] << "," << ppm->avetgap << ","
				<< proc_time << ","
				<< std::setw(2) << std::setfill('0') << t1->tm_hour
				<< std::setw(2) << std::setfill('0') << t1->tm_min
				<< std::setw(2) << std::setfill('0') << t1->tm_sec << std::endl;
			ofs.close();
		}
		{
			std::ofstream ofs(FILENAME_PROCESS, std::ios_base::app);
			ofs << i << "," << j << "," << k << ",err,";
			for (int i = 0; i < This->optlog_itr; i++)
			{
				ofs << This->optlog_err[i] << ",";
			}
			ofs << std::endl;
			ofs << i << "," << j << "," << k << ",minerr,";
			for (int i = 0; i < This->optlog_itr; i++)
			{
				ofs << This->optlog_minerr[i] << ",";
			}
			ofs << std::endl;
			ofs.close();

			This->optlog_itr = -1;
		}
		{
			std::ofstream ofs(FILENAME_PARAM, std::ios_base::app);
			ofs << i << "," << j << "," << k << ",,";
			for (int i = 0; i < c->Pcnt; i++)
			{
				ofs << c->Pbl[i] << ",";
			}
			ofs << ",";
			for (int i = 0; i < c->Pcnt; i++)
			{
				ofs << c->Pbr[i] << ",";
			}
			ofs << ",";
			for (int i = 0; i < c->Pcnt; i++)
			{
				ofs << c->Pa[i] << ",";
			}
			ofs << ",";
			for (int i = 0; i < c->Pcnt; i++)
			{
				ofs << c->Px[i] << ",";
			}
			ofs << ",";
			for (int i = 0; i < c->Pcnt; i++)
			{
				ofs << c->Py[i] << ",";
			}
			ofs << std::endl;
			ofs.close();
		}
		This->optlog_itr = -1;
	}
	break;
	}

	This->gwin->redraw();
	This->gwin_cp->redraw();
	This->gwin_gr->redraw();
}

void ControlPanel::cb_btn_batchproc(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	papermodel* ppm = &(This->ppm);
	crease* c = &(ppm->crs[0]);

	This->cb_optmat->value(1);
	This->btn_makePrmList2->do_callback();

#if 0
	int input_no = atoi(This->in_input_no->value());
	int target_no = atoi(This->in_target_no->value());
	int mask_no = atoi(This->in_mask_no->value());

	This->cnt_batchproc = (input_no*TARGET_DATA_SIZE + target_no)* TMASK_DATA_SIZE + mask_no;
	This->phase_batchproc = 0;

	// reset log file
	if (This->cnt_batchproc == 0) {
		std::ofstream ofs(FILENAME_RESULT);
		ofs << "initial,target,mask,iteration,minerr,ave dist,proc time,time" << std::endl;
		ofs.close();
		ofs.open(FILENAME_PROCESS);
		ofs << "initial,target,mask," << std::endl;
		ofs.close();
		ofs.open(FILENAME_PARAM);
		ofs << "initial,target,mask,,Pbl,,,,,,,,Pbr,,,,,,,,Pa,,,,,,,,Pk,,,,,,,,Pt,,,,,,,," << std::endl;
		ofs.close();
	}
	Fl::add_idle(ControlPanel::idle_batchproc, This);

#else

	std::ofstream ofs(FILENAME_RESULT);
	ofs << "initial,target,mask,iteration,minerr,ave dist,proc time,time" << std::endl;
	ofs.close();
	ofs.open(FILENAME_PROCESS);
	ofs << "initial,target,mask," << std::endl;
	ofs.close();
	ofs.open(FILENAME_PARAM);
	ofs << "initial,target,mask,,Pbl,,,,,,,,Pbr,,,,,,,,Pa,,,,,,,,Pk,,,,,,,,Pt,,,,,,,," << std::endl;
	ofs.close();

	for (int i = 0; i < INPUT_DATA_SIZE; i++) {
		for (int j = 0; j < TARGET_DATA_SIZE; j++) {
			if (i == j) {
				continue;
			}

			for (int k = 0; k < TMASK_DATA_SIZE; k++) {
				std::cout << "process index = " << i << ", " << j << ", " << k << std::endl;

				//
				// load initial shape
				//
				if (i < 10) { // load rulings
					c->load("input/P0.txt");
					c->loadm2m3("input/P0_m2m3.txt");
					c->loadPb((char *)fname_init[i].c_str());
					ppm->re_sidx = 0;
					ppm->re_eidx = c->Xcnt;
					int ret = 0;
					if (i < 5) { // Synthetic data
						ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, 30.0 / 180.0 * M_PI, 0);
					}
					else { // random data
						ret = c->calcR_TA(1/*flg_interpolate*/, &ppm->rp, ppm->re_sidx, ppm->re_eidx, This->prm2[i - 5][c->Pcnt * 2 + c->Pcnt / 2], 0);
					}
					if (ret != 0) {
						std::cout << "ERROR (" << i << ", " << j << ", " << k << ")" << std::endl;
					}
					ppm->set_postproc_type(PPTYPE_PRICURVE);
					ppm->postproc();
					ppm->set_postproc_type(PPTYPE_UNDEF);
				}
				else{ // load random data (torsion and folding angle)
					c->load((char*)fname_init[i].c_str());
					c->loadm2m3("input/P0_m2m3.txt");
					c->calcCPA_X(ppm->flg_interpolate, &ppm->rp);
					ppm->set_postproc_type(PPTYPE_PRICURVE);
					ppm->postproc();
					ppm->set_postproc_type(PPTYPE_UNDEF);
				}

				//
				// load target points & masks
				//
				ppm->loadTgtMask((char*)fname_target[j].c_str(), (char*)fname_tmask[k].c_str());

				//
				// optimize
				//
				clock_t start_clock, end_clock;
				start_clock = clock();

				This->btn_randrul2->do_callback();

				end_clock = clock();
				//std::cout << "process time: " << (double)(end_clock - start_clock) / CLOCKS_PER_SEC << std::endl;

				//
				// evaluate
				//
				ppm->loadTgtMask((char*)fname_target[j].c_str(), (char*)fname_tmask[0].c_str());
				ppm->calcAvetgapMat(); // calculate ppm->avetgap;

				//
				// save result
				//
				{
					time_t t0 = time(NULL);
					struct tm* t1 = localtime(&t0);

					std::ofstream ofs(FILENAME_RESULT, std::ios_base::app);
					ofs << i << "," << j << "," << k << ","
						<< This->optlog_itr << "," << This->optlog_minerr[This->optlog_itr-1] << "," << ppm->avetgap << ","
						<< (double)(end_clock - start_clock) / CLOCKS_PER_SEC << "," 
						<< std::setw(2) << std::setfill('0') << t1->tm_hour
						<< std::setw(2) << std::setfill('0') << t1->tm_min
						<< std::setw(2) << std::setfill('0') << t1->tm_sec << std::endl;
					ofs.close();
				}
				{
					std::ofstream ofs(FILENAME_PROCESS, std::ios_base::app);
					ofs << i << "," << j << "," << k << ",err,";
					for (int i = 0; i < This->optlog_itr; i++)
					{
						ofs << This->optlog_err[i] << ",";
					}
					ofs << std::endl;
					ofs << i << "," << j << "," << k << ",minerr,";
					for (int i = 0; i < This->optlog_itr; i++)
					{
						ofs << This->optlog_minerr[i] << ",";
					}
					ofs << std::endl;
					ofs.close();

					This->optlog_itr = -1;
				}
				{
					std::ofstream ofs(FILENAME_PARAM, std::ios_base::app);
					ofs << i << "," << j << "," << k << ",,";
					for (int i = 0; i < c->Pcnt; i++)
					{
						ofs << c->Pbl[i] << ",";
					}
					ofs << ",";
					for (int i = 0; i < c->Pcnt; i++)
					{
						ofs << c->Pbr[i] << ",";
					}
					ofs << ",";
					for (int i = 0; i < c->Pcnt; i++)
					{
						ofs << c->Pa[i] << ",";
					}
					ofs << ",";
					for (int i = 0; i < c->Pcnt; i++)
					{
						ofs << c->Px[i] << ",";
					}
					ofs << ",";
					for (int i = 0; i < c->Pcnt; i++)
					{
						ofs << c->Py[i] << ",";
					}
					ofs << std::endl;
					ofs.close();
				}
				This->optlog_itr = -1;
			} // k
		} // j
	} // i
#endif
}

void ControlPanel::cb_btn_batchstop(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	Fl::remove_idle(ControlPanel::idle_batchproc, This);
}
