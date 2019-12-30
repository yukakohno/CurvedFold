#define _USE_MATH_DEFINES
#include <math.h>
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

void ControlPanel::cb_btn_eval_gap( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	int divnum = This->vs_divnum->value();

	double errdata[ MAX_STP_CNT*7 ];
	memset( errdata, 0, sizeof(double)*MAX_STP_CNT*7 );
	int spcnt = ppm->checkGap( errdata, MAX_STP_CNT, 7, divnum );
	int cnt=0;
	double ave=0.0;
	printf("error:");
	for( int i=0; i<spcnt; i++ ){
		if( errdata[i*7] ){
			ave += errdata[i*7];
			cnt++;
		}
		printf(" %.2f", errdata[i*7]);
	}
	ave/=(double)cnt;
	printf("\naverage error = %f\n", ave);
}

void ControlPanel::cb_btn_eval_collision( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	int divnum = This->vs_divnum->value();

	int errdata[ MAX_PL_FCNT ];
	memset( errdata, 0, sizeof(int)*MAX_PL_FCNT );
	int spcnt = ppm->checkCollision( errdata, MAX_PL_FCNT, 1, divnum );

	printf("collision:");
	int total=0;
	for( int i=0; i<ppm->plcnt; i++ ){
		total += errdata[i];
		printf(" %d", errdata[i] );
	}
	printf("\n# collision face = %d\n", total );
}

void ControlPanel::cb_btn_eval_rulingcross( Fl_Widget *wgt, void *idx )
{
	ControlPanel *This = (ControlPanel *)idx;
	papermodel *ppm = &(This->ppm);
	crease *c0 = &(ppm->crs[0]);

	double errdata[ MAX_SPCNT*2 ];
	memset( errdata, 0, sizeof(double)*MAX_SPCNT*2 );
	c0->checkRulingCross( errdata, MAX_SPCNT, 2 );

	double area=0.0;
	printf("crossing area (left):");
	for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
		area += errdata[i*2];
		printf(" %.2f", errdata[i*2] );
	}
	printf("\ncrossing area (right):");
	for( int i=c0->Xsidx; i<=c0->Xeidx; i++ ){
		area += errdata[i*2+1];
		printf(" %.2f", errdata[i*2+1] );
	}
	printf("\ncrossing area = %f\n", area);
}
