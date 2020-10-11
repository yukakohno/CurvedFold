#include "ControlPanel.h"

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
