#ifndef CURVEFOLD
#define CURVEFOLD

#include "crease.h"
#include "curvedraw.h"
#include "curveintersect.h"

#define MAXCN 10	// max # of curves
#define MAX_CRS_CNT 10
#define MAX_CRS_OTYPE 5

#define MAX_PL_FCNT 1000	// MAX # of faces on polygon, > MAX_SPCNT*(MAXCN+1)
#define MAX_PL_ECNT 250		// MAX # of edges on polygon, > MAX_SPCNT*2+MAXCN*2

#define PPWIDTH  200	// ���̕��i�����l�j
#define PPHEIGHT 200	// ���̍����i�����l�j

typedef enum _postproc_type {
	PPTYPE_X,
	PPTYPE_UNDEF,
	PPTYPE_OPEN,		//�@ crease::calc****()					�����グ��
	PPTYPE_PRICURVE,	//�A crease::calc****()					�Ȑ��p�����[�^�C��
	PPTYPE_RTCURVE,		//�B crease:: calcXA_CP()				2D�ŋȐ��ʒu�p�x��ύX
	PPTYPE_ADDCREASE,	//�C ControlPanel:: cb_btn_proccurve(..)	"ADD CREASE"
	PPTYPE_FMOT,		//�D crease::calc****()					folding motion
	PPTYPE_PPEDGE		//�H GraphWindowCP::handle(int event)		���[�C��
} postproc_type;

//typedef enum _evaltype { EVTYPE_TORSION, EVTYPE_BETALR, EVTYPE_RULCROSS } evaltype;

class papermodel
{
public:
	postproc_type flg_postproc_type;
	int flg_interpolate;	// 0: use X*[] as is, 1: re-generate X*[] from P*[]
	rectify_params rp;

	// for B-spline, array size=1000
	double *hSx, *hSy, *hL, *ht1;

	// papar size & edge
	int pw,ph, psx,psy,pex,pey;

	// crease
	int crcnt;
	crease crs[MAX_CRS_CNT];

	// -----------------------------

	// curve
	int dccnt;
	curvedraw dcurve[MAXCN];
	int tccnt;
	curveintersect tcurve[MAXCN];

	// -----------------------------

	// Polygon
	int plcnt, plvcnt[MAX_PL_FCNT];
	double plx[MAX_PL_FCNT*4], ply[MAX_PL_FCNT*4], plz[MAX_PL_FCNT*4];
	double plx_cp[MAX_PL_FCNT*4], ply_cp[MAX_PL_FCNT*4];
	double plmat[MAX_PL_FCNT*16];

	// PolyEdge(side)
	int plecnt;
	double plex[MAX_PL_ECNT*2], pley[MAX_PL_ECNT*2], plez[MAX_PL_ECNT*2];

	// Obj
	int o_vcnt, o_fcnt, o_fvcnt[MAX_PL_FCNT];
	double o_vx[MAX_PL_FCNT*4], o_vy[MAX_PL_FCNT*4], o_vz[MAX_PL_FCNT*4];
	int o_fi[MAX_PL_FCNT][4];

	// -------------------------------------------------------------------

	papermodel();
	~papermodel();
	void init();
	void initX();
	void set_postproc_type( postproc_type ppt );

	// ------------------------- FILEIO -------------------------------------------

	int dumpsvg0( char *fname ); // �܂��
	int dumpsvg1( char *fname ); // �܂���{rulings
	int dumpObj(  char *fname, double scl ); // �܂��������3D�`��

	int dumpcv0( char *fname );
	int loadcv0( char *fname );

	int check180( char *fname );	// ���_����̊p�x
	int checkquadplane( char *fname );	// quad �̕��ʓx

	// ------------------------- procseq -------------------------------------------

	int postproc();
	int resetpostproc( int _cidx );

	// ------------------------- poly -------------------------------------------

	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d );
	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double x3_cp, double y3_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d, double x3_3d, double y3_3d, double z3_3d );
	int addPlymat( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double *mat );
	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1 );
	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1,
		double x2, double y2, double z2 );
	int addEdge1( double x0, double y0, double z0, double x1, double y1, double z1 );
	int calcRulingPly();	// 2D,3D�|���S���A�ϊ��s��쐬

	int makeObj3( double thickness );	// poly�x�[�X�A���݈��i�S���_���������j�A��ŏd�����_�폜

};

#endif
