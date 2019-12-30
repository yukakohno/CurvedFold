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

#define PPWIDTH  200	// 紙の幅（初期値）
#define PPHEIGHT 200	// 紙の高さ（初期値）

typedef enum _postproc_type {
	PPTYPE_X,
	PPTYPE_UNDEF,
	PPTYPE_OPEN,		//�@ crease::calc****()					立ち上げ時
	PPTYPE_PRICURVE,	//�A crease::calc****()					曲線パラメータ修正
	PPTYPE_RTCURVE,		//�B crease:: calcXA_CP()				2Dで曲線位置角度を変更
	PPTYPE_ADDCREASE,	//�C ControlPanel:: cb_btn_proccurve(..)	"ADD CREASE"
	PPTYPE_FMOT,		//�D crease::calc****()					folding motion
	PPTYPE_MOVECP0,		//�E GraphWindowCP::handle(int event)		control point 位置修正（マウス）
	PPTYPE_MOVECP1,		//�F GraphWindowCP::handle(int event)		control point 位置修正（マウス＋最適化）
	PPTYPE_OPTIMIZE,	//�G ControlPanel::cb_btn_rectifyC(..)	"OPTIMIZE"
	PPTYPE_PPEDGE		//�H GraphWindowCP::handle(int event)		紙端修正
} postproc_type;

typedef enum _evaltype { EVTYPE_TORSION, EVTYPE_BETALR, EVTYPE_RULCROSS } evaltype;

class papermodel
{
public:
	postproc_type flg_postproc_type;
	int flg_interpolate;	// 0: use X*[] as is, 1: re-generate X*[] from P*[]
	int flg_rectifyA, flg_rectifyT, flg_rectifyR;	// 1: rectify alpha, torsion, ruling length

	// for B-spline, array size=1000
	double *hSx, *hSy, *hL, *ht1;

	// papar size & edge
	int pw,ph, psx,psy,pex,pey;

	// crease
	int crcnt, lcrcnt, rcrcnt;
	crease crs[MAX_CRS_CNT], *lcrs[MAX_CRS_CNT], *rcrs[MAX_CRS_CNT];

	// -----------------------------

	// curve
	int dccnt;
	curvedraw dcurve[MAXCN];
	int fccnt;
	curveintersect fcurve[MAXCN];
	int tccnt;
	curveintersect tcurve[MAXCN];

	// rectify crease
	evaltype hCrs_evaltype;
	int flg_rectifyCrs;				// the crease to be optimized in optRulLen_SplineCP(), set by scroll bar
	int hCrs_cridx, hCrs_cpidx;		// the crease and CP to be optimized in optRulLen_SplineCP1(), selected by mouse
	double hCrs_cpdx, hCrs_cpdy;			// optRulLen_SplineCP1()
	double hCrs_mindiff[MAX_CRS_OTYPE];		// （暫定）最小値
	double hCrs_mgndiff[MAX_CRS_OTYPE];		// コスト増加の許容値
	double hCrs_thresdiff[MAX_CRS_OTYPE];	// これ以下になったら最適化終了
	double cverr[MAX_SPCNT], cvminerr[MAX_SPCNT];	// 表示用
	double hCrs_cvmaxerr[MAX_CRS_OTYPE];	// 表示用、↑にかける

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
	void initHCrs();
	void set_postproc_type( postproc_type ppt );

	// ------------------------- FILEIO -------------------------------------------

	int dumpsvg0( char *fname ); // 折り線
	int dumpsvg1( char *fname ); // 折り線＋rulings
	int dumpObj(  char *fname, double scl ); // 折り線部分の3D形状

	int dumpcv0( char *fname );
	int loadcv0( char *fname );
	//int dumpcv1( char *fname );
	//int loadcv1( char *fname );
	int dumpBsplineCP( char *fname );
	int loadBsplineCP( char *fname );

	int check180( char *fname );	// 頂点周りの角度
	int checkquatplane( char *fname );	// quad の平面度

	// ------------------------- procseq -------------------------------------------

	int postproc();
	int resetpostproc( int _cidx );

	// ------------------------- curve -------------------------------------------

	int CP2DC();
	int FC2Crs();
	int addCrease( int lcrs_idx, int rcrs_idx, double slen, double dlen );

	// ------------------------- rulLength -------------------------------------------

	// crease0 のruling長さを crease1 により切り取る
	int setRulLen_curve( double *cvx, double *cvy, double *rx, double *ry, int cvcnt, // crease0 &ruling
		double *cvx1, double *cvy1, int cvcnt1, // crease1 (fold curve)
		double *_Cx, double *_Cy, // control point of crease1 (output)
		double *rlen, ruling_end_type *rflg ); // length of crease 0 ruling (output)

	int setRulLen_SplineCP( double *cvx, double *cvy, double *rx, double *ry, int cvcnt, // crease0 &ruling
		double *_Cx, double *_Cy, // control point of crease1
		double *rlen, ruling_end_type *rflg ); // length of crease 0 ruling (output)
	int setRulLen_SplineCP( double *cvx, double *cvy, double *rx, double *ry, int si, int ei, // crease0 &ruling
		double *_Cx, double *_Cy, // control point of crease1
		double *rlen, ruling_end_type *rflg ); // length of crease 0 ruling (output)

	// ------------------------- optRulLength -------------------------------------------

	double calcAveDiff( double *tr0, double *tr1, double max, int cvcnt );
	double calcRulCross( double *xx, double *xy, double *rx, double *ry, double *rlen, int cvcnt );

	int set_rectifyCrs12( int cridx, int cpidx, double dx, double dy );
	int set_resetCrs12();

	int optRulLen_SplineCP( crease *c0, crease *c1 );
	int optRulLen_SplineCP1( crease *c0, crease *c1, int cpidx, double cpdx, double cpdy );

	// ------------------------- poly -------------------------------------------

	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d );
	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double x3_cp, double y3_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d, double x3_3d, double y3_3d, double z3_3d );

	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double x3_cp, double y3_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d, double *x3_3d, double *y3_3d, double *z3_3d );

	int addPlymat( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double *mat );
	int addPly( double *x_cp, double *y_cp, double *x_3d, double *y_3d, double *z_3d, int cnt );

	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1 );
	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1,
		double x2, double y2, double z2 );
	int addEdge1( double x0, double y0, double z0, double x1, double y1, double z1 );
	int calcRulingPly();	// 2D,3Dポリゴン、変換行列作成

	int makeObj3( double thickness );	// polyベース、厚み一定（全頂点同じ方向）、後で重複頂点削除

};

#endif
