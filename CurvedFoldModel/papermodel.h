#ifndef CURVEFOLD
#define CURVEFOLD

#include "crease.h"
#include "curvedraw.h"
#include "curveintersect.h"

#define MAXCN 10	// max # of curves
#define MAX_CRS_CNT 10
#define MAX_CRS_OTYPE 5
#define MAX_STP_CNT 300
#define MAX_TGT_CNT 100

#define MAX_PL_FCNT 1000	// MAX # of faces on polygon, > MAX_SPCNT*(MAXCN+1)
#define MAX_PL_ECNT 250		// MAX # of edges on polygon, > MAX_SPCNT*2+MAXCN*2

#define PPWIDTH  200	// 紙の幅（初期値）
#define PPHEIGHT 200	// 紙の高さ（初期値）

typedef enum _postproc_type {
	PPTYPE_X,
	PPTYPE_UNDEF,
	PPTYPE_OPEN,		//① crease::calc****()					立ち上げ時
	PPTYPE_PRICURVE,	//② crease::calc****()					曲線パラメータ修正
	PPTYPE_RTCURVE,		//③ crease:: calcXA_CP()				2Dで曲線位置角度を変更
	PPTYPE_ADDCREASE,	//④ ControlPanel:: cb_btn_proccurve(..)	"ADD CREASE"
	PPTYPE_FMOT,		//⑤ crease::calc****()					folding motion
	PPTYPE_MOVECP0,		//⑥ GraphWindowCP::handle(int event)		control point 位置修正（マウス）
	PPTYPE_MOVECP1,		//⑦ GraphWindowCP::handle(int event)		control point 位置修正（マウス＋最適化）
	PPTYPE_OPTIMIZE,	//⑧ ControlPanel::cb_btn_rectifyC(..)	"OPTIMIZE"
	PPTYPE_PPEDGE		//⑨ GraphWindowCP::handle(int event)		紙端修正
} postproc_type;

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

	// -----------------------------

	// Polygon
	int plcnt, plvcnt[MAX_PL_FCNT];
	double plx[MAX_PL_FCNT*4], ply[MAX_PL_FCNT*4], plz[MAX_PL_FCNT*4];
	double plx_cp[MAX_PL_FCNT*4], ply_cp[MAX_PL_FCNT*4];
	double plmat[MAX_PL_FCNT*16];
	int pl_cridx[MAX_PL_FCNT];

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

	int dumpsvg00( char *fname, int divnum ); // 折り線×divnum
	int dumpsvg0( char *fname ); // 折り線
	int dumpsvg1( char *fname ); // 折り線＋rulings
	int dumpObj(  char *fname, double scl ); // 折り線部分の3D形状

	int check180( char *fname );	// 頂点周りの角度
	int checkquadplane( char *fname );	// quad の平面度
	int checkGap( char *fname, int divnum );	// ピース間の隙間

	// ------------------------- evaluate -------------------------------------------

	int checkGap( double *errdata, int row, int col, int divnum );
	int checkCollision( int *errdata, int row, int col, int divnum );

	// ------------------------- procseq -------------------------------------------

	int postproc();
	int resetpostproc( int _cidx );

	// ------------------------- spherepart -------------------------------------------

	int setSpherePieceBoundary( int divnum, double crv, curve_type cype );

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

	// ------------------------- poly -------------------------------------------

	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d, int cridx );
	int addPly( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double x3_cp, double y3_cp,
		double x0_3d, double y0_3d, double z0_3d, double x1_3d, double y1_3d, double z1_3d,
		double x2_3d, double y2_3d, double z2_3d, double x3_3d, double y3_3d, double z3_3d, int cridx );
	int addPlymat( double x0_cp, double y0_cp, double x1_cp, double y1_cp, double x2_cp, double y2_cp, double *mat, int cridx );
	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1 );
	int addEdge( double x0, double y0, double z0, double x1, double y1, double z1,
		double x2, double y2, double z2 );
	int addEdge1( double x0, double y0, double z0, double x1, double y1, double z1 );
	int calcRulingPly();	// 2D,3Dポリゴン、変換行列作成

	int makeObj3( double thickness );	// polyベース、厚み一定（全頂点同じ方向）、後で重複頂点削除

	// ------------------------- stitch -------------------------------------------

	// stitch
	int stpcnt;
	double stx[2][MAX_STP_CNT], sty[2][MAX_STP_CNT], stz[2][MAX_STP_CNT];
	double stx_cp[2][MAX_STP_CNT], sty_cp[2][MAX_STP_CNT];
	double stx0[2][MAX_STP_CNT], sty0[2][MAX_STP_CNT], stz0[2][MAX_STP_CNT];
	double stx1[2][MAX_STP_CNT], sty1[2][MAX_STP_CNT], stz1[2][MAX_STP_CNT];
	double stxm[2][MAX_STP_CNT], stym[2][MAX_STP_CNT], stzm[2][MAX_STP_CNT];

	// length at ruling end
	int re_sidx, re_eidx;
	double Pre[MAX_CPCNT];
	double re_clen[MAX_SPCNT];
	int setP_re( int Pidx ); // crlen_rulend[] -> Pre[]
	//int interpolate_spline_RulEnd( double *_P, int _Pcnt, double *_X, int _Xcnt );

	int getMat_stitch( int divnum, int maxpcnt, double *mObj );
	int getMat_stitch1( int divnum, int maxpcnt, double *mObj );
	int getParam_stitch( int divnum, int rulingOnly );
	int getParam_calcRE();
	int getParam_optimize( int divnum );

	// ------------------------- optFold -------------------------------------------

	int tgcnt;
	double tgx[MAX_TGT_CNT], tgy[MAX_TGT_CNT], tgz[MAX_TGT_CNT];
	double ogx[MAX_TGT_CNT], ogy[MAX_TGT_CNT], ogz[MAX_TGT_CNT];
	double ogx_cp[MAX_TGT_CNT], ogy_cp[MAX_TGT_CNT];
	double tgap[MAX_TGT_CNT], avetgap;

	int saveTgt( char *fname );
	int loadTgt( char *fname );
	int getTgt2D3D();
	double calcRulCross( double *xx, double *xy, double *rx, double *ry, double *rlen, int cvcnt );
	int checkCollision();
	int checkRulAngle();
	int checkRulCross();
	int checkRulCross( double &_area );
	int checkRulCreaseCross();
	int optMat( int mode );
	int optFold();
	int optTorsion();
};

#endif
