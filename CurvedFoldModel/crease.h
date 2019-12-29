#ifndef CREASE
#define CREASE

#include <stdio.h>

#define MAX_CPCNT 10	// MAX # of control points on the crese
#define MAX_SPCNT 100	// MAX # of sample points on the crese

#define XCNT 40		// default # of sample points on the crese
#define XSPC 5.0	// default space between the sample points
#define CEMGN 3		// default end margin of the crease curve

#define CCNT 6 // # of control points of B spline

#define MIN_CURV 0.005

typedef enum _crease_end_type {
		CETYPE_UNDEF, CETYPE_TRIM,
		CETYPE_EGTOP0, CETYPE_EGRIGHT0, CETYPE_EGBOTTOM0, CETYPE_EGLEFT0,	// paper edge
		CETYPE_EGTOP1, CETYPE_EGRIGHT1, CETYPE_EGBOTTOM1, CETYPE_EGLEFT1	// paper edge
} crease_end_type;

typedef enum _ruling_end_type {
	RETYPE_UNDEF, RETYPE_CURVE, RETYPE_TRIM, RETYPE_EGTOP, RETYPE_EGRIGHT, RETYPE_EGBOTTOM, RETYPE_EGLEFT
} ruling_end_type;

// defined in curvedraw.h
//typedef enum _curve_type {
//	CTYPE_UNDEF, CTYPE_TRIM, CTYPE_FOLD, CTYPE_FOLD_LEFT, CTYPE_FOLD_RIGHT
//} curve_type;
//typedef enum _cvertex_type {
//	CVTYPE_UNDEF, CVTYPE_RUL_LEFT, CVTYPE_RUL_RIGHT, CVTYPE_FCURVE_START, CVTYPE_FCURVE_END, CVTYPE_PAPER_EDGE
//} cvertex_type;

class crease
{
public:
	int rl;				// -1:left, 1:right
	int flg_org;		// generated from 0:org* 1:CP*
	int flg_src_s1e1;	// 0: k2d, 1: kv*cos(alpha)

	// original curve shape (reference to curvedraw or curveintersect)
	int org_idx;
	int org_cnt;
	double *org_x, *org_y;

	// control points on the crease
	int Pcnt;
	double Px2d[MAX_CPCNT], Py2d[MAX_CPCNT]; // Px2d <- k2d
	double Px[MAX_CPCNT], Py[MAX_CPCNT], Pz[MAX_CPCNT]; // Px <- kv, Py <- tr
	double Pa[MAX_CPCNT];	// alpha

	// position & rotation (3D, CP)
	double m3[16], m2[9];

	// crease curve
	int Xcnt;

	// crease pattern
	double Xx2d[MAX_SPCNT], Xy2d[MAX_SPCNT];
	double Tx2d[MAX_SPCNT], Ty2d[MAX_SPCNT], d2d[MAX_SPCNT];
	double Nx2d[MAX_SPCNT], Ny2d[MAX_SPCNT], k2d[MAX_SPCNT];

	// 3d curve
	double Xx[MAX_SPCNT], Xy[MAX_SPCNT], Xz[MAX_SPCNT];
	double Tx[MAX_SPCNT], Ty[MAX_SPCNT], Tz[MAX_SPCNT], dx[MAX_SPCNT];
	double Nx[MAX_SPCNT], Ny[MAX_SPCNT], Nz[MAX_SPCNT], kv[MAX_SPCNT];
	double Bx[MAX_SPCNT], By[MAX_SPCNT], Bz[MAX_SPCNT], tr[MAX_SPCNT];

	// folding angles
	double alpha[MAX_SPCNT], cosa[MAX_SPCNT], sina[MAX_SPCNT], tana[MAX_SPCNT], da[MAX_SPCNT];

	// rulings
	double betal[MAX_SPCNT], betar[MAX_SPCNT], cotbl[MAX_SPCNT], cotbr[MAX_SPCNT];
	double cosbl[MAX_SPCNT], sinbl[MAX_SPCNT], cosbr[MAX_SPCNT], sinbr[MAX_SPCNT];
	double rlx[MAX_SPCNT], rly[MAX_SPCNT], rlz[MAX_SPCNT], rlx_cp[MAX_SPCNT], rly_cp[MAX_SPCNT];
	double rrx[MAX_SPCNT], rry[MAX_SPCNT], rrz[MAX_SPCNT], rrx_cp[MAX_SPCNT], rry_cp[MAX_SPCNT];
	ruling_end_type rlflg[MAX_SPCNT], rrflg[MAX_SPCNT];
	double rllen[MAX_SPCNT], rrlen[MAX_SPCNT], maxrlen;

	// end points of crease curves
	int Xsidx, Xeidx;
	double Xxs2d, Xys2d, Xxs, Xys, Xzs, Xxs2d0, Xys2d0;
	double Xxe2d, Xye2d, Xxe, Xye, Xze, Xxe2d0, Xye2d0;
	crease_end_type Xstype, Xetype;

	// used in dumpAll()
	int fileioflg[30];
	enum{ PX2, XX2, TX2, NX2, D2, K2, PX3, XX3, TX3, NX3, BX3, D3, K3, T3, PA, ALPHA, ALPHA0, DA, PB, BETA, BETA0, R3, R2 };

	crease();
	~crease();
	void init();
	void initX();
	int copy( crease *c );
	void init_rlflg();
	void init_rrflg();

	int org2CP();

	// ------------------------- FILEIO -------------------------------------------

	int load( char *fname );  // P*をファイル読み込み
	int dumpP( char *fname, int mode );
	int dumpP_( char *fname, int mode );
	int dumpX( char *fname );
	int loadm2m3( char *fname );
	int dumpm2m3( char *fname );
	int dumpAll( char *fname );

	int check180( char *fname );	// 頂点周りの角度
	int check180( FILE *fp );
	int checkquadplane( char *fname );	// quad の平面度
	int checkquadplane( FILE *fp );
	int checkCP3D( char *fname );

	// ------------------------- proseq -------------------------------------------

	int calcXA_CP( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR );	// a) 3D曲線と折り角度から、折り線を求める
	int calcCPA_X( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR );	// b) 折り線と折り角度から、3D曲線を求める
	int calcCPX_A( int flg_interpolate, int flg_rectifyT, int flg_rectifyA, int flg_rectifyR );	// c) 折り線と3D曲線から、折り角度を求める

	// ------------------------- interpolate -------------------------------------------

	static int interpolate_spline( double *P, int Pcnt, double *X, int Xcnt ); // スプライン補間

	int setP_k( int Pidx ); // kv -> Px[_Pcnt]
	int setP_t( int Pidx ); // tr -> Py[_Pcnt]
	int setP_a( int Pidx );  // alpha -> Pa[_Pcnt]
	int setP_k2( int Pidx ); // kt -> Px2d[_Pcnt]

	// P*->X* スプライン補間（sに関して等間隔）, mode=0:CP&3D 1:CP 2:3D
	// used only in file open
	int interpolate_spline2( int mode );

	// ------------------------- TNB -------------------------------------------

	int calcTN2d(); // X2d -> T2d,N2d,k2d
	int calcTN2d( int Xsidx0, int Xeidx0 ); // X2d -> T2d,N2d,k2d
	static int calcTN2d(int Xcnt0, double *Xx2d0, double *Xy2d0,
		double *Tx2d0, double *Ty2d0, double *d2d0,
		double *Nx2d0, double *Ny2d0, double *k2d0 );
	int calcXTN2d( double m2[9] ); // k2d -> X2d -> calcTN2d();
	// ★calcXTN2d()を繰り返すと誤差が累積しX2dが変化。

	int calcTNB( int flg_rectifyT ); // X -> T,N,B,kv,tr
	int calcTNB( int flg_rectifyT, int Xsidx0, int Xeidx0 );
	static int calcTNB( int Xcnt, double *Xx, double *Xy, double *Xz,
		double *Tx, double *Ty, double *Tz,	double *dx,
		double *Nx, double *Ny, double *Nz, double *kv,
		double *Bx, double *By, double *Bz, double *tr, int flg_rectifyT );
	int calcXTNB( int flg_rectifyT, double m3[16] ); // kv,tr -> (rectifyTau) -> X -> calcTNB();
	// ★calcXTNB()を繰り返すとtrの値が変化していく -> ruling方向に影響

	int calcTau( int flg_rectifyT ); // tr算出
	static int calcTau( int Xcnt, double *Xx, double *Xy, double *Xz,
		double *Tx, double *Ty, double *Tz, double *dx,
		double *Nx, double *Ny, double *Nz, double *kv,
		double *Bx, double *By, double *Bz, double *tr, int flg_rectifyT );

	// ------------------------- alpha -------------------------------------------

	int setA(double a); // alpha = a を設定
	int calcAlpha( int flg_rectifyA ); // kv, k2d -> alpha, sina, cosa, tana, da
	static int calcAlpha( int Xcnt, double *Xx, double *Xy, double *Xz, double *kv, double *k2d,
		double *cosa, double *sina, double *tana, double *alpha, double *da, int flg_rectifyA );
	int calcAK_K2D( int flg_rectifyA ); // alpha, kv -> (rectifyAlpha) -> k2d, sina, cosa, tana, da
	int calcAK2D_K( int flg_rectifyA ); // alpha, k2d -> (rectifyAlpha) -> kv
	int calcDA(); // alpha, X -> da

	// ------------------------- beta -------------------------------------------

	int calcRuling( int flg_rectifyR ); // ruling算出
	static int calcRuling( int cxcnt,
		double *Tx, double *Ty, double *Tz, double *Tx2d, double *Ty2d, 
		double *Nx, double *Ny, double *Nz,  double *Bx, double *By, double *Bz,					   
		double *da, double *kv, double *tr,
		double *betal, double *betar, double *sinbl, double *sinbr,
		double *cosbl, double *cosbr, double *cotbl, double *cotbr, double *sina, double *cosa, 
		double *rlx, double *rly, double *rlz, double *rlx_cp, double *rly_cp, double *rllen,
		double *rrx, double *rry, double *rrz, double *rrx_cp, double *rry_cp, double *rrlen, int rl );
	int calcRuling2D( int flg_rectifyR );
	static int calcRuling2D( int cxcnt, double *da, double *sina, double *kv, double *tr,
		double *betal, double *betar, double *sinbl, double *sinbr,
		double *cosbl, double *cosbr, double *cotbl, double *cotbr, int rl );
	int calcRuling3D();
	static int calcRuling3D( int cvcnt, double *Tx, double *Ty, double *Tz, double *Tx2d, double *Ty2d, 
		double *Nx, double *Ny, double *Nz,  double *Bx, double *By, double *Bz,
		double *sina, double *cosa, double *sinbl, double *sinbr, double *cosbl, double *cosbr,
		double *rlx, double *rly, double *rlz, double *rrx, double *rry, double *rrz,
		double *rlx_cp, double *rly_cp, double *rrx_cp, double *rry_cp, int rl );

	int calcRuling2DH( int flg_rectifyR ); // ruling horizontal

	// ------------------------- length -------------------------------------------

	void setRLen( double len );
	int calcRLen0(); // ruling長さ
	int calcRLen1_( double *rx_cp, double *ry_cp, double *rlen, char *logfile );
	int calcRLen1(); // ruling長さ

	int calcRLenP( int psx, int psy, int pex, int pey ); // ruling長さを紙の端まで
	static int calcRLenP_( double *cvx, double *cvy, double *rx, double *ry, int cvcnt,
		int psx, int psy, int pex, int pey, double *rlen, ruling_end_type *rflg, char *logfile );

	int calcRLenC( int dccnt, void *dc ); // 曲線トリムがあればruling長さ修正
	static int calcRLenC_( double *cvx, double *cvy, double *rx, double *ry, int cvcnt,
					double *cvx1, double *cvy1, int cvcnt1, double *rlen, ruling_end_type *rflg );

	int calcRLenHoseiR(); // 曲率小さい部分は長さ0に
	static int calcRLenHoseiR( int Xcnt, double *k2d, double *rlen, ruling_end_type *rflg ); // 曲率小さい部分は長さ0に

	int setEnds( int psx, int psy, int pex, int pey, int dccnt, void *dcurve );
	void initEnds();

	// ------------------------- rectify -------------------------------------------

	int gets1e1( int *_s1, int *_e1, int semax, int mgn, int flg_src_s1e1  );
	static int gets1e1( int Xcnt, int Xs, int Xe, double *val, double thres, int mgn,
		int *_s1, int *_e1, int semax );
	int rectifyTau();
	static int rectifyTau( int Xcnt, int *s1, int *e1, int scnt, double *tr, double *kv );
	int rectifyAlpha();
	int rectifyTauBezier();
	int rectifyTauBezier2();
	int rectifyAlphaBezier( int src_s1e1 ); // 0:k2d, 1:kv*cos(alpha)
	static int rectifyAlphaBezier( int Xcnt, int *s1, int *e1, int scnt, double *alpha );
	int rectifyTau2(); // 積分値が同じになるように
};

#endif