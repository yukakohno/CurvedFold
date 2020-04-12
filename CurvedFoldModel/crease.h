#ifndef CREASE
#define CREASE

#include <stdio.h>

#define MAX_FRAME 41
#define MAX_CPCNT 10	// MAX # of control points on the crese
#define MAX_SPCNT 100	// MAX # of sample points on the crese

#define XCNT_DEF 40		// default # of sample points on the crese
#define XSPC_DEF 5.0	// default space between the sample points
#define CEMGN 3			// default end margin of the crease curve

#define CCNT 6 // # of control points of B spline

#define FM_MODE 2	// 1: Control Point, 2: Sample Point

typedef enum _crease_end_type {
		CETYPE_UNDEF, CETYPE_TRIM,
		CETYPE_EGTOP0, CETYPE_EGRIGHT0, CETYPE_EGBOTTOM0, CETYPE_EGLEFT0,	// paper edge
		CETYPE_EGTOP1, CETYPE_EGRIGHT1, CETYPE_EGBOTTOM1, CETYPE_EGLEFT1	// paper edge
} crease_end_type;

typedef enum _ruling_end_type {
	RETYPE_UNDEF, RETYPE_CURVE, RETYPE_TRIM, RETYPE_EGTOP, RETYPE_EGRIGHT, RETYPE_EGBOTTOM, RETYPE_EGLEFT
} ruling_end_type;

#define MIN_CURV 0.005
#define RECT_ASIZE 10

typedef struct{
	int flg_src_s1e1;	// 0:k2d, 1:kv*cos(alpha)
	double kvthres;		// 曲率0とみなす閾値
	int s1mgn;			// 曲率0区間を広げる幅 >=0
	int s2mgn;			// 補間区間  -1: 極値まで, 0: 補間なし, >0: 補間幅
	int method;			// 補間方法  1: linear, 2:Bezier
	int src;			// 0: すべて算出  1: s1e1は既存, 2: s2e2は既存
	int s1[RECT_ASIZE], e1[RECT_ASIZE], s2[RECT_ASIZE], e2[RECT_ASIZE], scnt;
	double svlen0[RECT_ASIZE], svlen1[RECT_ASIZE];	// (s1-s2)*0.5
	double evlen0[RECT_ASIZE], evlen1[RECT_ASIZE];	// (e2-e1)*0.5
} rectify_param;

typedef struct{
	bool flg_rectifyT;
	rectify_param rectT;
	bool flg_rectifyA;
	rectify_param rectA;
	bool flg_rectifyR;
	double rectifyR_kvthres;
} rectify_params;

class crease
{
public:
	int rl;				// -1:left, 1:right
	int flg_org;		// generated from 0:org* 1:CP*
	//int flg_src_s1e1;	// 0: k2d, 1: kv*cos(alpha)
	int flg_xspc_def;	// 1: XSPC_DEF, 0: else

	// original curve shape (reference to curvedraw or curveintersect)
	int org_idx;
	int org_cnt;
	double *org_x, *org_y;

	// control points of B-spline curve
	double CPx[CCNT], CPy[CCNT];

	// control points on the crease
	int Pcnt;
	double Px2d[MAX_CPCNT], Py2d[MAX_CPCNT]; // Px2d <- k2d
	double Px[MAX_CPCNT], Py[MAX_CPCNT], Pz[MAX_CPCNT]; // Px <- kv, Py <- tr
	double Pa[MAX_CPCNT];	// alpha
	double Pbl[MAX_CPCNT], Pbr[MAX_CPCNT];

	// foldmotion
	double Px2d_org[MAX_CPCNT], Px_org[MAX_CPCNT], Py_org[MAX_CPCNT], Pa_org[MAX_CPCNT], m3_org[16];
	int FM_fidx_org, FM_fidx_min, FM_fidx_max;
	int flg_FM_Pt[MAX_FRAME], flg_FM_Pa[MAX_FRAME], flg_FM_m3[MAX_FRAME];
	double FM_Pt[MAX_FRAME][MAX_CPCNT];	// tr
	double FM_Pa[MAX_FRAME][MAX_CPCNT];	// alpha
	double FM_tr[MAX_FRAME][MAX_SPCNT];	// tr
	double FM_alpha[MAX_FRAME][MAX_SPCNT];	// alpha
	double FM_m3[MAX_FRAME][16];

	// position & rotation (3D, CP)
	double m3[16], m2[9];

	// crease curve
	int Xcnt, XcntOrg;
	double XspcOrg;

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

	int loadMotion( char *fname );	// return -1: error, 0: flg_usecp=0, 1: flg_usecp=1
	int dumpMotion( char *fname, int flg_usecp );
	int loadMotionFrame( int frm, char *fname );	// return -1: error, 0: flg_usecp=0, 1: flg_usecp=1
	int dumpMotionFrame( int frm, char *fname, int flg_usecp );

	int loadPb( char *fname );
	int dumpPb( char *fname );

	int check180( char *fname );	// 頂点周りの角度
	int check180( FILE *fp );
	int checkquadplane( char *fname );	// quad の平面度
	int checkquadplane( FILE *fp );
	int checkCP3D( char *fname );

	// ------------------------- evaluate -------------------------------------------

	int check180( double *errdata, int row, int col );
	int checkquadplane( double *errdata, int row, int col );
	int checkRulingCross( double *errdata, int row, int col );
	int checkRulingCross( int rl ); // -1:left, 1:right

	// ------------------------- proseq -------------------------------------------

	int calcXA_CP( int flg_interpolate, rectify_params *rp );
	int calcCPA_X( int flg_interpolate, rectify_params *rp );
	int calcCPX_A( int flg_interpolate, rectify_params *rp );
	int calcR_TA( int flg_interpolate, rectify_params *rp, int mini, int maxi, double ma );
	int calcR_TA2( rectify_params *rp, int mini, int maxi ); // calcR_TA( flg_interpolate=0 ) と同じなので不使用

	// ------------------------- interpolate -------------------------------------------

	static int interpolate_spline( double *P, int Pcnt, double *X, int Xcnt ); // スプライン補間
	static int interpolate_spline_RulAngle( double *_P, int _Pcnt,
		double *_B, int _Xcnt, double *cotB, double *cosB, double *sinB );

	int setP_k( int Pidx ); // kv -> Px[_Pcnt]
	int setP_t( int Pidx ); // tr -> Py[_Pcnt]
	int setP_a( int Pidx );  // alpha -> Pa[_Pcnt]
	int setP_k2( int Pidx ); // kt -> Px2d[_Pcnt]
	int setP_Bl( int Pidx ); // betal -> Pbl[_Pcnt]
	int setP_Br( int Pidx ); // betar -> Pbr[_Pcnt]

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

	int calcTNB(); // X -> T,N,B,kv,tr
	int calcTNB( int Xsidx0, int Xeidx0 );
	static int calcTNB( int Xcnt, double *Xx, double *Xy, double *Xz,
		double *Tx, double *Ty, double *Tz,	double *dx,
		double *Nx, double *Ny, double *Nz, double *kv,
		double *Bx, double *By, double *Bz, double *tr );
	int calcXTNB( double m3[16] ); // kv,tr -> (rectifyTau) -> X -> calcTNB();

	int calcTau(); // tr算出
	static int calcTau( int Xcnt, double *Xx, double *Xy, double *Xz,
		double *Tx, double *Ty, double *Tz, double *dx,
		double *Nx, double *Ny, double *Nz, double *kv,
		double *Bx, double *By, double *Bz, double *tr );

	// ------------------------- alpha -------------------------------------------

	int setA(double a); // alpha = a を設定
	int calcAlpha( int flg_rectifyA, rectify_param *rp ); // kv, k2d -> alpha, sina, cosa, tana, da
	static int calcAlpha( int Xcnt, double *Xx, double *Xy, double *Xz, double *kv, double *k2d,
		double *cosa, double *sina, double *tana, double *alpha, double *da,
		int flg_rectifyA, rectify_param *rp );
	int calcAK_K2D(); // alpha, kv -> k2d, sina, cosa, tana, da
	int calcAK2D_K(); // alpha, k2d -> kv
	int calcDA(); // alpha, X -> da

	int calcRul2TA( int besteffort, int mini, int maxi ); // papermodel->tr,alpha

	// ------------------------- beta -------------------------------------------

	int calcRuling( int flg_rectifyR, double rectifyR_kvthres ); // ruling算出
	static int calcRuling( int cxcnt,
		double *Tx, double *Ty, double *Tz, double *Tx2d, double *Ty2d, 
		double *Nx, double *Ny, double *Nz,  double *Bx, double *By, double *Bz,					   
		double *da, double *kv, double *tr,
		double *betal, double *betar, double *sinbl, double *sinbr,
		double *cosbl, double *cosbr, double *cotbl, double *cotbr, double *sina, double *cosa, 
		double *rlx, double *rly, double *rlz, double *rlx_cp, double *rly_cp, double *rllen,
		double *rrx, double *rry, double *rrz, double *rrx_cp, double *rry_cp, double *rrlen, int rl );
	int calcRuling2D( int flg_rectifyR, double rectifyR_kvthres );
	static int calcRuling2D( int cxcnt, double *da, double *sina, double *kv, double *tr,
		double *betal, double *betar, double *sinbl, double *sinbr,
		double *cosbl, double *cosbr, double *cotbl, double *cotbr, int rl );
	int calcRuling3D();
	static int calcRuling3D( int cvcnt, double *Tx, double *Ty, double *Tz, double *Tx2d, double *Ty2d, 
		double *Nx, double *Ny, double *Nz,  double *Bx, double *By, double *Bz,
		double *sina, double *cosa, double *sinbl, double *sinbr, double *cosbl, double *cosbr,
		double *rlx, double *rly, double *rlz, double *rrx, double *rry, double *rrz,
		double *rlx_cp, double *rly_cp, double *rrx_cp, double *rry_cp, int rl );

	int calcRuling2DH( int flg_rectifyR, double rectifyR_kvthres ); // ruling horizontal

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

	int calcRLenHoseiR( double rectifyR_kvthres ); // 曲率小さい部分は長さ0に
	static int calcRLenHoseiR( int Xcnt, double *k2d, double rectifyR_kvthres, double *rlen, ruling_end_type *rflg ); // 曲率小さい部分は長さ0に

	int setEnds( int psx, int psy, int pex, int pey, int dccnt, void *dcurve );
	void initEnds();

	// ------------------------- rectify -------------------------------------------

	static void setDefault( rectify_param *rp );

	int gets1e1( int flg_src_s1e1, double thres, int mgn, int *_s1, int *_e1, int semax );
	static int gets1e1( int Xcnt, int Xs, int Xe, double *val, double thres, int mgn,
		int *_s1, int *_e1, int semax );
	static int gets2e2( int Xcnt, int *s1, int *e1, int scnt, double *alpha, int *s2, int *e2 );

	int rectifyTau( double rectifyT_kvthres, int rectifyT_mgn );
	int rectifyTau2( int src_s1e1, double rectifyT_kvthres, int rectifyT_mgn ); // 積分値が同じになるように
	static int rectifyTau( int Xcnt, int *s1, int *e1, int scnt, double *tr, double *kv );
	int rectifyTauBezier( rectify_param *rp );
	int rectifyTauBezier2( rectify_param *rp );

	int rectifyAlpha( rectify_param *rp );
	int rectifyAlphaBezier( rectify_param *rp );
	static int rectifyAlphaBezier( int Xcnt, double *alpha0, rectify_param *rp );

	// ------------------------- RL -------------------------------------------
	int setL2R( crease *c ); // copy c->rl* to this->rr*, set Xsidx Xeidx
	int setR2L( crease *c ); // copy c->rr* to this->rl*, set Xsidx Xeidx
	int calcLeft();
	int calcRight();
	int calcLeftK();	// 
	int calcRightK();	// 

	// ------------------------- foldmotion -------------------------------------------

	void initFM();
	int updateFM( int flg_interpolateCP );

	// ------------------------- stitch -------------------------------------------

	int getSamplePoints2D( double space, double *x, double *y ); // 不使用、使う場合は要更新
	int getSamplePoints3D( double space, double *x, double *y, double *z, double *len_sp );
	int getVertexPoints3D( double space, double *x, double *y, double *z, double *len_sp );

};

#endif