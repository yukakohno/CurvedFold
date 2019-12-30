/*************************************************************
	Curvature.h
	PCAによる主軸方向の設定
	2次曲面a*x^2+b*y^2+h*x*y=0に近似
	※NPAL卒論プログラムから流用
*************************************************************/

#ifndef CURVATURE_H
#define CURVATURE_H

//#include "PC3D.h"

typedef struct{
	double m0[16];		// 固有ベクトル＋位置（回転前）
	double ev[3];		// 固有値
	double a,b,h;		// 二次曲面の係数，z = a*x*x + b*y*y + 2*h*x*y
	double sint,cost;	// pivot角度
	double m[16];		// 固有ベクトル＋位置
	double A,B;			// 二次曲面の係数，Z = A*X*X + B*Y*Y
	double err;			// 二次曲面近似誤差 090514追加
	double rho0,rho1;	// 曲率
	double shapeindex;	
} CurvInfo;

/*		３次元の点群 → 主成分分析
入力	double cen[3];			// 中心座標
		double x[], y[], z[];	// 点群座標
		int n;					// 点数
出力	double eval[3];			// 固有値
		double evec[4*4];		// 固有ベクトル，中心座標を含む変換行列
返り値	0:正常終了，1:異常終了
*/
extern int PCA3(double *x, double *y, double *z, int n, double *_cen,
				double *eval, double *evec);


/* 二次曲面 a*x^2 + b*y^2 + 2h*x*y - z = 0 に近似する
入力	double cen[3];			// 中心座標
		double x[], y[], z[];	// 点群座標
		int n;					// 点数
出力	double a,b,h;			// 二次曲面係数
返り値	0:正常終了，1:異常終了
*/
extern int QuadricSurface(double *px, double *py, double *pz, int n,
						  double *a, double *b, double *h);

extern double calcError(double *px, double *py, double *pz, int n, double a, double b, double h);

extern int Curvature(double *cen, double *nrm, double *_px, double *_py, double *_pz, int n,
			  CurvInfo *cinfo);

extern void dumpCurvInfo(CurvInfo *info);
extern void dumpCurvInfo(FILE *fp, CurvInfo *info);
extern void readCurvInfo(FILE *fp, CurvInfo *info);


extern int PCA2(double *x, double *y, int n, double *cen, double *eval, double *evec);


#endif // CURVATURE_H
