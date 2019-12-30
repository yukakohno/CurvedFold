/*************************************************************
	Curvature.cpp
	PCAによる主軸方向の設定
	2次曲面a*x^2+b*y^2+h*x*y=0に近似
	※NPAL井芹卒論プログラムから流用
*************************************************************/

// 2002 12/20
// 固有ベクトルを求める
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "Curvature.h"
#include "util.h"

#define PI 3.1415

// 単位ベクトルに入ってきたベクトルを変換する
double unitvec( double v[] )
{
   double d, eps=0.00001;
   if( (d=sqrt((v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))) < eps ) 
	   return 0.;
   v[0] /= d;  v[1] /= d;  v[2] /= d;  
   return d;
}
// ベクトルの大きさを計算する
double veclen( double v[] ){
	return( (double)sqrt((double)(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])) ); 
}
// 二つのベクトル間の角度を計算する
double vecangle( double u[], double v[] )
{
      double a,d1,d2;
      d1=veclen(u); d2=veclen(v);
      a = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
      return( 180.*(double)acos((double)(a/(d1*d2)))/3.141593 );
}
// 外積を計算する
// 引数３つ。u1[]とu2[]の外積を、u3[]に格納する。
void CrossProduct(double u1[3],double u2[3],double u3[3])
{
	u3[0] = u1[1]*u2[2] - u1[2]*u2[1];
	u3[1] = u1[2]*u2[0] - u1[0]*u2[2];
	u3[2] = u1[0]*u2[1] - u1[1]*u2[0];
}
// 固有値ベクトルの順番を合わせる
// 引数は４つ、固有値と固有ベクトル３つ
void SortLum(double lu[3],double ve1[3],double ve2[3],double ve3[3])
{
	double temp=0.0;
	double vtemp = 0.0;
	double vector[3][3];
	int i,j,k;

	// データ読み込み
	for(i = 0; i< 3; i++){
		vector[0][i] = ve1[i];
		vector[1][i] = ve2[i];
		vector[2][i] = ve3[i];
	}
	// バブルソートを使って、降順に並び替え
	for (i = 0; i < 2/*n - 1*/; i++) {
		for (j = 2/*n - 1*/; j > i; j--) {
			if (lu[j - 1] < lu[j]) {	// 前の要素の方が小さかったら 
				temp = lu[j];			// 交換する 
				lu[j] = lu[j - 1];
				lu[j - 1]= temp;
				for(k = 0; k < 3; k++){
					vtemp = vector[j][k];
					vector[j][k] = vector[j-1][k];
					vector[j-1][k] = vtemp;
				}
			}
		}	
	}
	// データ書き込み
	for(i = 0; i< 3; i++){
		ve1[i] = vector[0][i];
		ve2[i] = vector[1][i];
		ve3[i] = vector[2][i];
	}
}

//
//	対象行列 mat[3*3] から固有値 evel[3] 固有ベクトル evec[9] 算出
//
int Eigen3(double *mat, double *eval, double *evec)
{
	int  i;
	double v1[3],v2[3],v3[3];
	double a, b, c, d, e, f, g, h, j, k, ab, q, q2, q3, lum[3];
    double m, n1, n2, t;
    double r,r1,s,s1,eps=0.001;
    double x1,y1,x2,y2,x3,y3,x4,tmp,y4,x5,y5,x6,y6;

//	a	d	e
//		b	f
//			c
	a = mat[0];	b = mat[4];	c = mat[8];
	d = mat[1];	e = mat[2];	f = mat[5];

    g = a + b + c;
    h = -(a*a + b*b + c*c) + a*b + b*c + c*a - 3.*d*d - 3.*e*e - 3.*f*f;
    j = a*b + b*c + c*a - d*d - e*e - f*f;
    k = - a*b*c + c*d*d +b*e*e - 2.*d*e*f + a*f*f;
    m = 2.*g*g*g - 9.*g*j -27.*k;
    n1 = 4.*h*h*h + m*m;

    if( n1 >= 0. ){
		return -1;
	}
    n2 = sqrt(-n1);

    ab = pow( sqrt(m*m - n1),(1./3.) );  

	// 090616 ykohno追加
	if(!ab || !e){	// 分母
		return -1;
	}

	t = atan2(n2,m)/3.;
    r = ab*cos(t); 
	s = ab*sin(t);
    ab = 1./ab;  r1 = ab*cos(t);  s1 = -ab*sin(t);
    q = pow(2.,(1./3.));
	q2 = pow(2.,(2./3.));
	q3 = sqrt(3.);
    lum[0] = (2.*g + q2*r - 2.*q*h*r1)/6.;
    lum[1] = (4.*g - q*(q*r - 2.*h*r1 + q*q3*s + 2.*q3*h*s1))/12.;
    lum[2] = (4.*g + q*(-q*r + 2.*h*r1 + q*q3*s + 2.*q3*h*s1))/12.;
    
	x1 = b*e - d*f - e*(2.*g + q2*r - 2.*q*h*r1)/6.;
    y1 = -e*(q2*s - 2.*q*h*s1)/6.;
    x2 = c +(-2.*g - q2*r + 2.*q*h*r1)/6.;
    y2 = (-q2*s + 2.*q*h*s1)/6.;
    ab = x1*x1 + y1*y1;
    v1[0]= (f*(e*f*x1 - d*(x1*x2 + y1*y2)) - x2*ab)/(e*ab);
    v1[1]= (-e*f*x1 + d*x1*x2 + d*y1*y2)/ab;
    v1[2]= 1.;
    
	tmp = q*(q*r - 2.*h*r1 + q*q3*s + 2.*q3*h*s1)/12.;
    x3 = b - g/3. + tmp;
    y3 = (q*q3*r - q*s + 2.*h*(q3*r1 + s1))/(-6.*q2);
    x4 = c - g/3. + tmp;
    y4 = y3;
    ab = (e*x3 - d*f)*(e*x3 - d*f) + e*e*y3*y3;
    v2[0]= (f*(-e*e*f*x3 - d*d*f*x4 + d*e*(f*f + x3*x4 + y3*y4))
            + x4*ab) / (-e*ab);
    v2[1]= (-e*e*f*x3 -d*d*f*x4 + d*e*(f*f + x3*x4 + y3*y4))/ab;
    v2[2]= 1.;
    
	tmp = q*(q*r - 2.*h*r1 - q*q3*s - 2.*q3*h*s1)/12.;
    x5 = b - g/3. + tmp;
    y5 = (q*q3*r + 2.*q3*h*r1 + q*s - 2.*h*s1)/(6.*q2);
    x6 = c - g/3. + tmp;
    y6 = y5;
    ab = (e*x5 - d*f)*(e*x5 - d*f) + e*e*y5*y5;
    v3[0]= (f*(-e*e*f*x5 - d*d*f*x6 + d*e*(f*f + x5*x6 + y5*y6)) + x6*ab)/(-e*ab);
    v3[1]= (-e*e*f*x5 - d*d*f*x6 + d*e*(f*f + x5*x6 + y5*y6))/ab;
    v3[2]= 1.;
    
	// 単位ベクトルにする
	unitvec(v1); unitvec(v2); unitvec(v3);
	// 固有値の大きい順位、ベクトルを並べる
	SortLum(lum,v1,v2,v3);
	// 外積を計算する
	CrossProduct(v1,v2,v3);
#if 0
	// 固有ベクトルにする
	for( i=0; i<3; i++ ){
		v1[i] *= (double)sqrt((double)lum[0]);
		v2[i] *= (double)sqrt((double)lum[1]);
		v3[i] *= (double)sqrt((double)lum[2]);
    }
#endif

	// 固有値、固有ベクトルを格納する
	if(eval){
		eval[0] = lum[0];	eval[1] = lum[1];	eval[2] = lum[2];
	}
	if(evec){
		evec[0] = v1[0];	evec[1] = v1[1];	evec[2] = v1[2];
		evec[3] = v2[0];	evec[4] = v2[1];	evec[5] = v2[2];
		evec[6] = v3[0];	evec[7] = v3[1];	evec[8] = v3[2];
	}
	return 0;
}


/*		３次元の点群 → 主成分分析
入力	double cen[3];			// 中心座標
		double x[], y[], z[];	// 点群座標
		int n;					// 点数
出力	double eval[3];			// 固有値
		double evec[4*4];		// 固有ベクトル，中心座標を含む変換行列
返り値	0:正常終了，1:異常終了
*/
int PCA3(double *x, double *y, double *z, int n, double *_cen,
			double *eval, double *evec)
{
	int  i;
	double cen[3],v1[3],v2[3],v3[3];
	double a, b, c, d, e, f, g, h, j, k, ab, q, q2, q3, lum[3];
    double m, n1, n2, t;
    double r,r1,s,s1,eps=0.001;
    double x1,y1,x2,y2,x3,y3,x4,tmp,y4,x5,y5,x6,y6;

	if(!_cen){
		// 中心が設定されていない場合平均位置を中心とする
		cen[0] = cen[1] = cen[2] = 0;
		for(i = 0; i < n; i++){
			cen[0] += x[i];	cen[1] += y[i];	cen[2] += z[i];
		}
		cen[0] /= n;	cen[1] /= n;	cen[2] /= n;
	} else {
		cen[0] = _cen[0];	cen[1] = _cen[1];	cen[2] = _cen[2];
	}

    a = b = c = d = e = f = 0.;
    for( i=0; i < n; i++ ){
        a += (x[i] - cen[0])*(x[i] - cen[0]);
        b += (y[i] - cen[1])*(y[i] - cen[1]);
        c += (z[i] - cen[2])*(z[i] - cen[2]);
        d += (x[i] - cen[0])*(y[i] - cen[1]);
        f += (y[i] - cen[1])*(z[i] - cen[2]);
        e += (z[i] - cen[2])*(x[i] - cen[0]);

    }
    a /= (double)n; b /= (double)n; c /= (double)n;
    d /= (double)n; e /= (double)n; f /= (double)n;

    g = a + b + c;
    h = -(a*a + b*b + c*c) + a*b + b*c + c*a - 3.*d*d - 3.*e*e - 3.*f*f;
    j = a*b + b*c + c*a - d*d - e*e - f*f;
    k = - a*b*c + c*d*d +b*e*e - 2.*d*e*f + a*f*f;
    m = 2.*g*g*g - 9.*g*j -27.*k;
    n1 = 4.*h*h*h + m*m;

    if( n1 >= 0. ){
		return -1;
	}
    n2 = sqrt(-n1);

    ab = pow( sqrt(m*m - n1),(1./3.) );  

	// 090616 ykohno追加
	if(!ab || !e){	// 分母
		return -1;
	}

	t = atan2(n2,m)/3.;
    r = ab*cos(t); 
	s = ab*sin(t);
    ab = 1./ab;  r1 = ab*cos(t);  s1 = -ab*sin(t);
    q = pow(2.,(1./3.));
	q2 = pow(2.,(2./3.));
	q3 = sqrt(3.);
    lum[0] = (2.*g + q2*r - 2.*q*h*r1)/6.;
    lum[1] = (4.*g - q*(q*r - 2.*h*r1 + q*q3*s + 2.*q3*h*s1))/12.;
    lum[2] = (4.*g + q*(-q*r + 2.*h*r1 + q*q3*s + 2.*q3*h*s1))/12.;
    
	x1 = b*e - d*f - e*(2.*g + q2*r - 2.*q*h*r1)/6.;
    y1 = -e*(q2*s - 2.*q*h*s1)/6.;
    x2 = c +(-2.*g - q2*r + 2.*q*h*r1)/6.;
    y2 = (-q2*s + 2.*q*h*s1)/6.;
    ab = x1*x1 + y1*y1;
    v1[0]= (f*(e*f*x1 - d*(x1*x2 + y1*y2)) - x2*ab)/(e*ab);
    v1[1]= (-e*f*x1 + d*x1*x2 + d*y1*y2)/ab;
    v1[2]= 1.;
    
	tmp = q*(q*r - 2.*h*r1 + q*q3*s + 2.*q3*h*s1)/12.;
    x3 = b - g/3. + tmp;
    y3 = (q*q3*r - q*s + 2.*h*(q3*r1 + s1))/(-6.*q2);
    x4 = c - g/3. + tmp;
    y4 = y3;
    ab = (e*x3 - d*f)*(e*x3 - d*f) + e*e*y3*y3;
    v2[0]= (f*(-e*e*f*x3 - d*d*f*x4 + d*e*(f*f + x3*x4 + y3*y4))
            + x4*ab) / (-e*ab);
    v2[1]= (-e*e*f*x3 -d*d*f*x4 + d*e*(f*f + x3*x4 + y3*y4))/ab;
    v2[2]= 1.;
    
	tmp = q*(q*r - 2.*h*r1 - q*q3*s - 2.*q3*h*s1)/12.;
    x5 = b - g/3. + tmp;
    y5 = (q*q3*r + 2.*q3*h*r1 + q*s - 2.*h*s1)/(6.*q2);
    x6 = c - g/3. + tmp;
    y6 = y5;
    ab = (e*x5 - d*f)*(e*x5 - d*f) + e*e*y5*y5;
    v3[0]= (f*(-e*e*f*x5 - d*d*f*x6 + d*e*(f*f + x5*x6 + y5*y6)) + x6*ab)/(-e*ab);
    v3[1]= (-e*e*f*x5 - d*d*f*x6 + d*e*(f*f + x5*x6 + y5*y6))/ab;
    v3[2]= 1.;
    
	// 単位ベクトルにする
	unitvec(v1); unitvec(v2); unitvec(v3);
	// 固有値の大きい順位、ベクトルを並べる
	SortLum(lum,v1,v2,v3);
	// 外積を計算する
	CrossProduct(v1,v2,v3);
#if 0
	// 固有ベクトルにする
	for( i=0; i<3; i++ ){
		v1[i] *= (double)sqrt((double)lum[0]);
		v2[i] *= (double)sqrt((double)lum[1]);
		v3[i] *= (double)sqrt((double)lum[2]);
    }
#endif

	// 固有値、固有ベクトルを格納する
	if(eval){
		eval[0] = lum[0];	eval[1] = lum[1];	eval[2] = lum[2];
	}
	if(evec){
		evec[0] = v1[0];	evec[1] = v1[1];	evec[2] = v1[2];	evec[3] = 0;
		evec[4] = v2[0];	evec[5] = v2[1];	evec[6] = v2[2];	evec[7] = 0;
		evec[8] = v3[0];	evec[9] = v3[1];	evec[10] = v3[2];	evec[11] = 0;
		evec[12] = cen[0];	evec[13] = cen[1];	evec[14] = cen[2];	evec[15] = 1;
	}
	return 0;
}

/* Cramerの公式
入力	z[3]
		c[3*3]
出力	x[3]
返り値	0:正常終了，1:異常終了
*/
int Cramer(double *z,double *c, double *x)
{
	int i;
	double Total;
	double copy[9];
	double e = 0.1;
	
	/*	0	1	2
		3	4	5
		6	7	8	*/

	// 1/|A|の計算
	Total = c[0]*c[4]*c[8] + c[1]*c[5]*c[6] + c[2]*c[3]*c[7]
			- c[2]*c[4]*c[6] - c[1]*c[3]*c[8] - c[0]*c[5]*c[7];

	if(Total < e){
		return -1;
	}
	// 計算のための行列の初期化
	memcpy(copy,c,sizeof(double)*9);

	// x[i]の計算
	for(i = 0; i < 3; i++){
		// 行列を解計算のために変える
		copy[0+i] = z[0];	copy[3+i] = z[1];	copy[6+i] = z[2];

		x[i] = copy[0]*copy[4]*copy[8] + copy[1]*copy[5]*copy[6] + copy[2]*copy[3]*copy[7]
			- copy[2]*copy[4]*copy[6] - copy[1]*copy[3]*copy[8] - copy[0]*copy[5]*copy[7];

		// 行列を基に戻す
		copy[0+i] = c[0+i];	copy[3+i] = c[3+i];	copy[6+i] = c[6+i];
	}

	x[0] /= Total;
	x[1] /= Total;
	x[2] /= Total;
	return 0;
}

/* 二次曲面 a*x^2 + b*y^2 + 2h*x*y - z = 0 に近似する
入力	double cen[3];			// 中心座標
		double x[], y[], z[];	// 点群座標
		int n;					// 点数
出力	double a,b,h;			// 二次曲面係数
返り値	0:正常終了，1:異常終了
*/
int QuadricSurface(double *px, double *py, double *pz, int n,
			double *a, double *b, double *h)
{
	int i;
	double x,y,z, x4,x3y1,x2y2,x1y3,y4,x2z1,y2z1,xyz;
	double C[9], Z[3], X[3];

	memset(C, 0, sizeof(double)*9);
	memset(Z, 0, sizeof(double)*3);

	x4 = x3y1 = x2y2 = x1y3 = y4 = x2z1 = y2z1 = xyz = 0.0;
	for(i = 0; i < n; i++){
		// 点の座標
		x = px[i];	y = py[i];	z = pz[i];

		x4 = x4 + x*x*x*x;
		x3y1 = x3y1 + x*x*x*y;
		x2y2 = x2y2 + x*x*y*y;
		x1y3 = x1y3 + x*y*y*y;
		y4 = y4 + y*y*y*y;

		x2z1 = x2z1 + x*x*z;
		y2z1 = y2z1 + y*y*z;
		xyz = xyz + x*y*z;
	}
	// 行列Ｃ[3][3]
	C[0] = x4;		C[1] = 2*x3y1; 	C[2] = x2y2;
	C[3] = x3y1;	C[4] = 2*x2y2;	C[5] = x1y3;
	C[6] = x2y2;	C[7] = 2*x1y3;	C[8] = y4;
	// 行列Ｚ[3]
	Z[0] = x2z1;	Z[1] = xyz;		Z[2] = y2z1;

	// クラメルの公式を使って計算する
	if( Cramer(Z,C,X) < 0 ){
		return -1;
	}

	*a = X[0];
	*h = X[1];
	*b = X[2];

	return 0;
}

double calcError(double *px, double *py, double *pz, int n, double a, double b, double h)
{
	int i;
	double tmperr, err = 0.0;

	for(i = 0; i < n; i++){
		tmperr = (a*px[i]*px[i] + b*py[i]*py[i] + 2*h*px[i]*py[i]) - pz[i];
		err += tmperr * tmperr;
	}
	return err/(double)n;
}

int Curvature(double *cen, double *nrm, double *_px, double *_py, double *_pz, int n,
			  CurvInfo *cinfo)
{
/*	// 法線方向算出
入力	double cen[3];			// 中心座標，NULLの場合は重心
		double nrm[3];			// 法線方向，NULLの場合はPCAの第3固有ベクトル
		double x[], y[], z[];	// 点群座標
		int n;					// 点数
出力	info					// 曲率情報
返り値	0:正常終了，1:異常終了
*/
	int i, sts = 0;
	double *px, *py, *pz, eval[3],evec[16];
	double tmpx0,tmpy0,tmpz0,tmpx1,tmpy1,tmpz1;
	double a,b,h,root,err;

	px = new double[n];	memcpy(px,_px,sizeof(double)*n);
	py = new double[n];	memcpy(py,_py,sizeof(double)*n);
	pz = new double[n];	memcpy(pz,_pz,sizeof(double)*n);

	// 座標変換
	if(nrm && (nrm[0] || nrm[1] || nrm[2])){	// 法線方向設定済
		// 指定した法線方向使用

		// 固有値は使わないので暫定0
		eval[0] = eval[1] = eval[2] = 0;
		evec[3] = evec[7] = evec[11] = 0; evec[15] = 1;

		// 中心設定されてない場合は重心で
		if(cen){
			evec[12] = cen[0];
			evec[13] = cen[1];
			evec[14] = cen[2];
		} else{
			double c[3] = { 0.0, 0.0, 0.0 };
			for(i = 0; i < n; i++){
				c[0] += px[i];	c[1] += py[i];	c[2] += pz[i];
			}
			evec[12] = c[0]/n;
			evec[13] = c[1]/n;
			evec[14] = c[2]/n;
		}

		evec[8] = nrm[0];	evec[9] = nrm[1];	evec[10] = nrm[2];	// z軸
		// (0,1,0)×z軸 -> x軸
//		evec[0] = evec[10];	evec[1] = 0;	evec[2] = -evec[8];
		double tmpvec[3] = {0,1,0};
		CrossProduct(tmpvec, &(evec[8]), &(evec[0]));
		// z軸×x軸 -> y軸
//		evec[4] = evec[9] * evec[2] - evec[10] * evec[1];
//		evec[5] = evec[10] * evec[0] - evec[8] * evec[2];
//		evec[6] = evec[8] * evec[1] - evec[9] * evec[0];
		CrossProduct(&(evec[8]), &(evec[0]), &(evec[4]));

	} else{
		// PCAの法線方向使用
		sts = PCA3(px, py, pz, n, cen, eval, evec);
		if(sts < 0) goto end;
	}

	// 途中結果格納
	memcpy(cinfo->ev,eval,sizeof(double)*3);
	memcpy(cinfo->m0,evec,sizeof(double)*16);

	for(i = 0; i < n; i++){
		tmpx0 = px[i]-evec[12];
		tmpy0 = py[i]-evec[13];
		tmpz0 = pz[i]-evec[14];
		tmpx1 = evec[0]*tmpx0 + evec[1]*tmpy0 + evec[2]*tmpz0;
		tmpy1 = evec[4]*tmpx0 + evec[5]*tmpy0 + evec[6]*tmpz0;
		tmpz1 = evec[8]*tmpx0 + evec[9]*tmpy0 + evec[10]*tmpz0;
		px[i] = tmpx1;
		py[i] = tmpy1;
		pz[i] = tmpz1;
	}

	// 二次曲面への近似
	sts = QuadricSurface(px, py, pz, n, &a, &b, &h);
	err = calcError(px, py, pz, n, a, b, h);
	root = sqrt((a-b)*(a-b)+4*h*h);

	// find pivot angle [-PI/2, PI/2]
	double sin2t, cos2t, sint, cost;
	cos2t = (b-a)/root;	// [-1,1]
	sin2t = 2*h/root;	// [-1,1]
	if( sin2t >= 0.0 ) {
		// 半角の公式
		sint = sqrt((1.0-cos2t)*0.5);	// [0, PI/2]
		cost = sqrt((1.0+cos2t)*0.5);	// [0, PI/2]
	} else{
		sint = -sqrt((1.0-cos2t)*0.5);	// [-PI/2, 0]
		cost = sqrt((1.0+cos2t)*0.5);	// [0, PI/2]
	}

	// pivot point cloud, 座標変換してZ=AX^2+BY^2の形に
	for(i = 0; i < n; i++){
		tmpx0 = px[i]*cost - py[i]*sint;
		tmpy0 = px[i]*sint + py[i]*cost;
		px[i] = tmpx0;
		py[i] = tmpy0;
	}

	// pivot axis
	double pivot[16] = {cost,-sint,0,0, sint,cost,0,0, 0,0,1,0, 0,0,0,1};
//	double pivot[16] = {cost,sint,0,0, -sint,cost,0,0, 0,0,1,0, 0,0,0,1};
	mult_m44_n44_rot(pivot, evec);	// B=A*B, _rot:回転成分のみかける

#if 0
	// DEBUG: pivot前後で点群の回転方向，軸の回転方向が合っているか確認
	double *px0, *py0, *pz0;
	// backup p* -> p*0
	px0 = new double[n];	memcpy(px0,px,sizeof(double)*n);
	py0 = new double[n];	memcpy(py0,py,sizeof(double)*n);
	pz0 = new double[n];	memcpy(pz0,pz,sizeof(double)*n);
	// copy again _p* -> p*
	memcpy(px,_px,sizeof(double)*n);
	memcpy(py,_py,sizeof(double)*n);
	memcpy(pz,_pz,sizeof(double)*n);
	// 座標変換 p*
	for(i = 0; i < n; i++){
		tmpx0 = px[i]-evec[12];
		tmpy0 = py[i]-evec[13];
		tmpz0 = pz[i]-evec[14];
		tmpx1 = evec[0]*tmpx0 + evec[1]*tmpy0 + evec[2]*tmpz0;
		tmpy1 = evec[4]*tmpx0 + evec[5]*tmpy0 + evec[6]*tmpz0;
		tmpz1 = evec[8]*tmpx0 + evec[9]*tmpy0 + evec[10]*tmpz0;
		px[i] = tmpx1;
		py[i] = tmpy1;
		pz[i] = tmpz1;
	}
	for(i = 0; i < n; i++){
		if(px0[i] != px[i])
			printf("%d : inconsistent data x\n", i);
		if(py0[i] != py[i])
			printf("%d : inconsistent data y\n", i);
		if(pz0[i] != pz[i])
			printf("%d : inconsistent data z\n", i);
	}
	delete[] px0; px0 = NULL;
	delete[] py0; py0 = NULL;
	delete[] pz0; pz0 = NULL;
#endif

	// 結果格納
	memcpy(cinfo->m,evec,sizeof(double)*16);
	cinfo->a = a;
	cinfo->b = b;
	cinfo->h = h;
	cinfo->sint = sint;
	cinfo->cost = cost;
	cinfo->rho0 = (a+b) - root; // 2*A
	cinfo->rho1 = (a+b) + root; // 2*B
	cinfo->A = cinfo->rho0 * 0.5;
	cinfo->B = cinfo->rho1 * 0.5;
	cinfo->err = err;
	cinfo->shapeindex = 2.0/PI*atan(-(a+b)/root);

#if 1 // debug
	memcpy(_px,px,sizeof(double)*n);
	memcpy(_py,py,sizeof(double)*n);
	memcpy(_pz,pz,sizeof(double)*n);
#endif

end:
	if(px){ delete[] px; px = NULL; }
	if(py){ delete[] py; py = NULL; }
	if(pz){ delete[] pz; pz = NULL; }
	return sts;
}

void dumpCurvInfo(CurvInfo *cinfo)
{
	printf("X(%.2f,%.2f,%.2f), Y(%.2f,%.2f,%.2f), Z(%.2f,%.2f,%.2f)\n",
		cinfo->m0[0],cinfo->m0[1],cinfo->m0[2], cinfo->m0[4],cinfo->m0[5],cinfo->m0[6], cinfo->m0[8],cinfo->m0[9],cinfo->m0[10]);
	printf("Z = %f x^2 + %f y^2 + %f xy\n", cinfo->a, cinfo->b, cinfo->h*2);
	printf("pivot angle= %f or %f, sint= %f, cost= %f\n",
		asin(cinfo->sint)*180.0/PI, acos(cinfo->cost)*180.0/PI, cinfo->sint, cinfo->cost);
	printf("X(%.2f,%.2f,%.2f), Y(%.2f,%.2f,%.2f), Z(%.2f,%.2f,%.2f)\n",
		cinfo->m[0],cinfo->m[1],cinfo->m[2], cinfo->m[4],cinfo->m[5],cinfo->m[6], cinfo->m[8],cinfo->m[9],cinfo->m[10]);
	printf("Z = %f X^2 + %f Y^2\n", cinfo->A, cinfo->B);
	printf("curvature: %f,%f\tshape index: %f\n", cinfo->rho0, cinfo->rho1, cinfo->shapeindex);
	printf("gause: %f\taverage: %f\n", cinfo->rho0*cinfo->rho1, (cinfo->rho0+cinfo->rho1)*0.5 );
}

void dumpCurvInfo(FILE *fp, CurvInfo *cinfo)
{
/*
// 曲率情報構造体
// 出力するのは一部のみ
typedef struct{
	double m0[16];		// 固有ベクトル＋位置（回転前）
	double ev[3];		// 固有値
	double a,b,h;		// 二次曲面の係数，z = a*x*x + b*y*y + 2*h*x*y
	double sint,cost;	// pivot角度
	double m[16];		// 固有ベクトル＋位置
	double A,B;			// 二次曲面の係数，Z = A*X*X + B*Y*Y
	double rho0,rho1;	// 曲率
	double shapeindex;	
} CurvInfo;
*/
	fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,",
		cinfo->m[0],cinfo->m[1],cinfo->m[2], cinfo->m[4],cinfo->m[5],cinfo->m[6], cinfo->m[8],cinfo->m[9],cinfo->m[10]);
	fprintf(fp,"%f,%f,%f,%f,%f\n", cinfo->A, cinfo->B, cinfo->rho0, cinfo->rho1, cinfo->shapeindex);
}

void readCurvInfo(FILE *fp, CurvInfo *cinfo)
{
	char buf[1024];

	if(fgets(buf,1024,fp)){
		sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", 
			&cinfo->m[0],&cinfo->m[1],&cinfo->m[2],
			&cinfo->m[4],&cinfo->m[5],&cinfo->m[6],
			&cinfo->m[8],&cinfo->m[9],&cinfo->m[10],
			&cinfo->A, &cinfo->B, &cinfo->rho0, &cinfo->rho1, &cinfo->shapeindex);
	}
}


int PCA2(double *x, double *y, int n, double *_cen,
		 double *eval, double *evec)
{
	int  i;
	double cen[2], eps=0.0000001, a,b,c,D, len, lum0,lum1, v00,v01,v10,v11;

	if(!_cen){
		// 中心が設定されていない場合平均位置を中心とする
		cen[0] = cen[1] =0;
		for(i = 0; i < n; i++){
			cen[0] += x[i];	cen[1] += y[i];
		}
		cen[0] /= n;	cen[1] /= n;
	} else {
		cen[0] = _cen[0];	cen[1] = _cen[1];
	}

    a = b = c = 0.;
    for( i = 0; i < n; i++ ){
        a += (x[i] - cen[0])*(x[i] - cen[0]);
        b += (y[i] - cen[1])*(y[i] - cen[1]);
        c += (x[i] - cen[0])*(y[i] - cen[1]);
    }
    a /= (double)n; b /= (double)n; c /= (double)n;

	// {{a,c}{c,b}} の固有値・固有ベクトル -> eval,evec
	D = (a-b)*(a-b)+4*c*c;
	if(D < eps){
		goto ERR;
	}
	D = sqrt(D);
	lum0 = (a+b+D)/2.0;
	lum1 = (a+b-D)/2.0;

	v00 = b-c - lum0;	v01 = a-c - lum0;
	len = sqrt(v00*v00+v01*v01);
	v00 /= len;	v01 /= len;

	v10 = b-c - lum1;	v11 = a-c - lum1;
	len = sqrt(v10*v10+v11*v11);
	v10 /= len;	v11 /= len;

	// 固有値、固有ベクトルを格納する
	if(eval){
		eval[0] = lum0;	eval[1] = lum1;
	}
	if(evec){
		evec[0] = v00;		evec[1] = v01;		evec[2] = 0;
		evec[3] = v10;		evec[4] = v11;		evec[5] = 0;
		evec[6] = cen[0];	evec[7] = cen[1];	evec[8] = 1;
	}
	return 0;
ERR:
	return -1;
}
