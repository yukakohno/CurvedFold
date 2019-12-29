// 簡略化した３次スプライン曲線の生成方法
// http://www5d.biglobe.ne.jp/stssk/maze/spline.html

#define MaxSplineSize 100

class Spline {
	int num;
	double a[MaxSplineSize+1], b[MaxSplineSize+1], c[MaxSplineSize+1], d[MaxSplineSize+1];
public:
	Spline() { num=0; }
	void init(double *sp, int num);
	double culc(double t);
	double culc_d(double t);
	double culc_d2(double t);
	double culc_d3(double t);
};
