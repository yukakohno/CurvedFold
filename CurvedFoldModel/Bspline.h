#ifndef BSPLINE
#define BSPLINE

class Bspline
{
public:
	Bspline() {}
	~Bspline() {}

	// toption 0:“™ŠÔŠu, 1:“_ŠÔ‹——£‚É‰‚¶‚Ä
	int calcT(int tcnt, double tmax, int toption, double *t, double *vx=NULL, double *vy=NULL, double *vz=NULL);

	double baseN(int i,int k,double t,int *nv);

	int calcBSpline(double *cx, double *cy, double *cz, int ccnt,
		int *nv, int nvcnt, int deg, int pcnt, double *t,
		double *vx, double *vy, double *vz );
	int calcCP( double *vx, double *vy, double *vz, int pcnt,
		int *nv, int nvcnt, int deg, double *t, int ccnt,
		double *cx, double *cy, double *cz );

	// ”÷•ª’lZo
	int calcBSDiff(double *cx, double *cy, double *cz, int ccnt,
		int *nv, int nvcnt, int deg, int pcnt, double *t,
		double *vx, double *vy, double *vz, int diffn=0 ); // diffn ”÷•ª‰ñ”

#if 0
	int calcBSDiff2(double *cx, double *cy, double *cz, int ccnt,
		int *nv, int nvcnt, int deg, int pcnt, double *t,
		double *vx, double *vy, double *vz );
#endif
};

#endif
