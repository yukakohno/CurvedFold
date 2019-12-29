#ifndef BEZIER
#define BEZIER

#define BTCNT 100

class Bezier
{
public:
	Bezier() {}
	~Bezier() {}

	double tt[BTCNT+1], tx[BTCNT+1], ty[BTCNT+1];

	double x0,y0, x1,y1, x2,y2, x3,y3;

	int set( double x0, double y0, double dx0, double dy0, double len0,
		double x1, double y1, double dx1, double dy1, double len1,
		char *fname );

	double gety( double x );

};

#endif
