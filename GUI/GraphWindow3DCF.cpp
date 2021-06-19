#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

//#include <GL/glu.h>
#include <GL/glut.h>
//#include <glut.h>
#include <FL/Fl.h>
#include <FL/gl.h>
#include <FL/Enumerations.H>
#include <FL/fl_draw.H>
#include "GraphWindow3DCF.h"
#include "ControlPanel.h"
#include "../CurvedFoldModel/util.h"

/* OpenCV */
#include <opencv2/opencv.hpp>
//using namespace std;
//using namespace cv;

GraphWindow3DCF::GraphWindow3DCF(int X, int Y, int W, int H, const char *L)
: GraphWindow3D(X, Y, W, H, L)
{
	GraphWindow3DCF(X, Y, W, H);
}

GraphWindow3DCF::GraphWindow3DCF(int X, int Y, int W, int H)
: GraphWindow3D(X, Y, W, H)
{
	init();
	sbw = this->w();
	sbh = this->h();
	sbstep = sbw * 3;
	selectbuffer = new unsigned char[sbstep * sbh];
}

GraphWindow3DCF::~GraphWindow3DCF()
{
	init();

	delete[] selectbuffer;
}

void GraphWindow3DCF::initObject()
{
	prevquat2[0] = prevquat2[1] = prevquat2[2] = 0; prevquat2[3] = 1;
	currquat2[0] = currquat2[1] = currquat2[2] = 0; currquat2[3] = 1;
	dispquat2[0] = dispquat2[1] = dispquat2[2] = 0; dispquat2[3] = 1;
	prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
	currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
	disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
	unit_m44( mObject );
}

void GraphWindow3DCF::initTexture()
{
	/* テクスチャの読み込みに使う配列 */
	GLubyte texture[TEXHEIGHT][TEXWIDTH][4];

	/* テクスチャ画像の読み込み */
#if 1
	cv::Mat tex = cv::imread(TEXFNAME, 1);
	for(int j=0 ; j<tex.rows; j++){
		for(int i=0 ; i<tex.cols; i++){
			texture[j][i][0] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[0];
			texture[j][i][1] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[1];
			texture[j][i][2] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[2];
			texture[j][i][3] = (GLubyte)255;
		}
	}
#elif 0
	IplImage *tex = cvLoadImage( "texture/curve1_cc2.png", CV_LOAD_IMAGE_COLOR );
	if( tex ){
		for(int j=0 ; j<TEXHEIGHT; j++){
			for(int i=0 ; i<TEXWIDTH; i++){
				texture[j][i][2] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep]);
				texture[j][i][1] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep +1]);
				texture[j][i][0] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep +2]);
				int thres = 20;
				if( abs(texture[j][i][0]-texture[j][i][1]) < thres
					&& abs(texture[j][i][1]-texture[j][i][2]) < thres
					&& abs(texture[j][i][2]-texture[j][i][0]) < thres )
				{
					texture[j][i][3] = (GLubyte)63;
				} else {
					texture[j][i][3] = (GLubyte)255;
				}
			}
		}
	}
	cvReleaseImage(&tex);
#else
	for(int j=0 ; j<TEXHEIGHT; j++){
		for(int i=0 ; i<TEXWIDTH; i++){
			if( (j/16)%2 == (i/16)%2 ){
			//if( (j/(TEXHEIGHT/2))%2 == (i/(TEXWIDTH/2))%2 ){
				texture[j][i][0] = (GLubyte)214;
				texture[j][i][1] = (GLubyte)214;
				texture[j][i][2] = (GLubyte)214;
				texture[j][i][3] = (GLubyte)255;
			} else {
				texture[j][i][0] = (GLubyte)239;
				texture[j][i][1] = (GLubyte)239;
				texture[j][i][2] = (GLubyte)239;
				texture[j][i][3] = (GLubyte)255;
			}
		}
	}
#endif

	/* テクスチャ画像はバイト単位に詰め込まれている */
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	/* テクスチャを拡大・縮小する方法の指定 */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	/* テクスチャの繰り返し方法の指定 */
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	/* テクスチャ環境 */
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	/* テクスチャの割り当て */
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXWIDTH, TEXHEIGHT,
		0, GL_RGBA, GL_UNSIGNED_BYTE, texture);

	/* アルファテストの判別関数 */
	//glAlphaFunc(GL_GREATER, 0.5);
}

void GraphWindow3DCF::init()
{
	GraphWindow3D::init();

	disp_axis = 1;
	disp_axis2 = 0;
	disp_X = 1;
	disp_TNB = 0;
	disp_R = 0;
	disp_maxrlen = 0;
	disp_PLY = 1;
	disp_PTN = 1;
	disp_CP = 0;	ppos = 0;	pprm = 0;
	disp_ONE = 0;
	disp_PRI = 1;
	disp_stitch = 0;
	disp_TGT = 1;
	disp_modTGT = -2;
	disp_LIN_SMOOTH = 0;

	flg_addcurve = 0;

	divtype = 1;
	divnum = 1;

	ppm = NULL;

	push_button2 = 0;
	push_x2 = push_y2 = 0;
	initObject();

	vcnt = 0;
	vx = vy = vz = NULL;

	selectidx = -1;
	resetRotation(-1);
}

void GraphWindow3DCF::resetObjTrans()
{
	prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
	currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
	disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
	mObject[12] = disptrans2[0];
	mObject[13] = disptrans2[1];
	mObject[14] = disptrans2[2];
	//redraw();
}

void GraphWindow3DCF::resetObjRot()
{
	quat2.trackball(prevquat2, 0.0, 0.0, 0.0, 0.0);
	quat2.trackball(currquat2, 0.0, 0.0, 0.0, 0.0);
	quat2.trackball(dispquat2, 0.0, 0.0, 0.0, 0.0); 
	quat2.build_rotmatrix(mObject, dispquat2);
	mObject[12] = disptrans2[0];
	mObject[13] = disptrans2[1];
	mObject[14] = disptrans2[2];
	//redraw();
}

extern double norm( double *x, double *y, double *z );

void GraphWindow3DCF::draw3DCurveFold()
{
	int i,j;

	// 0:RED, 1:YLW, 2:GRN, 3*CYN, 4:BLU, 5:PPL, 6:DRED, 7:DGRN, 8:DBLU, 9:GRAY
	float color[10][3] = {{1,0,0},{1,1,0},{0,1,0},{0,1,1},{0,0,1},{1,0,1},{.5,0,0},{0,.5,0},{0,0,.5},{.5,.5,.5}};
	float gcolor[12][3] = {{1,0,0},{1,.4,0},{1,.7,0}, {1,1,0},{.7,1,0},{.4,1,0},
	{0,1,0},{0,1,.5},{0,.5,1}, {0,0,1},{.5,0,1},{1,0,1}};
	GLfloat mat_col[10][4] = {{1,0,0,1},{1,1,0,1},{0,1,0,1},{0,1,1,1},{0,0,1,1},{1,0,1,1},{.5,0,0,1},{0,.5,0,1},{0,0,.5,1},{.5,.5,.5,1}};
	int Pcnt = 0, Xcnt = 0;
	if (ppm == NULL) {
		return;
	}
	Pcnt = ppm->crs[0].Pcnt;
	Xcnt = ppm->crs[0].Xcnt;

	glPushMatrix();

	// 原点を見るように, カメラ位置を移動
	glTranslated(0.,0.,-n_plane-1000);

	// 視点設定
	SetRTS();

	if(disp_axis)
		drawAxis(); 

	for( int dc=0; dc<divnum; dc++ ){
		glPushMatrix();
		glRotated( (double)dc * (360.0/(double)divnum), 0.0, 0.0, 1.0 );
#if 1
		if( disp_stitch && dc==0 ){
			for( int i=0; i<ppm->stpcnt; i++ ){
				glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
				glPushMatrix();
				glTranslatef( ppm->stx0[0][i], ppm->sty0[0][i], ppm->stz0[0][i] );
				glutWireCube( 1.0 );
				glPopMatrix();
				glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
				glPushMatrix();
				glTranslatef( ppm->stx1[0][i], ppm->sty1[0][i], ppm->stz1[0][i] );
				glutWireCube( 1.0 );
				glPopMatrix();
				glColor3f( 1.0, 1.0, 0.0 ); // YELLOW
				glPushMatrix();
				glTranslatef( ppm->stxm[0][i], ppm->stym[0][i], ppm->stzm[0][i] );
				glutSolidSphere( 1.0, 8, 8 );
				glPopMatrix();
				glColor3f( 0.0, 0.0, 0.0 ); // BLACK
				glBegin(GL_LINE_STRIP);
				glVertex3f( ppm->stx0[0][i], ppm->sty0[0][i], ppm->stz0[0][i] );
				glVertex3f( ppm->stxm[0][i], ppm->stym[0][i], ppm->stzm[0][i] );
				glVertex3f( ppm->stx1[0][i], ppm->sty1[0][i], ppm->stz1[0][i] );
				glEnd();
			}
			for( int i=0; i<ppm->stpcnt; i++ ){
				glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
				glPushMatrix();
				glTranslatef( ppm->stx0[1][i], ppm->sty0[1][i], ppm->stz0[1][i] );
				glutWireCube( 1.0 );
				glPopMatrix();
				glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
				glPushMatrix();
				glTranslatef( ppm->stx1[1][i], ppm->sty1[1][i], ppm->stz1[1][i] );
				glutWireCube( 1.0 );
				glPopMatrix();
				glColor3f( 1.0, 1.0, 0.0 ); // YELLOW
				glPushMatrix();
				glTranslatef( ppm->stxm[1][i], ppm->stym[1][i], ppm->stzm[1][i] );
				glutSolidSphere( 1.0, 8, 8 );
				glPopMatrix();
				glColor3f( 0.0, 0.0, 0.0 ); // BLACK
				glBegin(GL_LINE_STRIP);
				glVertex3f( ppm->stx0[1][i], ppm->sty0[1][i], ppm->stz0[1][i] );
				glVertex3f( ppm->stxm[1][i], ppm->stym[1][i], ppm->stzm[1][i] );
				glVertex3f( ppm->stx1[1][i], ppm->sty1[1][i], ppm->stz1[1][i] );
				glEnd();
			}
		}
#endif
		glPushMatrix();
		SetRTS2();
		
		if( dc==0 && disp_axis2 )
		{
			glPushMatrix();
			crease *c = &(ppm->crs[0]);
			glMultMatrixd( c->m3 );
			glScaled( 0.2, 0.2, 0.2 );
			drawAxis();
			glPopMatrix();
		}
/*-----------------------------------------------------------------*/

		if( !disp_ONE || disp_ONE && dc==0 ){

	//
	// control points
	//
	if( disp_CP && ppos>-1 && pprm>-1 ){
		switch(pprm){	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR
			case P_CV2D:	glColor3f( 1.0, 0.0, 1.0 );	break;
			case P_CV3D:	glColor3f( 0.0, 1.0, 0.0 );	break;
			case P_TRSN:	glColor3f( 0.0, 0.0, 1.0 );	break;
			case P_FLDA:	glColor3f( 1.0, 0.0, 0.0 );	break;
			case P_RULL:	glColor3f( 0.32, 0.12, 0.46 );	break; // Dark Violet: 522076
			case P_RULR:	glColor3f( 0.9, 0.0, 0.4 );	break; // Raspberry: E30B5D
			case P_TRSN1:	glColor3f( 0.0, 0.0, 1.0 );	break;
			case P_FLDA1:	glColor3f( 1.0, 0.0, 0.0 );	break;
			case P_RULL1:	glColor3f( 0.32, 0.12, 0.46 );	break; // Dark Violet: 522076
			case P_RULR1:	glColor3f( 0.9, 0.0, 0.4 );	break; // Raspberry: E30B5D
		}
		crease *c = &(ppm->crs[0]);
		int ii = (int) ( (double)ppos/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		glPushMatrix();
		glTranslatef( c->Xx[ii], c->Xy[ii], c->Xz[ii] );
		glutSolidSphere( 3.0, 8, 8 );
		glPopMatrix();

		int prev, next;
		if( ppos==0 ){ prev=0; } else { prev=ppos-1; }
		if( ppos==Pcnt-1 ){ next=Pcnt-1; } else { next=ppos+1; }
		int i0 = (int) ( (double)prev/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		int i1 = (int) ( (double)next/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for( i=i0; i<i1; i++ ){
			glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
		}
		glEnd();

		glColor3f( 0.2, 0.2, 0.2 );
		for( int jj=0; jj<Pcnt; jj++ ){
			if( jj==ppos ){
				continue;
			}
			ii = (int) ( (double)jj/(double)(Pcnt-1) * (double)(Xcnt-1) ) ;
			glPushMatrix();
			glTranslatef( c->Xx[ii], c->Xy[ii], c->Xz[ii] );
			glutSolidSphere( 2.5, 8, 8 );
			glPopMatrix();
		}
	}

	//
	// sample points for stitch
	//
	if( disp_stitch ){
		glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
		for( int i=0; i<ppm->stpcnt; i++ ){
			glPushMatrix();
			glTranslatef( ppm->stx[0][i], ppm->sty[0][i], ppm->stz[0][i] );
			glutSolidSphere( 1.0, 8, 8 );
			glPopMatrix();
		}
		glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
		for( int i=0; i<ppm->stpcnt; i++ ){
			glPushMatrix();
			glTranslatef( ppm->stx[1][i], ppm->sty[1][i], ppm->stz[1][i] );
			glutSolidSphere( 1.0, 8, 8 );
			glPopMatrix();
		}
	}

	//
	// curve
	//
	if( disp_X && ppm ){
		glLineWidth(2.0);
		glColor3f( 0.2, 0.2, 0.2 );
		//glColor3f( 0.9, 0.7, 0.0 );
		for( j=0; j<ppm->crcnt; j++ ){
			//if( disp_PRI && j>0 ){
			//	break;
			//}
			crease *c = &(ppm->crs[j]);
			glBegin(GL_LINE_STRIP);
			if( c->Xxs!=0.0 || c->Xys!=0.0 || c->Xzs!=0.0 ){
				glVertex3f(c->Xxs, c->Xys, c->Xzs);
			}
			for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
			}
			if( c->Xxe!=0.0 || c->Xye!=0.0 || c->Xze!=0.0 ){
				glVertex3f(c->Xxe, c->Xye, c->Xze);
			}
			glEnd();
		}
	}

	//
	// TNB
	//
	if( disp_TNB && ppm ){
		glLineWidth(2.0);
		glBegin(GL_LINES);
		double len=5.0;
		for( j=0; j<ppm->crcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = &(ppm->crs[j]);
			for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
				glColor3f( 1.0, 0.0, 0.0 );
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
				glVertex3f(c->Xx[i] + c->Tx[i]*len, c->Xy[i] + c->Ty[i]*len, c->Xz[i] + c->Tz[i]*len );
				glColor3f( 0.0, 1.0, 0.0 );
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
				glVertex3f(c->Xx[i] + c->Nx[i]*len, c->Xy[i] + c->Ny[i]*len, c->Xz[i] + c->Nz[i]*len );
				glColor3f( 0.0, 0.0, 1.0 );
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
				glVertex3f(c->Xx[i] + c->Bx[i]*len, c->Xy[i] + c->By[i]*len, c->Xz[i] + c->Bz[i]*len );
			}
		}
		glEnd();
	}

	//
	// Ruling
	//
	if( disp_R && ppm ){
		glLineWidth(2.0);
		glBegin(GL_LINES);
		for( j=0; j<ppm->lcrcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = ppm->lcrs[j];
			if( j%2==0 ){
				//glColor3f( 0.6, 0.4, 0.9 ); // Lavender: B57EDC
				glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
			} else {
				glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
			}
			//glColor3f( 1.0, 1.0, 0.0 );
			for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
#ifdef EVERYOTHER
				if( i%2==0 ) continue;
#endif
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
				glVertex3f(c->Xx[i] + c->rlx[i]*c->rllen[i],
					c->Xy[i] + c->rly[i]*c->rllen[i],
					c->Xz[i] + c->rlz[i]*c->rllen[i] );
			}
		}
		for( j=0; j<ppm->rcrcnt; j++ ){
			if( disp_PRI && j>0 ){
				break;
			}
			crease *c = ppm->rcrs[j];
			if( j%2==0 ){
				glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
			} else {
				//glColor3f( 0.6, 0.4, 0.9 ); // Lavender: B57EDC
				glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
			}
			//glColor3f( 1.0, 1.0, 0.0 );
			for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
#ifdef EVERYOTHER
				if( i%2==0 ) continue;
#endif
				glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
				glVertex3f(c->Xx[i] + c->rrx[i]*c->rrlen[i],
					c->Xy[i] + c->rry[i]*c->rrlen[i],
					c->Xz[i] + c->rrz[i]*c->rrlen[i] );			}
		}
		glEnd();
	}

		} // if( !disp_ONE || disp_ONE && dc==0 )

	//
	// Polygons
	//
	//if( disp_PLY && ppm ){
	if( disp_PLY && ppm || disp_ONE && dc!=0 ){
		glLineWidth(1.0);
		glColor3f( 0.0, 0.0, 0.0 );
		for( i=0; i<ppm->plcnt; i++ ){
			if( disp_PRI && ppm->pl_cridx[i]>0 ){
				continue;
			}
			glBegin(GL_LINE_STRIP);
			for( j=0; j<ppm->plvcnt[i]; j++ ){
				glVertex3f(ppm->plx[i*4+j], ppm->ply[i*4+j], ppm->plz[i*4+j] );
			}
			glVertex3f(ppm->plx[i*4], ppm->ply[i*4], ppm->plz[i*4] );
			glEnd();
		}
	}

	if( !disp_ONE || disp_ONE && dc==0 ){

	if( disp_PTN && ppm )
	{
		//double dangle = (double)dc*2.0*M_PI/(double)divnum -15.0/180.0*M_PI; // for origami-sphere
		double dangle = (double)dc * 2.0 * M_PI / (double)divnum;
		double cosd = cos(-dangle);
		double sind = sin(-dangle);

		//glEnable(GL_ALPHA_TEST);
		glEnable(GL_NORMALIZE);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_TEXTURE_2D);	/* テクスチャマッピング開始 */
		if( disp_POLY_OFFSET ){
			glEnable(GL_POLYGON_OFFSET_FILL);
		}
		glPolygonOffset(1.0, 1.0);

#if 1	// CP上のポリゴンを座標変換して表示
		glColor3f( 1.0, 0.0, 1.0 );
		for( i=0; i<ppm->plcnt; i++ ){
			if( disp_PRI && ppm->pl_cridx[i]>0 ){
				continue;
			}
			glPushMatrix();
			glMultMatrixd( &(ppm->plmat[i*16]) );
			glBegin(GL_POLYGON);
			glNormal3f( 0.0, 0.0, -1.0 );
			for( j=0; j<ppm->plvcnt[i]; j++ ){
				//glTexCoord2d( (ppm->plx_cp[i*4+j]-ppm->psx)/(double)PPWIDTH, (ppm->ply_cp[i*4+j]-ppm->psy)/(double)PPHEIGHT );
				double x0 = ppm->plx_cp[i*4+j]-ppm->psx;
				double y0 = ppm->ply_cp[i*4+j]-ppm->psy;
				double tx = cosd*x0+sind*y0;
				double ty = -sind*x0+cosd*y0;
				glTexCoord2d( tx/(double)PPWIDTH, ty/(double)PPHEIGHT );
				//glTexCoord2d( tx*8.0/(double)PPWIDTH, ty*8.0/(double)PPHEIGHT );
				glVertex3f(ppm->plx_cp[i*4+j], ppm->ply_cp[i*4+j], 0.0 );
			}
			glEnd();
			glPopMatrix();
		}
#else	// 3Dポリゴンを表示
		glColor3f( 1.0, 0.0, 1.0 );
		for( i=0; i<ppm->plcnt; i++ ){
			if( disp_PRI && ppm->pl_cridx[i]>0 ){
				continue;
			}
			glBegin(GL_POLYGON);
			for( j=0; j<ppm->plvcnt[i]; j++ ){
				glTexCoord2d( (ppm->plx_cp[i*4+j]-ppm->psx)/(double)PPWIDTH, (ppm->ply_cp[i*4+j]-ppm->psy)/(double)PPHEIGHT );
				glVertex3f(ppm->plx[i*4+j], ppm->ply[i*4+j], ppm->plz[i*4+j] );
			}
			glTexCoord2d( (ppm->plx_cp[i*4]-ppm->psx)/(double)PPWIDTH, (ppm->ply_cp[i*4]-ppm->psy)/(double)PPHEIGHT );
			glVertex3f(ppm->plx[i*4], ppm->ply[i*4], ppm->plz[i*4] );
			glEnd();
		}
#endif
		//glDisable(GL_ALPHA_TEST);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_NORMALIZE);
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

#define DEMO_TARGET
	//
	// Target points
	//
	if( disp_TGT ){
		for( i=0; i<ppm->tgcnt; i++ )
		{
			glPushMatrix();
			glTranslatef( ppm->tgx[i], ppm->tgy[i], ppm->tgz[i] );
#ifdef DEMO_TARGET
			glColor3f(0.0, 0.0, 1.0);
			if (i == disp_modTGT) {
				glColor3f(1.0, 0.5, 0.0);
			}
			glutSolidSphere(3.0, 8, 8);
#else
			glColor3f( 0.5, 0.5, 0.5 );
			glutSolidSphere( 1.5, 8, 8 );
#endif
			glPopMatrix();

			glPushMatrix();
			glTranslatef( ppm->ogx[i], ppm->ogy[i], ppm->ogz[i] );
#ifdef DEMO_TARGET
			glColor3f(1.0, 0.0, 0.0);
			if (i == disp_modTGT) {
				glColor3f(1.0, 1.0, 0.0);
			}
			glutSolidSphere(3.0, 8, 8);
#else
			float r,g,b;
			if( ppm->tgap[i] >= 0 ){
				r = ppm->tgap[i]*0.05; r = r < 1.0 ? r : 1.0;
				g = 1.0 - r;
				b = 0.0;
			} else {
				b = -ppm->tgap[i]*0.05; b = b < 1.0 ? b : 1.0;
				g = 1.0 - b;
				r = 0.0;
			}
			glColor3f( r, g, b );
			glutSolidSphere( 1.5, 8, 8 );
#endif
			glPopMatrix();

			glLineWidth(2.0);
			glBegin(GL_LINES);
			glColor3f( 0.0, 0.0, 0.0 );
			//glColor3f( 0.5, 0.5, 0.5 );
			if (i == disp_modTGT) {
				glColor3f(1.0, 0.5, 0.0);
			}
			glVertex3f( ppm->tgx[i], ppm->tgy[i], ppm->tgz[i] );
#ifndef DEMO_TARGET
			glColor3f( r, g, b );
#endif
			if (i == disp_modTGT) {
				glColor3f(1.0, 1.0, 0.0);
			}
			glVertex3f( ppm->ogx[i], ppm->ogy[i], ppm->ogz[i] );
			glEnd();
		}
	}

#if 1
	//
	// Select target points
	//
	if (disp_TGT && -1 < disp_modTGT && disp_modTGT < ppm->tgcnt && push_x2 > -1 && push_y2 > -1 && ppm)
	{
		glDrawBuffer(GL_FRONT);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//glClearColor(1.f, 1.f, 1.f, 1.f);

		for (i = 0; i < ppm->tgcnt; i++)
		{
			glPushMatrix();
			glTranslatef(ppm->tgx[i], ppm->tgy[i], ppm->tgz[i]);
			glColor3ub((GLubyte)i, (GLubyte)i, (GLubyte)i);
			glutSolidSphere(3.0, 8, 8);
			glPopMatrix();
		}

		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, sbw, sbh, GL_RGB, GL_UNSIGNED_BYTE, selectbuffer);

		//int selectidx = selectbuffer[(sbh - 1 - push_y2) * sbstep + push_x2 * 3];
		//std::cout << "selectidx = " << selectidx << ", push_x2 = " << push_x2 << ", push_y2 = " << push_y2 << std::endl;
#if 0
		cv::Mat image(cv::Size(sbw, sbh), CV_8UC3);
		for (int j = 0; j < sbh; j++) {
			for (int i = 0; i < sbw; i++) {
				image.at<cv::Vec3b>(j, i)[0] = selectbuffer[(sbh - 1 - j) * sbstep + i * 3 + 2];
				image.at<cv::Vec3b>(j, i)[1] = selectbuffer[(sbh - 1 - j) * sbstep + i * 3 + 1];
				image.at<cv::Vec3b>(j, i)[2] = selectbuffer[(sbh - 1 - j) * sbstep + i * 3];
			}
		}
		imwrite("buffer.bmp", image);
#endif
		glDrawBuffer(GL_BACK);
	}
#endif

		} // if( !disp_ONE || disp_ONE && dc==0 )
/*-----------------------------------------------------------------*/

		glPopMatrix();
		glPopMatrix();
	} // dc

	// pointcloud
	glPointSize(1.0);
	glColor3f( 0.5, 0.5, 0.5 );
	glBegin(GL_POINTS);
	for (i = 0; i < vcnt; i++){
		glVertex3f(vx[i], vy[i], vz[i]);
	}
	glEnd();

	glPopMatrix();
}

void GraphWindow3DCF::SetRTS2()
{
	//glTranslatef(disptrans2[0], disptrans2[1], disptrans2[2]);

	//float mCameraInv[16];
	//memcpy( mCameraInv, mCamera, sizeof(float)*16 );
	//inv_m44( mCameraInv );
	//glMultMatrixf(mCameraInv);
	glMultMatrixf(mObject);
	//glMultMatrixf(mCamera);
}

void GraphWindow3DCF::draw()
{
	if (!valid())
	{
		initTexture();
		SetLightMat();	// 光源環境設定

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if(projection)
			glFrustum(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);// l,r, b,t, n,f 
		else
			glOrtho(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);

		glClearColor (1.f, 1.f, 1.f, 1.f);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport (0, 0, this->w(), this->h()); // 書かなくてもOK

	if( disp_LIN_SMOOTH ){
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_LINE_SMOOTH);
		//glEnable(GL_POLYGON_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	}

	draw3DCurveFold();

	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);

	glFlush();
	//	glutPostRedisplay();
}

int GraphWindow3DCF::handle(int event)
{
	int i, key;

	switch(event){
	case FL_FOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYDOWN\n");
#endif
		key = Fl::event_key();
		if(key=='r' || key=='R' ){ // reset rotation
			quat2.trackball(prevquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(currquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(dispquat2, 0.0, 0.0, 0.0, 0.0); 
			quat2.build_rotmatrix(mObject, dispquat2);
		} else if( key=='t' || key=='T' ){ // reset translation
			prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
			currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
			disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
		}
		redraw();
		break;
	case FL_KEYUP:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYUP\n");
#endif
		break;
	case FL_PUSH:
		push_button2 = Fl::event_button();
		sts_sft2 = Fl::event_shift();
		sts_alt2 = Fl::event_alt();
		sts_ctrl2 = Fl::event_ctrl();
		push_x2 = Fl::event_x();
		push_y2 = Fl::event_y();

#ifdef DEBUG_MOUSE_EVENT
		printf("FL_PUSH\n");
		printf("sts_sft:%d, sts_alt:%d, sts_ctrl:%d, (%d,%d)\n",
			sts_sft, sts_alt, sts_ctrl, push_x, push_y);
#endif
		if (push_button2 == FL_LEFT_MOUSE && -1 < disp_modTGT && disp_modTGT < ppm->tgcnt) {
			redraw();
			selectidx = selectbuffer[(sbh - 1 - push_y2) * sbstep + push_x2 * 3];
			std::cout << "selectidx = " << selectidx << ", push_x2 = " << push_x2 << ", push_y2 = " << push_y2 << std::endl;
			if (selectidx != 255) {
				disp_modTGT = ((ControlPanel*)cwin)->gwin_cp->disp_modTGT = selectidx;
				prevtrans2[0] = ppm->tgx[selectidx];
				prevtrans2[1] = ppm->tgy[selectidx];
				prevtrans2[2] = ppm->tgz[selectidx];
				redraw();
				((ControlPanel*)cwin)->gwin_cp->redraw();
				return 1; // ドラッグイベントを有効にする
			}
		}
		if(push_button2 == FL_RIGHT_MOUSE && (sts_sft2 || sts_alt2 || sts_ctrl2)){
			return 1; // ドラッグイベントを有効にする
		}else{
			return GraphWindow3D::handle(event);
			//return GraphWindow::handle(event);
		}
		break;
	case FL_DRAG:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_DRAG\n");
		printf("mouse dragged: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if (push_button2 == FL_LEFT_MOUSE
			&& -1 < disp_modTGT && disp_modTGT < ppm->tgcnt
			&& -1 < selectidx && selectidx < 255) {

			// マウス移動量 -> カメラ座標での物体移動量
			currtrans2[0] = Fl::event_x() - push_x2;
			currtrans2[1] = -(Fl::event_y() - push_y2);
			currtrans2[2] = 0;

			// 回転に応じて移動方向を修正
			mult_m44_v4(mCamera, currtrans2);

			// スケールに応じて移動量修正
			currtrans2[0] /= dispscale;
			currtrans2[1] /= dispscale;
			currtrans2[2] /= dispscale;

			// 前の移動量に加算
			disptrans2[0] = prevtrans2[0] + currtrans2[0];
			disptrans2[1] = prevtrans2[1] + currtrans2[1];
			disptrans2[2] = prevtrans2[2] + currtrans2[2];
			ppm->tgx[selectidx] = disptrans2[0];
			ppm->tgy[selectidx] = disptrans2[1];
			ppm->tgz[selectidx] = disptrans2[2];
			redraw();
		}
		if( push_button2 == FL_RIGHT_MOUSE ){
			if(sts_sft2){	// rotation
				// ドラッグした量をcurrquatに保存
				quat2.trackball(currquat2, (2.0*push_x2 - w_w)/w_w, (w_h - 2.0*push_y2)/w_h,
					(2.0*Fl::event_x() - w_w)/w_w, (w_h - 2.0*Fl::event_y())/w_h);
				// prevquat（前回分）とcurrquat（今回分）を足してdispquat（表示用）に格納
				float tmpq0[4], tmpq1[4], camq[4], objq[4], curq[4];
				mat_quat( mCamera, camq );
				mat_quat( mObject, objq );
				camq[3] = -camq[3];
				objq[3] = -objq[3];
				quat2.add_quat(currquat2, camq, tmpq0);
				quat2.add_quat(tmpq0, objq, tmpq1);
				camq[3] = -camq[3];
				objq[3] = -objq[3];
				quat2.add_quat(camq, tmpq1, tmpq0);
				quat2.add_quat(objq, tmpq0, curq);
				quat2.add_quat(prevquat2, curq, dispquat2);

				quat2.build_rotmatrix(mObject, dispquat2);
				mObject[12] = disptrans2[0];
				mObject[13] = disptrans2[1];
				mObject[14] = disptrans2[2];

			} else if(sts_ctrl2){ // translation
				// マウス移動量 -> カメラ座標での物体移動量
				currtrans2[0] = Fl::event_x() - push_x2;
				currtrans2[1] = -(Fl::event_y() - push_y2);
				currtrans2[2] = 0;

				// 回転に応じて移動方向を修正
				mult_m44_v4(mCamera, currtrans2);

				// スケールに応じて移動量修正
				currtrans2[0] /= dispscale;
				currtrans2[1] /= dispscale;
				currtrans2[2] /= dispscale;

				// 前の移動量に加算
				disptrans2[0] = prevtrans2[0] + currtrans2[0];
				disptrans2[1] = prevtrans2[1] + currtrans2[1];
				disptrans2[2] = prevtrans2[2] + currtrans2[2];
				mObject[12] = disptrans2[0];
				mObject[13] = disptrans2[1];
				mObject[14] = disptrans2[2];
			}
			redraw();
		}
		break;
	case FL_RELEASE:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_RELEASE\n");
		printf("mouse released: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if (push_button2 == FL_LEFT_MOUSE && -1 < disp_modTGT && disp_modTGT < ppm->tgcnt) {
			selectidx = -1;
			for (i = 0; i < 3; i++) prevtrans2[i] = 0;
			redraw();
		}
		if( push_button2 == FL_RIGHT_MOUSE && (sts_sft2 || sts_alt2 || sts_ctrl2) )
		{
			if(sts_sft2){
				for(i=0;i<4;i++) prevquat2[i] = dispquat2[i];
			} else if(sts_ctrl2){
				for(i=0;i<3;i++) prevtrans2[i] = disptrans2[i];
			}

			redraw();

			// キー状態を解除
			sts_sft2 = 0;
			sts_ctrl2 = 0;
			push_button2 = 0;
			push_x2 = 0;
			push_y2 = 0;

			return 1;
		}

		// キー状態を解除
		sts_sft2 = 0;
		sts_ctrl2 = 0;
		push_button2 = 0;
		push_x2 = 0;
		push_y2 = 0;

		break;
	case FL_MOVE:
		break;
	default:
		break;
	}

	return GraphWindow3D::handle(event);
}

void GraphWindow3DCF::exportImage(char *filename)
{
	int imgw, imgh, imgstep;
	unsigned char *img=NULL;

	draw(); // redraw immediately

	imgw = this->w();	// may need to be 4*X
	imgh = this->h();
	img = new unsigned char[imgw*imgh*3];
	imgstep = imgw*3;

	glReadPixels(0,0, imgw, imgh, GL_RGB, GL_UNSIGNED_BYTE, img);

	cv::Mat image( cv::Size(imgw, imgh), CV_8UC3 );
	for (int j = 0; j < imgh; j++) {
		for (int i = 0; i < imgw; i++) {
			image.at<cv::Vec3b>(j, i)[0] = img[(imgh-1-j)*imgstep + i*3+2];
			image.at<cv::Vec3b>(j, i)[1] = img[(imgh-1-j)*imgstep + i*3+1];
			image.at<cv::Vec3b>(j, i)[2] = img[(imgh-1-j)*imgstep + i*3];
		}
	}
	imwrite( filename, image);

	delete[] img; img = NULL;
}
