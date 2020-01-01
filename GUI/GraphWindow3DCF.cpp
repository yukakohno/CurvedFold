#include <stdio.h>
#include <stdlib.h>
#include <iostream>

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
}

GraphWindow3DCF::~GraphWindow3DCF()
{
	init();
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
#else
	IplImage *tex = cvLoadImage( TEXFNAME, CV_LOAD_IMAGE_COLOR );
	if( tex ){
		for(int j=0 ; j<TEXHEIGHT; j++){
			for(int i=0 ; i<TEXWIDTH; i++){
				texture[j][i][2] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep]);
				texture[j][i][1] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep +1]);
				texture[j][i][0] = (GLubyte)(tex->imageData[i*tex->nChannels +j*tex->widthStep +2]);
				texture[j][i][3] = (GLubyte)255;
			}
		}
	}
	cvReleaseImage(&tex);
#endif

	/* テクスチャ画像はバイト単位に詰め込まれている */
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	/* テクスチャを拡大・縮小する方法の指定 */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	/* テクスチャの繰り返し方法の指定 */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

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

	disp_axis = 0;
	disp_X = 1;
	disp_TNB = 0;
	disp_R = 1;
	disp_maxrlen = 0;
	disp_PLY = 0;
	disp_PTN = 1;
	disp_CP = 0;	ppos = 0;	pprm = 0;

	flg_addcurve = 0;

	ppm = NULL;

	projection=0;

	push_button2 = 0;
	push_x2 = push_y2 = 0;
	initObject();

	vcnt = 0;
	vx = vy = vz = NULL;

	resetRotation(-1);
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

	// 原点を見るように, カメラ位置を(0,0,n_plane*2)に
	// ＝これから描画するオブジェクトを(0.,0.,-n_plane*2)平行移動
	//	glTranslated(0.,0.,-n_plane*2);
	glTranslated(0.,0.,-n_plane-3000);

	// 視点設定
	SetRTS();

	if(disp_axis)
		drawAxis();

	glPushMatrix();
	SetRTS2();

	//
	// control points
	//
	if( disp_CP && ppos>-1 && pprm>-1 ){
		switch(pprm){	// P_CV2D, P_CV3D, P_TRSN, P_FLDA, P_RULL, P_RULR
			case P_CV2D:	glColor3f( 1.0, 0.0, 1.0 );	break;
			case P_CV3D:	glColor3f( 0.0, 1.0, 0.0 );	break;
			case P_TRSN:	glColor3f( 0.0, 0.0, 1.0 );	break;
			case P_FLDA:	glColor3f( 1.0, 0.0, 0.0 );	break;
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
		glLineWidth(4.0);
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
	// curve
	//
	if( disp_X && ppm ){
		glLineWidth(2.0);
		glColor3f( 0.2, 0.2, 0.2 );
		for( j=0; j<ppm->crcnt; j++ ){
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
		for( j=0; j<ppm->crcnt; j++ ){
			crease *c = &(ppm->crs[j]);
			if( c->rl<=0 ){
				//glColor3f( 0.6, 0.4, 0.9 ); // Lavender: B57EDC
				glColor3f( 0.32, 0.12, 0.46 ); // Dark Violet: 522076
				//glColor3f( 1.0, 1.0, 0.0 );
				for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
					glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
					glVertex3f(c->Xx[i] + c->rlx[i]*c->rllen[i],
						c->Xy[i] + c->rly[i]*c->rllen[i],
						c->Xz[i] + c->rlz[i]*c->rllen[i] );
				}
			}
			if(c->rl>=0) {
				glColor3f( 0.9, 0.0, 0.4 ); // Raspberry: E30B5D
				//glColor3f( 1.0, 1.0, 0.0 );
				for( i=c->Xsidx; i<c->Xeidx+1; i++ ){
					glVertex3f(c->Xx[i], c->Xy[i], c->Xz[i]);
					glVertex3f(c->Xx[i] + c->rrx[i]*c->rrlen[i],
						c->Xy[i] + c->rry[i]*c->rrlen[i],
						c->Xz[i] + c->rrz[i]*c->rrlen[i] );
				}
			}
		}
		glEnd();
	}

	//
	// Polygons
	//
	if( disp_PLY && ppm ){
		glLineWidth(1.0);
		glColor3f( 0.0, 0.0, 0.0 );
		for( i=0; i<ppm->plcnt; i++ ){
			glBegin(GL_LINE_STRIP);
			for( j=0; j<ppm->plvcnt[i]; j++ ){
				glVertex3f(ppm->plx[i*4+j], ppm->ply[i*4+j], ppm->plz[i*4+j] );
			}
			glVertex3f(ppm->plx[i*4], ppm->ply[i*4], ppm->plz[i*4] );
			glEnd();
		}
	}

	if( disp_PTN && ppm ){
		static const GLfloat color[] = { 1.0, 1.0, 1.0, 1.0 };  /* 材質 (色) */
		static const GLfloat lightcol[] = { 1.0, 1.0, 1.0, 1.0 }; /* 直接光強度 */
		static const GLfloat lightamb[] = { 0.5, 0.5, 0.5, 1.0 }; /* 環境光強度 */

		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
		glEnable(GL_ALPHA_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_NORMALIZE);
		glEnable(GL_LIGHT0);
		/* 光源の位置を設定 */
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightcol);
		glLightfv(GL_LIGHT0, GL_SPECULAR, lightcol);
		glLightfv(GL_LIGHT0, GL_AMBIENT, lightamb);

		/* テクスチャマッピング開始 */
		glEnable(GL_TEXTURE_2D);

#if 1	// CP上のポリゴンを座標変換して表示
		glColor3f( 1.0, 0.0, 1.0 );
		for( i=0; i<ppm->plcnt; i++ ){
			glPushMatrix();
			glMultMatrixd( &(ppm->plmat[i*16]) );
			glBegin(GL_POLYGON);
			glNormal3f( ppm->plmat[i*16+8], ppm->plmat[i*16+9], ppm->plmat[i*16+10] );
			for( j=0; j<ppm->plvcnt[i]; j++ ){
				glTexCoord2d( (ppm->plx_cp[i*4+j]-ppm->psx)/(double)PPWIDTH, (ppm->ply_cp[i*4+j]-ppm->psy)/(double)PPHEIGHT );
				glVertex3f(ppm->plx_cp[i*4+j], ppm->ply_cp[i*4+j], 0.0 );
			}
			glTexCoord2d( (ppm->plx_cp[i*4]-ppm->psx)/(double)PPWIDTH, (ppm->ply_cp[i*4]-ppm->psy)/(double)PPHEIGHT );
			glVertex3f(ppm->plx_cp[i*4], ppm->ply_cp[i*4], 0.0 );
			glEnd();
			glPopMatrix();
		}
#else	// 3Dポリゴンを表示
		glColor3f( 1.0, 0.0, 1.0 );
		for( i=0; i<ppm->plcnt; i++ ){
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
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_NORMALIZE);
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
	}
	glPopMatrix();

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
	glTranslatef(disptrans2[0], disptrans2[1], disptrans2[2]);
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
		glEnable(GL_ALPHA_TEST);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if(projection)
			glFrustum(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);// l,r, b,t, n,f 
		else
			glOrtho(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);

		glClearColor (1.f, 1.f, 1.f, 1.f);

		static const GLfloat lightpos[] = { 1.0, 0.0, -1.0, 0.0 }; /* 位置　　　 */
		glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport (0, 0, this->w(), this->h()); // 書かなくてもOK

	draw3DCurveFold();

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
		if(key == 'R'){ // reset rotation
			quat2.trackball(prevquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(currquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(dispquat2, 0.0, 0.0, 0.0, 0.0); 
			quat2.build_rotmatrix(mObject, dispquat2);
		} else if(key == 'T'){ // reset translation
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
			}
			redraw();
		}
		break;
	case FL_RELEASE:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_RELEASE\n");
		printf("mouse released: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
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
