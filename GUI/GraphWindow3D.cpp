#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES 
#include <math.h>

#include <FL/Fl.h>
#include <FL/gl.h>
//#include <GL/glu.h>
#include <GL/glut.h>
#include <FL/Enumerations.H>
#include <FL/fl_draw.H>
#include "GraphWindow3D.h"
#include "../CurvedFoldModel/util.h"
#include "common.h"

GraphWindow3D::GraphWindow3D(int X, int Y, int W, int H, const char *L)
: Fl_Gl_Window(X, Y, W, H, L) {
	w_w = W;
	w_h = H;
	init();

	printf("rotation : shift + left mouse drag\n");
	printf("scaling : alt + left mouse drag\n");
	printf("translation : ctrl + left mouse drag\n");
}

GraphWindow3D::GraphWindow3D(int X, int Y, int W, int H)
: Fl_Gl_Window(X, Y, W, H) {
	w_w = W;
	w_h = H;
	init();

	printf("rotation : shift + left mouse drag\n");
	printf("scaling : alt + left mouse drag\n");
	printf("translation : ctrl + left mouse drag\n");
}

void GraphWindow3D::init()
{
	n_plane = 20;
	f_plane = 4000;
	sts_alt = sts_sft = sts_ctrl = 0;
	push_button = 0;
	push_x = push_y = 0;
	prevquat[0] = prevquat[1] = prevquat[2] = 0; prevquat[3] = 1;
	currquat[0] = currquat[1] = currquat[2] = 0; currquat[3] = 1;
	dispquat[0] = dispquat[1] = dispquat[2] = 0; dispquat[3] = 1;
	mCamera[0]=1;	mCamera[1]=0;	mCamera[2]=0;	mCamera[3]=0;
	mCamera[4]=0;	mCamera[5]=1;	mCamera[6]=0;	mCamera[7]=0;
	mCamera[8]=0;	mCamera[9]=0;	mCamera[10]=1;	mCamera[11]=0;
	mCamera[12]=0;	mCamera[13]=0;	mCamera[14]=0;	mCamera[15]=1;
	prevscale = 1;	prevtrans[0] = prevtrans[1] = prevtrans[2] = 0; prevtrans[3] = 1;
	currscale = 1;	currtrans[0] = currtrans[1] = currtrans[2] = 0; currtrans[3] = 1;
	dispscale = 1;	disptrans[0] = disptrans[1] = disptrans[2] = 0; disptrans[3] = 1;
}

void GraphWindow3D::SetLightMat()
{
	GLfloat light_ambient[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 0.5, 0.5, 0.5, 1.0 };
	//	GLfloat light_position[] = {0.0, 0.0, 0.0, 1.0 }; // positional light, origin = camera pos
	//	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0 }; // directional light, from right-top-front
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	//	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);
	//	glEnable(GL_LIGHTING); // enable before 3D drawings

	GLfloat mat_ambient[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_specular[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat mat_shininess[] = { 10.0 };
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	//	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
}

void GraphWindow3D::SetRTS()
{
	glScalef(dispscale,dispscale,dispscale);
	glMultMatrixf(mCamera);
	//	glTranslatef(disptrans[0]/dispscale, disptrans[1]/dispscale, disptrans[2]/dispscale);
	glTranslatef(disptrans[0], disptrans[1], disptrans[2]);
}

void GraphWindow3D::drawAxis()
{
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);
	glBegin(GL_LINES);
	glColor3f(1.0,0.0,0.0);	glVertex3f(0.0,0.0,0.0);	glVertex3f(100.0,0.0,0.0);
	glColor3f(0.0,1.0,0.0);	glVertex3f(0.0,0.0,0.0);	glVertex3f(0.0,100.0,0.0);
	glColor3f(0.0,0.0,1.0);	glVertex3f(0.0,0.0,0.0);	glVertex3f(0.0,0.0,100.0);
	glEnd();
}

void GraphWindow3D::draw() // to be overloaded
{
	if (!valid())
	{
		SetLightMat();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if(projection)
			glFrustum(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane); 
		else
			glOrtho(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);
		glClearColor (1.0, 1.0, 1.0, 1.0);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport (0, 0, this->w(), this->h());

	//	draw3DShapes();

	glFlush();
}

// up = 0:y-up, 1:z-up, 2:x-up
void GraphWindow3D::resetRotation(int up)
{
	float a[3], phi;

	switch(up){
		case -1: // y-down
			a[0] = 1; a[1] = 0; a[2] = 0;
			phi = (float)M_PI;
			quat.axis_to_quat(a, phi, prevquat);
			quat.axis_to_quat(a, phi, currquat);
			quat.axis_to_quat(a, phi, dispquat);
			break;
		case 0: // y-up
			prevquat[0] = prevquat[1] = prevquat[2] = 0; prevquat[3] = 1;
			currquat[0] = currquat[1] = currquat[2] = 0; currquat[3] = 1;
			dispquat[0] = dispquat[1] = dispquat[2] = 0; dispquat[3] = 1;
			break;
		case 1: // z-up
			a[0] = 1; a[1] = 0; a[2] = 0;
			phi = (float)(M_PI/2.0);
			quat.axis_to_quat(a, phi, prevquat);
			quat.axis_to_quat(a, phi, currquat);
			quat.axis_to_quat(a, phi, dispquat);
			break;
		case 2: // x-up
			a[0] = 0; a[1] = 0; a[2] = 1;
			phi = (float)(-M_PI/2.0);
			quat.axis_to_quat(a, phi, prevquat);
			quat.axis_to_quat(a, phi, currquat);
			quat.axis_to_quat(a, phi, dispquat);
			break;
		case 3: // y-left
			a[0] = 1; a[1] = -1; a[2] = 1;
			phi = (float)-M_PI*2./3.;
			quat.axis_to_quat(a, phi, prevquat);
			quat.axis_to_quat(a, phi, currquat);
			quat.axis_to_quat(a, phi, dispquat);
			break;
	}
	quat.build_rotmatrix(mCamera, dispquat);
}

void GraphWindow3D::resetScale()
{
	prevscale = currscale = dispscale = 1;
}

void GraphWindow3D::resetTranslation()
{
	prevtrans[0] = prevtrans[1] = prevtrans[2] = 0; prevtrans[3] = 1;
	currtrans[0] = currtrans[1] = currtrans[2] = 0; currtrans[3] = 1;
	disptrans[0] = disptrans[1] = disptrans[2] = 0; disptrans[3] = 1;
}

int GraphWindow3D::handle(int event)
{
	int i, key;

	switch(event){
	case FL_FOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYDOWN\n");
#endif
		key = Fl::event_key();
		if(key == 'r'){ // reset rotation
			quat.trackball(prevquat, 0.0, 0.0, 0.0, 0.0);
			quat.trackball(currquat, 0.0, 0.0, 0.0, 0.0);
			quat.trackball(dispquat, 0.0, 0.0, 0.0, 0.0);
			quat.build_rotmatrix(mCamera, dispquat);
		} else if(key == 's'){ // reset scale
			prevscale = currscale = dispscale = 1;
		} else if(key == 't'){ // reset translation
			prevtrans[0] = prevtrans[1] = prevtrans[2] = 0; prevtrans[3] = 1;
			currtrans[0] = currtrans[1] = currtrans[2] = 0; currtrans[3] = 1;
			disptrans[0] = disptrans[1] = disptrans[2] = 0; disptrans[3] = 1;
		}
		redraw();
		break;
	case FL_KEYUP:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYUP\n");
#endif
		break;
	case FL_PUSH:
		push_button = Fl::event_button();
		sts_sft = Fl::event_shift();
		sts_alt = Fl::event_alt();
		sts_ctrl = Fl::event_ctrl();
		push_x = Fl::event_x();
		push_y = Fl::event_y();
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_PUSH\n");
		printf("sts_sft:%d, sts_alt:%d, sts_ctrl:%d, (%d,%d)\n",
			sts_sft, sts_alt, sts_ctrl, push_x, push_y);
#endif
		if(push_button == FL_LEFT_MOUSE && (sts_sft || sts_alt || sts_ctrl)){
			return 1; // ドラッグイベントを有効にする
		}else{
			return Fl_Gl_Window::handle(event);
		}
		break;
	case FL_DRAG:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_DRAG\n");
		printf("mouse dragged: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if( push_button == FL_LEFT_MOUSE ){
			if(sts_sft){
				quat.trackball(currquat, (2.0*push_x - w_w)/w_w, (w_h - 2.0*push_y)/w_h,
					(2.0*Fl::event_x() - w_w)/w_w, (w_h - 2.0*Fl::event_y())/w_h);
				quat.add_quat(currquat, prevquat, dispquat);
				quat.build_rotmatrix(mCamera, dispquat);
			} else if(sts_alt){
				currscale = 1 - (Fl::event_y() - push_y) * 0.1;
				if(fabs(currscale) > 0.001){ // avoid zero
					dispscale = prevscale * currscale;
				}
			} else if(sts_ctrl){
				currtrans[0] = Fl::event_x() - push_x;
				currtrans[1] = -(Fl::event_y() - push_y);
				currtrans[2] = 0;
				mult_m44_v4(mCamera, currtrans);
				currtrans[0] /= dispscale;
				currtrans[1] /= dispscale;
				currtrans[2] /= dispscale;
				disptrans[0] = prevtrans[0] + currtrans[0];
				disptrans[1] = prevtrans[1] + currtrans[1];
				disptrans[2] = prevtrans[2] + currtrans[2];
			}
			redraw();
		}
		break;
	case FL_RELEASE:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_RELEASE\n");
		printf("mouse released: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if( push_button == FL_LEFT_MOUSE && (sts_sft || sts_alt || sts_ctrl) )
		{
			if(sts_sft){
				for(i=0;i<4;i++) prevquat[i] = dispquat[i];
			} else if(sts_alt){
				prevscale = dispscale;
			} else if(sts_ctrl){
				for(i=0;i<3;i++) prevtrans[i] = disptrans[i];
			}

			redraw();

			sts_sft = 0;
			sts_alt = 0;
			sts_ctrl = 0;
			push_button = 0;
			push_x = 0;
			push_y = 0;

			return 1;
		}

		sts_sft = 0;
		sts_alt = 0;
		sts_ctrl = 0;
		push_button = 0;
		push_x = 0;
		push_y = 0;

		break;
	case FL_MOVE:
		//printf("FL_MOVE\n");
		break;
	default:
		break;
	}
	return Fl_Gl_Window::handle(event);
}
