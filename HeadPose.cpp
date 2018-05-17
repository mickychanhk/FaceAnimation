// HeadPose.cpp : Defines the entry point for the console application.
//
//#include <glut.h>
#include <glew.h>
#include "stdafx.h"

#include "opencv\cv.h"
#include "opencv\highgui.h"
//using namespace cv;

#include <vector>
#include <iostream>
#include <fstream>
#include <dlib\optimization.h>
#include <dlib/opencv.h>
#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/image_processing/render_face_detections.h>
#include <dlib/image_processing.h>
#include <dlib/gui_widgets.h>
#include <dlib/image_io.h>
using namespace std;
//#if defined(__APPLE__)
//#  include <OpenGL/gl.h>
//#  include <OpenGL/glu.h>
//#elif defined(__linux__) || defined(__MINGW32__) || defined(WIN32)
//#  include <GL/gl.h>
//#  include <GL/glu.h>
//#else
//#include <gl\GL.h>
//#include <gl\GLU.h>

//#endif

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include "glm.h"
#include "OGL_OCV_common.h"
#include "BlendShape.h"
#include "Tensor.hpp"
#include "nnls.h"
using namespace Eigen;

void loadNext();
void loadWithPoints(cv::Mat& ip, cv::Mat& img, vector<cv::Point2f> imagepoints);

const GLfloat light_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 0.0f, 0.0f, 1.0f, 0.0f };

const GLfloat mat_ambient[] = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };

double rot[9] = { 0 };
GLuint textureID;
cv::Mat backPxls;
vector<double> rv(3), tv(3);
cv::Mat rvec(rv), tvec(tv);
cv::Mat camMatrix;
OpenCVGLTexture imgTex, imgWithDrawing;
cv::VideoCapture webcam(0);
vector<GLMmodel*> empty_face;
GLMmodel* average_face;
vector<GLMmodel*> blendshapes;
vector<vector<GLMmodel*>> rank3Tensor;
vector<vector<double>> wIdX, wIdY, wIdZ;
vector<vector<double>>  wExpX, wExpY, wExpZ;
dlib::shape_predictor sp;
dlib::frontal_face_detector detector;
int chaningRadio = 1;
int faceDesign = 0;
cv::Mat op;
cv::Mat op_min;
vector<cv::Point2f> imgPoint;
vector<cv::Point2f> prevImgPoint, prevLandmarkPts;
int numOfuser = 2;
int numOfBlendshape = 47;
int numOfVertex = 11510;
Tensor3<double> coreTensor;
vector<Tensor2<double>> UIdAndUExp;
Tensor3<double> tempWExpTensor;
Tensor3<double> tempWIdTensor;
Tensor3<double> initTensor;
double calDistanceDiff(vector<cv::Point2f> curPoints, vector<cv::Point2f> lastPoints) {
	double variance = 0.0;
	double sum = 0.0;
	std::vector<double> diffs;
	if (curPoints.size() == lastPoints.size()) {
		for (int i = 0; i < curPoints.size(); i++) {
			double diff = std::sqrt(std::pow(curPoints[i].x - lastPoints[i].x, 2.0) + std::pow(curPoints[i].y - lastPoints[i].y, 2.0));
			sum += diff;
			diffs.push_back(diff);
		}
		double mean = sum / diffs.size();
		for (int i = 0; i < curPoints.size(); i++) {
			variance += std::pow(diffs[i] - mean, 2);
		}
		return variance / diffs.size();
	}
	return variance;
}
vector<double> findWid(vector<int> vertexPoint, vector<double> wexp, vector<double> wid) {
	cv::Mat_<double> transformationMatrix = (cv::Mat_<double>(3, 4) << 0.1, 0.1, 0.1, 0.1
		, 0.1, 0.1, 0.1, 0.1,
		0.1, 0.1, 0.1, 0.1);
	vector<double> widtest(numOfuser), wexptest(numOfBlendshape);
	widtest.push_back(1.0);
	widtest.push_back(0);
	wexptest.push_back(1);
	for (int i = 1; i < numOfBlendshape; i++) {
		wexptest.push_back(0);
	}
	initTensor = initTensor.modeProduct(widtest, 0).modeProduct(wexptest, 1);

	for (int k = 0; k < imgPoint.size(); k++) {
		/*double pX = initTensor(0, 0, 3 * vertexPoint[k]);
		double pY = initTensor(0, 0, 3 * vertexPoint[k] + 1);
		double pZ = initTensor(0, 0, 3 * vertexPoint[k] + 2);*/
		double pX = average_face->vertices[3 * vertexPoint[k]];
		double pY = average_face->vertices[3 * vertexPoint[k] + 1];
		double pZ = average_face->vertices[3 * vertexPoint[k] + 2];
	}
	//double rotss[9] = { 0 };
	//cv::Mat rotMss(3, 3, CV_64FC1, rotss);
	//// convert the vector to matrix
	//Rodrigues(rvec, rotMss);
	//double* _r = rotMss.ptr<double>();
	//double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
	//	_r[3],_r[4],_r[5],tv[1],
	//	_r[6],_r[7],_r[8],tv[2] };
	//cv::Matx34d P(_pm);
	//cv::Mat KP = camMatrix * cv::Mat(P);
	////double y0 = orignialFormula(wid, wexp, vertexPoint, KP);
	//int idx = 0;

	//if (imgPoint.size() != 0) {
	//	/*arma::vec b(imgPoint.size());
	//	arma::mat C(imgPoint.size(), numOfBlendshape);*/
	//	int numLM = 68;

	//	//for (int i = 0; i < numOfuser; i++) {
	//	for (int itr = 0; itr < 3; itr++) {
	//		MatrixXd A(numLM * 2, numOfuser);
	//		VectorXd b(numLM * 2);
	//		vector<double> widList;
	//		for (int j = 0; j < numOfuser; j++) {
	//			for (int k = 0; k < imgPoint.size(); k++) {
	//				/*double pX = rank3Tensor[j][0]->vertices[3 * vertexPoint[k]];
	//				double pY = rank3Tensor[j][0]->vertices[3 * vertexPoint[k] + 1];
	//				double pZ = rank3Tensor[j][0]->vertices[3 * vertexPoint[k] + 2];*/
	//				double pX = initTensor(j, 0, 3 * vertexPoint[k])* wexp[0] * wid[j];
	//				double pY = initTensor(j, 0, 3 * vertexPoint[k] + 1) * wexp[0] * wid[j];
	//				double pZ = initTensor(j, 0, 3 * vertexPoint[k] + 2) * wexp[0] * wid[j];
	//				cv::Mat_<double> X = (cv::Mat_<double>(4, 1)
	//					<< pX,
	//					pY,
	//					pZ, 1.0);//convert point to vector
	//				cv::Mat_<double> opt_m = KP * X; // KP = camera Matrix * Pose Matrix
	//												 //cout << "point x" << opt_m(0) << "point y" << opt_m(1) << "point z" << opt_m(2) << endl;
	//				double finalX = opt_m(0) / opt_m(2);
	//				double finalY = opt_m(1) / opt_m(2);
	//				A(k, j) = finalX;
	//				A(k + vertexPoint.size(), j) = finalY;
	//				b(k) = imgPoint[k].x;
	//				b(k + vertexPoint.size()) = imgPoint[k].y;

	//			}
	//		}
	//		VectorXd x(numOfuser);
	//		//Solving Ax = b issues with non-negative least square
	//		NNLS<MatrixXd>::solve(A, b, x);
	//		//x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
	//		double totalSum = 0;
	//		for (int i = 0; i < x.size(); i++) {
	//			widList.push_back(x(i));
	//			//widList.push_back(solution(i));
	//			//totalSum += fabs(solution(i));
	//		}
	//		wid = widList;
	//	}
	//	//for (int i = 0; i < x.size(); i++) {
	//	//	widList.push_back(x(i));
	//	//	/*widList.push_back(fabs(solution(i)));
	//	//	totalSum += fabs(solution(i));*/
	//	//}
	//	/*for (int i = 0; i < solution.size(); i++) {
	//		widList[i] = widList[i] / totalSum;
	//	}*/
	//}
	return wid;
}

vector<double> findWExp(vector<int> vertexPoint, vector<double> wexp, vector<double> wid) {
	double rotss[9] = { 0 };
	cv::Mat rotMss(3, 3, CV_64FC1, rotss);
	// convert the vector to matrix
	Rodrigues(rvec, rotMss);
	double* _r = rotMss.ptr<double>();
	double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
		_r[3],_r[4],_r[5],tv[1],
		_r[6],_r[7],_r[8],tv[2] };
	cv::Matx34d P(_pm);
	cv::Mat KP = camMatrix * cv::Mat(P);
	int idx = 0;

	if (imgPoint.size() != 0) {
		int numBS = 47;
		int numLM = 68;
		for (int itr = 0; itr < 3; itr++) {
			vector<double> wexpList;
			MatrixXd A(numLM * 2, numBS);
			VectorXd b(numLM * 2);
			//for (int i = 0; i < numOfuser; i++) {
			for (int j = 0; j < numBS; j++) {
				for (int k = 0; k < imgPoint.size(); k++) {
					/*	double pX = rank3Tensor[0][j]->vertices[3 * vertexPoint[k]];
						double pY = rank3Tensor[0][j]->vertices[3 * vertexPoint[k] + 1];
						double pZ = rank3Tensor[0][j]->vertices[3 * vertexPoint[k] + 2];*/
					double pX = initTensor(0, j, 3 * vertexPoint[k]) * wexp[j] * wid[0];// -tempWExpTensor(0, 0, vertexPoint[k]);
					double pY = initTensor(0, j, 3 * vertexPoint[k] + 1) * wexp[j] * wid[0];// -tempWExpTensor(0, 0, vertexPoint[k] + 1);
					double pZ = initTensor(0, j, 3 * vertexPoint[k] + 2) * wexp[j] * wid[0];// -tempWExpTensor(0, 0, vertexPoint[k] + 2);
					cv::Mat_<double> X = (cv::Mat_<double>(4, 1)
						<< pX,
						pY,
						pZ, 1.0);//convert point to vector
					cv::Mat_<double> opt_m = KP * X; // KP = camera Matrix * Pose Matrix
					double finalX = opt_m(0) / opt_m(2);
					double finalY = opt_m(1) / opt_m(2);
					A(k, j) = finalX;
					A(k + vertexPoint.size(), j) = finalY;
					//tempA[k + vertexPoint.size()][j] = finalY;
					if (j == 0) {
						b(k) = imgPoint[k].x;
						b(k + vertexPoint.size()) = imgPoint[k].y;
						//	tempB[k + vertexPoint.size()] = imgPoint[k].y;;
					}
				}
			}

			VectorXd x(numBS);
			//Solving Ax = b issues with non-negative least square
			NNLS<MatrixXd>::solve(A, b, x);
			//x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
			double totalSum = 0;
			for (int i = 0; i < x.size(); i++) {
				//wexpList.push_back(solution(i));
				//dlib::clamp_function(x(i), 0, 1);
				wexpList.push_back(x(i));
				//totalSum += fabs(solution(i));
			}
			wexp = wexpList;
		}
		/*for (int i = 0; i < solution.size(); i++) {
			wexpList[i] = wexpList[i] / totalSum;
		}*/
	}
	return wexp;
}

double roundFloat(double var)
{
	// we use array of chars to store number
	// as a string.
	char str[80];

	// Print in string the value of var 
	// with two decimal point
	sprintf(str, "%.6lf", var);
	double lf2 = atof(str);
	// scan string value in var 
	//sscanf(str, "%f", &var);

	return lf2;
}

void saveOpenGLBuffer() {
	static unsigned int opengl_buffer_num = 0;

	int vPort[4]; glGetIntegerv(GL_VIEWPORT, vPort);
	cv::Mat_<cv::Vec3b> opengl_image(vPort[3], vPort[2]);
	{
		cv::Mat_<cv::Vec4b> opengl_image_4b(vPort[3], vPort[2]);
		glReadPixels(0, 0, vPort[2], vPort[3], GL_BGRA, GL_UNSIGNED_BYTE, opengl_image_4b.data);
		cv::flip(opengl_image_4b, opengl_image_4b, 0);
		mixChannels(&opengl_image_4b, 1, &opengl_image, 1, &(cv::Vec6i(0, 0, 1, 1, 2, 2)[0]), 3);
	}
	stringstream ss; ss << "opengl_buffer_" << opengl_buffer_num++ << ".jpg";
	imwrite(ss.str(), opengl_image);
}
void resize(int width, int height)
{
	const float ar = (float)width / (float)height;

	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
	gluPerspective(47, 1.0, 0.01, 1000.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
int __w = 640, __h = 480;
void key(unsigned char key, int x, int y)
{

	switch (key)
	{
	case 27:
	case 'Q':
	case 'q':
		exit(0);
		break;
	case 'w':
	case 'W':
		__w++;
		__w = __w % 251;
		break;
	case 'h':
	case 'H':
		__h++;
		__h = __h % 250;
		break;
	case ' ':
		saveOpenGLBuffer();
		//loadNext();
		break;
	case 'z':
		break;
	case 'x':
		break;
	case 'c':
		//printf("current weight point : %f %f %f ", weightX, weightY, weightZ);
		break;
		//case 'x':
		//	weightX--;
		//	break;
	default:
		break;
	}

	glutPostRedisplay();
}
void idle(void)
{
	//loadNext();
	glutPostRedisplay();
}
void myGLinit() {
	//    glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);


	glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

	glEnable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}
void drawAxes() {

	//Z = red
	glPushMatrix();
	glRotated(180, 0, 1, 0);
	glColor4d(1, 0, 0, 0.5);
	//	glutSolidCylinder(0.05,1,15,20);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
	glEnd();
	glTranslated(0, 0, 1);
	glScaled(.1, .1, .1);
	glutSolidTetrahedron();
	glPopMatrix();

	//Y = green
	glPushMatrix();
	glRotated(-90, 1, 0, 0);
	glColor4d(0, 1, 0, 0.5);
	//	glutSolidCylinder(0.05,1,15,20);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
	glEnd();
	glTranslated(0, 0, 1);
	glScaled(.1, .1, .1);
	glutSolidTetrahedron();
	glPopMatrix();

	//X = blue
	glPushMatrix();
	glRotated(-90, 0, 1, 0);
	glColor4d(0, 0, 1, 0.5);
	//	glutSolidCylinder(0.05,1,15,20);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
	glEnd();
	glTranslated(0, 0, 1);
	glScaled(.1, .1, .1);
	glutSolidTetrahedron();
	glPopMatrix();
}

void display(void)
{
	// draw the image in the back
	int vPort[4]; glGetIntegerv(GL_VIEWPORT, vPort);
	glEnable2D();
	//drawOpenCVImageInGL(imgTex);
	//glTranslated(vPort[2] / 2.0, 0, 1);
	//drawOpenCVImageInGL(imgWithDrawing);
	glDisable2D();

	glClear(GL_DEPTH_BUFFER_BIT); // we want to draw stuff over the image

	// draw only on left part
	glViewport(0, 0, vPort[2], vPort[3]);

	glPushMatrix();

	gluLookAt(0, 0, 0, 0, 0, 1, 0, -1, 0);

#pragma region vertex index  
	vector<cv::Point3f > modelPoints;
	vector<int> vertexPoint;
	vertexPoint.push_back(2179);	 // l eye (v 8083) shape.part(44).x(), shape.part(44).y()));	 // l eye (v 314)
	vertexPoint.push_back(6886);	//			       shape.part(45).x(), shape.part(45).y()));	// l eye (v 314)
	vertexPoint.push_back(2044);	//			       shape.part(46).x(), shape.part(46).y()));	// l eye (v 314)
	vertexPoint.push_back(6825);	//			       shape.part(47).x(), shape.part(47).y()));	// l eye (v 314)
	vertexPoint.push_back(7228);	//			       shape.part(48).x(), shape.part(48).y()));	// l eye (v 314)
	vertexPoint.push_back(7246);	//			       shape.part(49).x(), shape.part(49).y()));	// l eye (v 314)

	vertexPoint.push_back(576);	 // r eye (v 2343)  shape.part(36).x(), shape.part(36).y()));	 // r eye (v 0)
	vertexPoint.push_back(618);	//					shape.part(37).x(), shape.part(37).y()));	 // r eye (v 0)
	vertexPoint.push_back(613);	//					shape.part(38).x(), shape.part(38).y()));	 // r eye (v 0)
	vertexPoint.push_back(754);	//					shape.part(39).x(), shape.part(39).y()));	 // r eye (v 0)
	vertexPoint.push_back(4352);	//					shape.part(40).x(), shape.part(40).y()));	 // r eye (v 0)
	vertexPoint.push_back(4333);	//					shape.part(41).x(), shape.part(41).y()));	 // r eye (v 0)

	vertexPoint.push_back(4241);	 //nose (v 1264) shape.part(27).x(), shape.part(27).y()));	 //nose (v 1879)
	vertexPoint.push_back(9198);	//				 shape.part(28).x(), shape.part(28).y()));	 //nose (v 1879)
	vertexPoint.push_back(6293);	//				 shape.part(29).x(), shape.part(29).y()));	 //nose (v 1879)
	vertexPoint.push_back(8177);	//				 shape.part(30).x(), shape.part(30).y()));	 //nose (v 1879)
	vertexPoint.push_back(10491);	//				 shape.part(31).x(), shape.part(31).y()));	 //nose (v 1879)
	vertexPoint.push_back(355);	//				 shape.part(32).x(), shape.part(32).y()));	 //nose (v 1879)
	vertexPoint.push_back(3563);	//				 shape.part(33).x(), shape.part(33).y()));	 //nose (v 1879)
	vertexPoint.push_back(1780);	//				 shape.part(34).x(), shape.part(34).y()));	 //nose (v 1879)
	vertexPoint.push_back(6462);	//				 shape.part(35).x(), shape.part(35).y()));	 //nose (v 1879)

	vertexPoint.push_back(1601);	 // l mouth shape.part(53).x(), shape.part(51).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(2812);	 // l mouth shape.part(54).x(), shape.part(52).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6220);	 // l mouth shape.part(55).x(), shape.part(53).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(1594);	 // l mouth shape.part(56).x(), shape.part(54).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6155);	 // l mouth shape.part(57).x(), shape.part(55).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6076);	 // l mouth shape.part(58).x(), shape.part(56).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6074);	 // l mouth shape.part(59).x(), shape.part(57).y()));	 // l mouth (v 1502)

	vertexPoint.push_back(6201);	//	 shape.part(62).x(), shape.part(62).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6092);	//	 shape.part(63).x(), shape.part(63).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(1649);	//	 shape.part(64).x(), shape.part(64).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6056);	//	 shape.part(65).x(), shape.part(65).y()));	 // l mouth (v 1502)
	vertexPoint.push_back(6064);	//	 shape.part(66).x(), shape.part(66).y()));	 // l mouth (v 1502)

	vertexPoint.push_back(3211);		 // r mouth shape.part(50).x(), shape.part(50).y()));	 // r mouth (v 695)
	vertexPoint.push_back(3239);		 // r mouth shape.part(49).x(), shape.part(49).y()));	 // r mouth (v 695)
	vertexPoint.push_back(3184);		 // r mouth shape.part(48).x(), shape.part(48).y()));	 // r mouth (v 695)
	vertexPoint.push_back(3194);		 // r mouth shape.part(59).x(), shape.part(59).y()));	 // r mouth (v 695)
	vertexPoint.push_back(10325);		 // r mouth shape.part(58).x(), shape.part(58).y()));	 // r mouth (v 695)

	vertexPoint.push_back(3210);		 // r mouth shape.part(61).x(), shape.part(61).y()));	 // r mouth (v 695)
	vertexPoint.push_back(254);		 // r mouth shape.part(60).x(), shape.part(60).y()));	 // r mouth (v 695)
	vertexPoint.push_back(10326);		 // r mouth shape.part(67).x(), shape.part(67).y()));	 // r mouth (v 695)

	vertexPoint.push_back(6880);	 // l contour shape.part(16).x(), shape.part(16).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(7078);	 // l contour shape.part(15).x(), shape.part(15).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(2436);	 // l contour shape.part(14).x(), shape.part(14).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(2557);	 // l contour shape.part(13).x(), shape.part(13).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(6748);	 // l contour shape.part(12).x(), shape.part(12).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(9088);	 // l contour shape.part(11).x(), shape.part(11).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(6504);	 // l contour shape.part(10).x(), shape.part(10).y()));	 	// l coutour (v 2011)
	vertexPoint.push_back(6541);	 // l contour shape.part(9).x(), shape.part(9).y()));	 	    // l coutour (v 2011)
	vertexPoint.push_back(9053);	 // l contour shape.part(8).x(), shape.part(8).y()));	 	    // l coutour (v 2011)

	vertexPoint.push_back(646);	 // r contour (v 3484) shape.part(0).x(), shape.part(0).y()));	 	// r coutour (v 1138
	vertexPoint.push_back(4194);	 // r contour (v 3484) shape.part(1).x(), shape.part(1).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(3917);	 // r contour (v 3484) shape.part(2).x(), shape.part(2).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(3638);	 // r contour (v 3484) shape.part(3).x(), shape.part(3).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(10592);	 // r contour (v 3484) shape.part(4).x(), shape.part(4).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(11486);	 // r contour (v 3484) shape.part(5).x(), shape.part(5).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(10537);	 // r contour (v 3484) shape.part(6).x(), shape.part(6).y()));	 	// r coutour (v 1138)
	vertexPoint.push_back(477);	 // r contour (v 3484) shape.part(7).x(), shape.part(7).y()));	 	// r coutour (v 1138)

	vertexPoint.push_back(4272);	 // r eyebrow(v 1138) shape.part(17).x(), shape.part(17).y()));	 	// r eyebrow(v 1138)
	vertexPoint.push_back(4296);	 // r eyebrow(v 1138) shape.part(18).x(), shape.part(18).y()));	 	// r eyebrow(v 1138)
	vertexPoint.push_back(4451);	 // r eyebrow(v 1138) shape.part(19).x(), shape.part(19).y()));	 	// r eyebrow(v 1138)
	vertexPoint.push_back(4443);	 // r eyebrow(v 1138) shape.part(20).x(), shape.part(20).y()));	 	// r eyebrow(v 1138)
	vertexPoint.push_back(4295);	 // r eyebrow(v 1138) shape.part(21).x(), shape.part(21).y()));	 	// r eyebrow(v 1138)

	vertexPoint.push_back(1045);	 // l eyebrow(v 1138) shape.part(22).x(), shape.part(22).y()));
	vertexPoint.push_back(2217);	 // l eyebrow(v 1138) shape.part(23).x(), shape.part(23).y()));
	vertexPoint.push_back(7347);	 // l eyebrow(v 1138) shape.part(24).x(), shape.part(24).y()));
	vertexPoint.push_back(7338);	 // l eyebrow(v 1138) shape.part(25).x(), shape.part(25).y()));
	vertexPoint.push_back(2151);	 // l eyebrow(v 1138) shape.part(26).x(), shape.part(26).y()));
#pragma endregion vertex index




//	GLMmodel* original_obj = new GLMmodel();
	//Store the neutral face vertex for reset the neutral face after drawing
	vector<GLfloat> originalVertex;
	vector<double> wexpList(numOfBlendshape), widList(numOfuser);
	wexpList[0] = 0;
	wexpList[1] = 1;
	for (int i = 2; i < numOfBlendshape; i++) {
		wexpList[i] = 0;
	}
	widList[0] = 0.99;
	for (int i = 1; i < numOfuser; i++) {
		widList[i] = 0;
	}
	vector<double> wid, wexp;
	//wid = findWid(vertexPoint, wexpList, widList);
	//wexp = findWExp(vertexPoint, wexpList, widList);

	////double sum = 0;
	////for (int i = 0; i < wid.size(); i++) {
	////	sum += wid[i];
	////}
	////cout << sum << endl;
	////cout << sum << endl;
	////cout << opt_wid << endl;
	for (int i = 1; i <= numOfVertex; i++) {
		originalVertex.push_back(average_face->vertices[3 * i]);
		originalVertex.push_back(average_face->vertices[3 * i + 1]);
		originalVertex.push_back(average_face->vertices[3 * i + 2]);
		average_face->vertices[3 * i] = 0;
		average_face->vertices[3 * i + 1] = 0;
		average_face->vertices[3 * i + 2] = 0;
	}

	if (wexpList.size() > 0) {
		for (int i = 0; i < numOfuser; i++) {
			for (int j = 0; j < numOfBlendshape; j++) {
				for (int k = 0; k < numOfVertex; k++) {
					//double faceX = rank3Tensor[i][j]->vertices[3 * (k)];// -rank3Tensor[i][0]->vertices[3 * k];
					//double faceY = rank3Tensor[i][j]->vertices[3 * (k)+1];// -rank3Tensor[i][0]->vertices[3 * k + 1];
					//double faceZ = rank3Tensor[i][j]->vertices[3 * (k)+2];// -rank3Tensor[i][0]->vertices[3 * k + 2];
					double faceX = initTensor(i, j, 3 * k) * widList[i] * wexpList[j];
					double faceY = initTensor(i, j, 3 * k + 1) * widList[i] * wexpList[j];
					double faceZ = initTensor(i, j, 3 * k + 2) * widList[i] * wexpList[j];
					//The first 3 item in the vertices is useless, so we start from 3 to set x, 4 to set y, 5 to set z in the vertices array
					int idx = k + 1;
					average_face->vertices[3 * idx] += faceX;
					average_face->vertices[3 * idx + 1] += faceY;
					average_face->vertices[3 * idx + 2] += faceZ;
				}
			}
		}
	}
	vector<int> calculatePose;
	calculatePose.push_back(1967);
	calculatePose.push_back(542);
	calculatePose.push_back(6328);
	calculatePose.push_back(1594);	 // l mouth (v 6601) or line 6605
	calculatePose.push_back(187);		 // r mouth (v 728) or line 732 
	calculatePose.push_back(6878);	 // l contour (v 8699)
	calculatePose.push_back(685);	 // r contour (v 3484)
	for (int i = 0; i < calculatePose.size(); i++) {
		double x = roundFloat(average_face->vertices[3 * calculatePose[i]]);
		double y = roundFloat(average_face->vertices[3 * calculatePose[i] + 1]);
		double z = roundFloat(average_face->vertices[3 * calculatePose[i] + 2]);
		modelPoints.push_back(cv::Point3d(
			x,
			y,
			z
		));
	}
	op_min = cv::Mat(modelPoints);
	loadNext();
	// draw the 3D head model
	glColor4f(1, 1, 1, 0.75);

	cv::Vec3d tvv(tv[0], tv[1], tv[2]);
	glTranslated(tvv[0], tvv[1], tvv[2]);
	//glTranslated(0, 0, 5);
	// control the head rotation
	double _d[16] = { rot[0],rot[1],rot[2],0,
		rot[3],rot[4],rot[5],0,
		rot[6],rot[7],rot[8],0,
		0,	   0,	  0		,1 };
	glMultMatrixd(_d);
	glmDraw(average_face, GLM_SMOOTH);

	//Reset the neutral face to the orignial
	for (int i = 0; i < originalVertex.size(); i++) {
		average_face->vertices[3 + i] = originalVertex[i];
	}

	//----------Axes
	//glScaled(.5, .5,.5);
	//drawAxes();
	//----------End axes

	glPopMatrix();
	// restore to looking at complete viewport
	glutSwapBuffers();
}

void init_opengl(int argc, char** argv) {
	glutInitWindowSize(640, 640);
	glutInitWindowPosition(40, 40);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH); // | GLUT_MULTISAMPLE
	glutCreateWindow("head pose");

	myGLinit();

	glutReshapeFunc(resize);
	glutDisplayFunc(display);
	glutKeyboardFunc(key);
	//glutSpecialFunc(special);
	glutIdleFunc(idle);
}

int start_opengl() {

	glutMainLoop();

	return 1;
}


void loadNext() {
	static int counter = 1;

	//printf("load %d\n", counter);

	const char* workingDir = "./";
	//Micky 

	cv::Mat temp;
	//temp = cv::imread("Stephen-640x480.jpg");
	webcam >> temp;
	cv::Mat dst = temp;
	//cv::resize(temp, dst, cv::Size(352, 240));
	dlib::cv_image <dlib::bgr_pixel > cimg(dst);
	std::vector<dlib::rectangle> dets = detector(cimg);
	if (dets.size() != 0) {
		dlib::full_object_detection shape = sp(cimg, dets[0]);
		//char buf[256] = { 0 };
		//sprintf(buf, "%sAngelina_Jolie/Angelina_Jolie_%04d.txt", workingDir, counter);

		vector<cv::Point2f> imagePoints;
		/*ifstream inputfile(buf);*/
		/*if (inputfile.fail()) {
			cerr << "can't read " << buf << endl; return;
		}*/
		vector<cv::Point2f> tempPt;
		tempPt.push_back(cv::Point(shape.part(24).x(), shape.part(24).y())); // l eye (v 314)
		tempPt.push_back(cv::Point(shape.part(19).x(), shape.part(19).y()));// r eye (v 0)
		tempPt.push_back(cv::Point(shape.part(33).x(), shape.part(33).y()));//nose (v 1879)
		tempPt.push_back(cv::Point(shape.part(54).x(), shape.part(54).y()));// l mouth (v 1502)
		tempPt.push_back(cv::Point(shape.part(48).x(), shape.part(48).y())); // r mouth (v 695)
		tempPt.push_back(cv::Point(shape.part(16).x(), shape.part(16).y()));// l coutour (v 2011)
		tempPt.push_back(cv::Point(shape.part(0).x(), shape.part(0).y()));// r coutour (v 1138

		imagePoints.push_back(cv::Point(shape.part(42).x(), shape.part(42).y()));	 // l eye (v 314)
		imagePoints.push_back(cv::Point(shape.part(43).x(), shape.part(43).y()));	// l eye (v 314)
		imagePoints.push_back(cv::Point(shape.part(44).x(), shape.part(44).y()));	// l eye (v 314)
		imagePoints.push_back(cv::Point(shape.part(45).x(), shape.part(45).y()));	// l eye (v 314)
		imagePoints.push_back(cv::Point(shape.part(46).x(), shape.part(46).y()));	// l eye (v 314)
		imagePoints.push_back(cv::Point(shape.part(47).x(), shape.part(47).y()));	// l eye (v 314)


		imagePoints.push_back(cv::Point(shape.part(36).x(), shape.part(36).y()));	 // r eye (v 0)
		imagePoints.push_back(cv::Point(shape.part(37).x(), shape.part(37).y()));	 // r eye (v 0)
		imagePoints.push_back(cv::Point(shape.part(38).x(), shape.part(38).y()));	 // r eye (v 0)
		imagePoints.push_back(cv::Point(shape.part(39).x(), shape.part(39).y()));	 // r eye (v 0)
		imagePoints.push_back(cv::Point(shape.part(40).x(), shape.part(40).y()));	 // r eye (v 0)
		imagePoints.push_back(cv::Point(shape.part(41).x(), shape.part(41).y()));	 // r eye (v 0)


		imagePoints.push_back(cv::Point(shape.part(27).x(), shape.part(27).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(28).x(), shape.part(28).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(29).x(), shape.part(29).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(30).x(), shape.part(30).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(31).x(), shape.part(31).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(32).x(), shape.part(32).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(33).x(), shape.part(33).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(34).x(), shape.part(34).y()));	 //nose (v 1879)
		imagePoints.push_back(cv::Point(shape.part(35).x(), shape.part(35).y()));	 //nose (v 1879)


		imagePoints.push_back(cv::Point(shape.part(51).x(), shape.part(51).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(52).x(), shape.part(52).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(53).x(), shape.part(53).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(54).x(), shape.part(54).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(55).x(), shape.part(55).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(56).x(), shape.part(56).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(57).x(), shape.part(57).y()));	 // l mouth (v 1502)

		imagePoints.push_back(cv::Point(shape.part(62).x(), shape.part(62).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(63).x(), shape.part(63).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(64).x(), shape.part(64).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(65).x(), shape.part(65).y()));	 // l mouth (v 1502)
		imagePoints.push_back(cv::Point(shape.part(66).x(), shape.part(66).y()));	 // l mouth (v 1502)

		imagePoints.push_back(cv::Point(shape.part(50).x(), shape.part(50).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(49).x(), shape.part(49).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(48).x(), shape.part(48).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(59).x(), shape.part(59).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(58).x(), shape.part(58).y()));	 // r mouth (v 695)

		imagePoints.push_back(cv::Point(shape.part(61).x(), shape.part(61).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(60).x(), shape.part(60).y()));	 // r mouth (v 695)
		imagePoints.push_back(cv::Point(shape.part(67).x(), shape.part(67).y()));	 // r mouth (v 695)


		imagePoints.push_back(cv::Point(shape.part(16).x(), shape.part(16).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(15).x(), shape.part(15).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(14).x(), shape.part(14).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(13).x(), shape.part(13).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(12).x(), shape.part(12).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(11).x(), shape.part(11).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(10).x(), shape.part(10).y()));	 	// l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(9).x(), shape.part(9).y()));	 	    // l coutour (v 2011)
		imagePoints.push_back(cv::Point(shape.part(8).x(), shape.part(8).y()));	 	    // l coutour (v 2011)

		imagePoints.push_back(cv::Point(shape.part(0).x(), shape.part(0).y()));	 	// r coutour (v 1138
		imagePoints.push_back(cv::Point(shape.part(1).x(), shape.part(1).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(2).x(), shape.part(2).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(3).x(), shape.part(3).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(4).x(), shape.part(4).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(5).x(), shape.part(5).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(6).x(), shape.part(6).y()));	 	// r coutour (v 1138)
		imagePoints.push_back(cv::Point(shape.part(7).x(), shape.part(7).y()));	 	// r coutour (v 1138)

		imagePoints.push_back(cv::Point(shape.part(17).x(), shape.part(17).y()));	 	// r eyebrow(v 1138)
		imagePoints.push_back(cv::Point(shape.part(18).x(), shape.part(18).y()));	 	// r eyebrow(v 1138)
		imagePoints.push_back(cv::Point(shape.part(19).x(), shape.part(19).y()));	 	// r eyebrow(v 1138)
		imagePoints.push_back(cv::Point(shape.part(20).x(), shape.part(20).y()));	 	// r eyebrow(v 1138)
		imagePoints.push_back(cv::Point(shape.part(21).x(), shape.part(21).y()));	 	// r eyebrow(v 1138)

		imagePoints.push_back(cv::Point(shape.part(22).x(), shape.part(22).y()));	 	// l eyebrow (v 1138)
		imagePoints.push_back(cv::Point(shape.part(23).x(), shape.part(23).y()));	 	// l eyebrow (v 1138)
		imagePoints.push_back(cv::Point(shape.part(24).x(), shape.part(24).y()));	 	// l eyebrow (v 1138)
		imagePoints.push_back(cv::Point(shape.part(25).x(), shape.part(25).y()));	 	// l eyebrow (v 1138)
		imagePoints.push_back(cv::Point(shape.part(26).x(), shape.part(26).y()));	 	// l eyebrow (v 1138)
		double varianceForLandmark = calDistanceDiff(imagePoints, prevLandmarkPts);
		if (varianceForLandmark < 0.5) {
			std::swap(prevLandmarkPts, imagePoints);
		}
		double variance = calDistanceDiff(tempPt, prevImgPoint);
		if (variance < 0.5) {
			std::swap(prevImgPoint, tempPt);
		}
		cv::Mat ip = cv::Mat(imagePoints);
		imgPoint = imagePoints;

		/*sprintf(buf, "%sAngelina_Jolie/Angelina_Jolie_%04d.jpg", workingDir, counter);*/
		cv::Mat ip_temp = cv::Mat(tempPt);
		/*cv::Mat img = cv::imread(buf);*/

		//imgTex.set(dst); //TODO: what if different size??

		// paint 2D feature points
		//for (unsigned int i = 0; i < imagePoints.size(); i++) circle(dst, imagePoints[i], 1, cv::Scalar(255, 0, 255), CV_FILLED);
		for (unsigned int i = 0; i < shape.num_parts(); i++) {
			circle(dst, imgPoint[i], 1, cv::Scalar(255, 0, 255), CV_FILLED);
			//circle(dst, cv::Point2f(shape.part(i).x(),shape.part(i).y()), 1, cv::Scalar(255, 0, 255), CV_FILLED);
		}
		loadWithPoints(ip_temp, dst, imagePoints);
		cv::imshow("micky", dst);
		std::swap(prevLandmarkPts, imagePoints);
		swap(prevImgPoint, tempPt);
		//imgWithDrawing.set(dst);
	}
	counter = (counter + 1);
}

void loadWithPoints(cv::Mat& ip, cv::Mat& img, vector<cv::Point2f> imagePoints) {
	int max_d = MAX(img.rows, img.cols);
	camMatrix = (cv::Mat_<double>(3, 3) <<
		max_d, 0, img.cols / 2.0,
		0, max_d, img.rows / 2.0,
		0, 0, 1.0);
	//cout << "using cam matrix " << endl << camMatrix << endl;

	double _dc[] = { 0,0,0,0 };

	//input 3D point, 2D landmark to get the head pose estimation with:
	//1. Head position (tvec)
	//2. Head rotation (face left, face right, face top , face bottom) (rvec)
	solvePnP(op_min, ip, camMatrix, cv::Mat(1, 4, CV_64FC1, _dc), rvec, tvec, false, CV_EPNP);
	cv::Mat rotM(3, 3, CV_64FC1, rot);
	// convert the vector to matrix
	Rodrigues(rvec, rotM);
	double* _r = rotM.ptr<double>();
	//printf("rotation mat: \n %.3f %.3f %.3f\n%.3f %.3f %.3f\n%.3f %.3f %.3f\n",
	//	_r[0], _r[1], _r[2], _r[3], _r[4], _r[5], _r[6], _r[7], _r[8]);

	//printf("trans vec: \n %.3f %.3f %.3f\n", tv[0], tv[1], tv[2]);

	double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
					  _r[3],_r[4],_r[5],tv[1],
					  _r[6],_r[7],_r[8],tv[2] };

	cv::Matx34d P(_pm);
	cv::Mat KP = camMatrix * cv::Mat(P);
	////cout << "KP " << endl << KP << endl;

	//reproject object points - check validity of found projection matrix
	//op - vertex point in 3D mesh, ip - 2D landmark of the images
	for (int i = 0; i < op_min.rows; i++) {

		cv::Mat_<double> X = (cv::Mat_<double>(4, 1) << op_min.at<float>(i, 0), op_min.at<float>(i, 1), op_min.at<float>(i, 2), 1.0);
		//cout << "object point " << X << endl;
		cv::Mat_<double> opt_p = KP * X;
		//cout << "object point " << opt_p << endl;
		double TwoDopX = opt_p(0) / opt_p(2);
		double TwoDopY = opt_p(1) / opt_p(2);
		/*double imageX = imagePoints[i].x;
		double imageY = imagePoints[i].y;*/
		//cout << " opX :" << TwoDopX << " opY :" << TwoDopY << endl;
		/*	double eimX = findEimX(opt_p(0), imageX, opt_p(2),max_d, img.cols / 2.0);
			double eimY = findEimX(opt_p(1), imageY, opt_p(2), max_d, img.rows / 2.0);
			cout << "eimX :" << eimX << " eimY :" << eimY << endl;
			cout << "ein " << eimX + eimY << endl;
			cout << " opX :" << TwoDopX << " opY :" << TwoDopY << endl;
			cout << " imgX :" << imageX << " imgY :" << imageY << endl;*/
			//double diffX = imageX - TwoDopX;
			//double diffY = imageY - TwoDopX;
			/*double dist = sqrt(pow((imageX - TwoDopX), 2) + pow((imageY - TwoDopY), 2));
			cout << "distance :" << dist << " X:" << TwoDopX << " Y:" << TwoDopY  << endl;*/
		cv::Point2f opt_p_img(TwoDopX, TwoDopY);
		//cout << "object point reproj " << opt_p_img << endl;
		circle(img, opt_p_img, 4, cv::Scalar(0, 0, 255), 1);
	}
	//transpose for following the head rotation is same as the web cam face rotation otherwise it will move with opposite direction
	rotM = rotM.t();// transpose to conform with majorness of opengl matrix
}


int main(int argc, char** argv)
{
	cout << "read blendshape" << endl;
	//calculated_average_face = glmReadOBJ("shape_0_1.obj");
	average_face = glmReadOBJ("shape_0_1.obj");
	/*for (int j = 0; j < numOfuser; j++) {
		vector<GLMmodel*> model;
		for (int i = 0; i < numOfBlendshape; i++) {
			char *theString = "FWH/";
			stringstream ss;
			ss << theString << j << "/" << "shape_" << i << ".obj";
			char *cstr = new char[ss.str().length() + 1];
			strcpy(cstr, ss.str().c_str());
			model.push_back(glmReadOBJ(cstr));
			if (j == 0) {
				char *theString2 = "AverageFace/";
				stringstream sss;
				sss << theString2 << "shape_" << i << ".obj";
				char *cstrr = new char[sss.str().length() + 1];
				strcpy(cstrr, sss.str().c_str());
				empty_face.push_back(glmReadOBJ(cstrr));
			}
		}
		rank3Tensor.push_back(model);
	}*/
	for (int i = 0; i < 7; i++) {
		prevImgPoint.push_back(cv::Point2f(0, 0));
	}
	for (int i = 0; i < 68; i++) {
		prevLandmarkPts.push_back(cv::Point2f(0, 0));
	}
	vector<BlendShape> shapes;

	const int nShapes = 2;			// 150 identity
	const int nExprs = 47;				// 46 expressions + 1 neutral
	const int nVerts = 11510;			// 11510 vertices for each mesh

	const string path = "C:\\Users\\ccha5314\\Desktop\\FacewareHouse_allData\\";
	const string foldername = "Tester_";
	const string bsfolder = "Blendshape";
	const string filename = "shape.bs";

	shapes.resize(nShapes);
	for (int i = 0; i < nShapes; i++) {
		stringstream ss;
		ss << path << foldername << (i + 1) << "\\" << bsfolder + "\\" + filename;

		shapes[i].read(ss.str());
	}
	int nCoords = nVerts * 3;

	// create an order 3 tensor for the blend shapes	
	Tensor3<double> t(nShapes, nExprs, nCoords);
	//cout << nShapes << " " << nExprs << " " << nCoords << endl;
	// fill in the data
	for (int i = 0; i < shapes.size(); i++) {
		const BlendShape& bsi = shapes[i];
		//cout << bsi.expressionCount() << endl;
		for (int j = 0; j < bsi.expressionCount(); j++) {
			const BlendShape::shape_t& bsij = bsi.expression(j);

			for (int k = 0, cidx = 0; k < nVerts; k++, cidx += 3) {
				const BlendShape::vert_t& v = bsij[k];

				t(i, j, cidx) = v.x;
				t(i, j, cidx + 1) = v.y;
				t(i, j, cidx + 2) = v.z;
			}
		}
	}
	initTensor = t;
	int ms[2] = { 0, 1 };		// only the first two modes
	int ds[2] = { numOfuser, 47 };	// pick 50 for identity and 25 for expression
	vector<int> modes(ms, ms + 2);
	vector<int> dims(ds, ds + 2);
	auto comp2 = t.svd(modes, dims);
	vector<double> wexptest, widtest;
	widtest.push_back(1.0);
	widtest.push_back(0);
	wexptest.push_back(1);
	for (int i = 1; i < numOfBlendshape; i++) {
		wexptest.push_back(0);
	}
	auto result = initTensor.modeProduct(widtest, 0).modeProduct(wexptest, 1);

	coreTensor = std::get<0>(comp2);
	UIdAndUExp = std::get<1>(comp2);
	/*dataTensor = coreTensor.modeProduct(UIdAndUExp[0], 0).modeProduct(UIdAndUExp[1], 1);*/
	vector<cv::Point3d> modelPoint;
	vector<int> vertexPoint;
	//tempPt.push_back(cv::Point(shape.part(45).x(), shape.part(45).y())); // l eye (v 314)
	//tempPt.push_back(cv::Point(shape.part(37).x(), shape.part(37).y()));// r eye (v 0)
	//tempPt.push_back(cv::Point(shape.part(33).x(), shape.part(33).y()));//nose (v 1879)
	//tempPt.push_back(cv::Point(shape.part(54).x(), shape.part(54).y()));// l mouth (v 1502)
	//tempPt.push_back(cv::Point(shape.part(48).x(), shape.part(48).y())); // r mouth (v 695)
	//tempPt.push_back(cv::Point(shape.part(16).x(), shape.part(16).y()));// l coutour (v 2011)
	//tempPt.push_back(cv::Point(shape.part(0).x(), shape.part(0).y()));// r coutour (v 1138
	vertexPoint.push_back(7347); // left eye
	vertexPoint.push_back(4451); //right eye
	vertexPoint.push_back(6328); // nose
	vertexPoint.push_back(1594);	 // l mouth (v 6601) or line 6605
	vertexPoint.push_back(187);		 // r mouth (v 728) or line 732 
	vertexPoint.push_back(6878);	 // l contour (v 8699)
	vertexPoint.push_back(685);	 // r contour (v 3484)

	for (int i = 0; i < vertexPoint.size(); i++) {
		double x = roundFloat(average_face->vertices[3 * vertexPoint[i]]);
		double y = roundFloat(average_face->vertices[3 * vertexPoint[i] + 1]);
		double z = roundFloat(average_face->vertices[3 * vertexPoint[i] + 2]);
		modelPoint.push_back(cv::Point3d(
			x,
			y,
			z
		));
	}
	op_min = cv::Mat(modelPoint);
	cout << "finished blendshape :" << endl;


	dlib::deserialize("shape_predictor_68_face_landmarks.dat") >> sp;
	detector = dlib::get_frontal_face_detector();
	//op = cv::Mat(modelPoints);
	//cv::Scalar m = mean(cv::Mat(modelPoints));
	//op = op - m;

	//assert(norm(mean(op)) < 1e-05); //make sure is centered

	//op = op + Scalar(0,0,25);

	//cout << "model points " << op << endl;
	rvec = cv::Mat(rv);
	double _d[9] = { 1,0,0,
					0,-1,0,
					0,0,-1 };
	Rodrigues(cv::Mat(3, 3, CV_64FC1, _d), rvec);
	tv[0] = 0; tv[1] = 0; tv[2] = 1;
	tvec = cv::Mat(tv);

	camMatrix = cv::Mat(3, 3, CV_64FC1);

	init_opengl(argc, argv); // get GL context, for loading textures

	// prepare OpenCV-OpenGL images
	imgTex = MakeOpenCVGLTexture(cv::Mat());
	imgWithDrawing = MakeOpenCVGLTexture(cv::Mat());

	loadNext();

	start_opengl();

	return 0;
}


