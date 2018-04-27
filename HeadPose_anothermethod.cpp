//// HeadPose.cpp : Defines the entry point for the console application.
////
////#include <glut.h>
//#include <glew.h>
//#include "stdafx.h"
//
//#include "opencv\cv.h"
//#include "opencv\highgui.h"
////using namespace cv;
//
//#include <vector>
//#include <iostream>
//#include <fstream>
//#include <dlib\optimization.h>
//#include <dlib/opencv.h>
//#include <dlib/image_processing/frontal_face_detector.h>
//#include <dlib/image_processing/render_face_detections.h>
//#include <dlib/image_processing.h>
//#include <dlib/gui_widgets.h>
//#include <dlib/image_io.h>
//using namespace std;
////#if defined(__APPLE__)
////#  include <OpenGL/gl.h>
////#  include <OpenGL/glu.h>
////#elif defined(__linux__) || defined(__MINGW32__) || defined(WIN32)
////#  include <GL/gl.h>
////#  include <GL/glu.h>
////#else
////#include <gl\GL.h>
////#include <gl\GLU.h>
//
////#endif
//
//#include <Eigen/Dense>
//#include <Eigen/Eigen>
//#include "glm.h"
//#include "OGL_OCV_common.h"
//#include "BlendShape.h"
//#include "Tensor.hpp"
//#include "nnls.h"
//using namespace Eigen;
//
//void loadNext();
//void loadWithPoints(cv::Mat& ip, cv::Mat& img, vector<cv::Point2f> imagepoints);
//
//const GLfloat light_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
//const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
//const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
//const GLfloat light_position[] = { 0.0f, 0.0f, 1.0f, 0.0f };
//
//const GLfloat mat_ambient[] = { 0.7f, 0.7f, 0.7f, 1.0f };
//const GLfloat mat_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
//const GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
//const GLfloat high_shininess[] = { 100.0f };
//
//double rot[9] = { 0 };
//GLuint textureID;
//cv::Mat backPxls;
//vector<double> rv(3), tv(3);
//cv::Mat rvec(rv), tvec(tv);
//cv::Mat camMatrix;
//OpenCVGLTexture imgTex, imgWithDrawing;
//cv::VideoCapture webcam(0);
//vector<GLMmodel*> empty_face;
//GLMmodel* average_face;
//vector<GLMmodel*> blendshapes;
//vector<vector<GLMmodel*>> rank3Tensor;
//vector<vector<double>> wIdX, wIdY, wIdZ;
//vector<vector<double>>  wExpX, wExpY, wExpZ;
//dlib::shape_predictor sp;
//dlib::frontal_face_detector detector;
//int chaningRadio = 1;
//int faceDesign = 0;
//cv::Mat op;
//cv::Mat op_min;
//vector<cv::Point2f> imgPoint;
//int numOfuser = 2;
//int numOfBlendshape = 47;
//int numOfVertex = 11510;
//Tensor3<double> coreTensor;
//vector<Tensor2<double>> widAndWExp;
//Tensor3<double> tempWExpTensor;
//Tensor3<double> tempWIdTensor;
//
//vector<double> findWid(vector<int> vertexPoint) {
//	double rotss[9] = { 0 };
//	cv::Mat rotMss(3, 3, CV_64FC1, rotss);
//	// convert the vector to matrix
//	Rodrigues(rvec, rotMss);
//	double* _r = rotMss.ptr<double>();
//	double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
//		_r[3],_r[4],_r[5],tv[1],
//		_r[6],_r[7],_r[8],tv[2] };
//	cv::Matx34d P(_pm);
//	cv::Mat KP = camMatrix * cv::Mat(P);
//	//double y0 = orignialFormula(wid, wexp, vertexPoint, KP);
//	int idx = 0;
//	vector<double> widList;
//	if (imgPoint.size() != 0) {
//		/*arma::vec b(imgPoint.size());
//		arma::mat C(imgPoint.size(), numOfBlendshape);*/
//		int numLM = 68;
//		MatrixXd A(numLM * 2, numOfuser);
//		VectorXd b(numLM * 2);
//		//for (int i = 0; i < numOfuser; i++) {
//
//		for (int j = 0; j < numOfuser; j++) {
//			for (int k = 0; k < imgPoint.size(); k++) {
//				double pX = rank3Tensor[j][0]->vertices[3 * vertexPoint[k]];
//				double pY = rank3Tensor[j][0]->vertices[3 * vertexPoint[k] + 1];
//				double pZ = rank3Tensor[j][0]->vertices[3 * vertexPoint[k] + 2];
//				/*double pX = tempWIdTensor(j, 0, vertexPoint[k]);
//				double pY = tempWIdTensor(j, 0, vertexPoint[k] + 1);
//				double pZ = tempWIdTensor(j, 0, vertexPoint[k] + 2);*/
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
//
//			}
//		}
//		VectorXd x(numOfuser);
//		//Solving Ax = b issues with non-negative least square
//		NNLS<MatrixXd>::solve(A, b, x);
//		//VectorXd solution = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
//		double totalSum = 0;
//		for (int i = 0; i < x.size(); i++) {
//			widList.push_back(x(i));
//			//widList.push_back(solution(i));
//			//totalSum += fabs(solution(i));
//		}
//		//for (int i = 0; i < x.size(); i++) {
//		//	widList.push_back(x(i));
//		//	/*widList.push_back(fabs(solution(i)));
//		//	totalSum += fabs(solution(i));*/
//		//}
//		/*for (int i = 0; i < solution.size(); i++) {
//		widList[i] = widList[i] / totalSum;
//		}*/
//	}
//	return widList;
//}
//
//vector<double> findWExp(vector<int> vertexPoint) {
//	double rotss[9] = { 0 };
//	cv::Mat rotMss(3, 3, CV_64FC1, rotss);
//	// convert the vector to matrix
//	Rodrigues(rvec, rotMss);
//	double* _r = rotMss.ptr<double>();
//	double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
//		_r[3],_r[4],_r[5],tv[1],
//		_r[6],_r[7],_r[8],tv[2] };
//	cv::Matx34d P(_pm);
//	cv::Mat KP = camMatrix * cv::Mat(P);
//	int idx = 0;
//	vector<double> wexpList;
//	if (imgPoint.size() != 0) {
//		int numBS = 47;
//		int numLM = 68;
//		MatrixXd A(numLM * 2, numBS);
//		VectorXd b(numLM * 2);
//		//for (int i = 0; i < numOfuser; i++) {
//		for (int j = 0; j < numBS; j++) {
//			for (int k = 0; k < imgPoint.size(); k++) {
//				double pX = rank3Tensor[0][j]->vertices[3 * vertexPoint[k]];
//				double pY = rank3Tensor[0][j]->vertices[3 * vertexPoint[k] + 1];
//				double pZ = rank3Tensor[0][j]->vertices[3 * vertexPoint[k] + 2];
//				//double pX = tempWExpTensor(0, j, vertexPoint[k]);// -tempWExpTensor(0, 0, vertexPoint[k]);
//				//double pY = tempWExpTensor(0, j, vertexPoint[k] + 1);// -tempWExpTensor(0, 0, vertexPoint[k] + 1);
//				//double pZ = tempWExpTensor(0, j, vertexPoint[k] + 2);// -tempWExpTensor(0, 0, vertexPoint[k] + 2);
//				cv::Mat_<double> X = (cv::Mat_<double>(4, 1)
//					<< pX,
//					pY,
//					pZ, 1.0);//convert point to vector
//				cv::Mat_<double> opt_m = KP * X; // KP = camera Matrix * Pose Matrix
//				double finalX = opt_m(0) / opt_m(2);
//				double finalY = opt_m(1) / opt_m(2);
//				A(k, j) = finalX;
//				A(k + vertexPoint.size(), j) = finalY;
//				//tempA[k + vertexPoint.size()][j] = finalY;
//				if (j == 0) {
//					b(k) = imgPoint[k].x;
//					b(k + vertexPoint.size()) = imgPoint[k].y;
//					//	tempB[k + vertexPoint.size()] = imgPoint[k].y;;
//				}
//			}
//		}
//
//		VectorXd x(numBS);
//		//Solving Ax = b issues with non-negative least square
//		NNLS<MatrixXd>::solve(A, b, x);
//		//VectorXd solution = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
//		double totalSum = 0;
//		for (int i = 0; i < x.size(); i++) {
//			//wexpList.push_back(solution(i));
//			//dlib::clamp_function(x(i), 0, 1);
//			wexpList.push_back(x(i));
//			//totalSum += fabs(solution(i));
//		}
//		/*for (int i = 0; i < solution.size(); i++) {
//		wexpList[i] = wexpList[i] / totalSum;
//		}*/
//	}
//	return wexpList;
//}
//
//double roundFloat(double var)
//{
//	// we use array of chars to store number
//	// as a string.
//	char str[80];
//
//	// Print in string the value of var 
//	// with two decimal point
//	sprintf(str, "%.6lf", var);
//	double lf2 = atof(str);
//	// scan string value in var 
//	//sscanf(str, "%f", &var);
//
//	return lf2;
//}
//
//void saveOpenGLBuffer() {
//	static unsigned int opengl_buffer_num = 0;
//
//	int vPort[4]; glGetIntegerv(GL_VIEWPORT, vPort);
//	cv::Mat_<cv::Vec3b> opengl_image(vPort[3], vPort[2]);
//	{
//		cv::Mat_<cv::Vec4b> opengl_image_4b(vPort[3], vPort[2]);
//		glReadPixels(0, 0, vPort[2], vPort[3], GL_BGRA, GL_UNSIGNED_BYTE, opengl_image_4b.data);
//		cv::flip(opengl_image_4b, opengl_image_4b, 0);
//		mixChannels(&opengl_image_4b, 1, &opengl_image, 1, &(cv::Vec6i(0, 0, 1, 1, 2, 2)[0]), 3);
//	}
//	stringstream ss; ss << "opengl_buffer_" << opengl_buffer_num++ << ".jpg";
//	imwrite(ss.str(), opengl_image);
//}
//void resize(int width, int height)
//{
//	const float ar = (float)width / (float)height;
//
//	glViewport(0, 0, width, height);
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	//glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
//	gluPerspective(47, 1.0, 0.01, 1000.0);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//}
//int __w = 640, __h = 480;
//void key(unsigned char key, int x, int y)
//{
//
//	switch (key)
//	{
//	case 27:
//	case 'Q':
//	case 'q':
//		exit(0);
//		break;
//	case 'w':
//	case 'W':
//		__w++;
//		__w = __w % 251;
//		break;
//	case 'h':
//	case 'H':
//		__h++;
//		__h = __h % 250;
//		break;
//	case ' ':
//		saveOpenGLBuffer();
//		//loadNext();
//		break;
//	case 'z':
//		break;
//	case 'x':
//		break;
//	case 'c':
//		//printf("current weight point : %f %f %f ", weightX, weightY, weightZ);
//		break;
//		//case 'x':
//		//	weightX--;
//		//	break;
//	default:
//		break;
//	}
//
//	glutPostRedisplay();
//}
//void idle(void)
//{
//	//loadNext();
//	glutPostRedisplay();
//}
//void myGLinit() {
//	//    glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;
//
//	glEnable(GL_CULL_FACE);
//	glCullFace(GL_BACK);
//
//
//	glShadeModel(GL_SMOOTH);
//
//	glEnable(GL_DEPTH_TEST);
//	glDepthFunc(GL_LEQUAL);
//
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//
//	glEnable(GL_LIGHT0);
//	glEnable(GL_NORMALIZE);
//	glEnable(GL_COLOR_MATERIAL);
//	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
//
//	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
//	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
//	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
//	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//
//	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
//	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
//
//	glEnable(GL_LIGHTING);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//
//}
//void drawAxes() {
//
//	//Z = red
//	glPushMatrix();
//	glRotated(180, 0, 1, 0);
//	glColor4d(1, 0, 0, 0.5);
//	//	glutSolidCylinder(0.05,1,15,20);
//	glBegin(GL_LINES);
//	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
//	glEnd();
//	glTranslated(0, 0, 1);
//	glScaled(.1, .1, .1);
//	glutSolidTetrahedron();
//	glPopMatrix();
//
//	//Y = green
//	glPushMatrix();
//	glRotated(-90, 1, 0, 0);
//	glColor4d(0, 1, 0, 0.5);
//	//	glutSolidCylinder(0.05,1,15,20);
//	glBegin(GL_LINES);
//	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
//	glEnd();
//	glTranslated(0, 0, 1);
//	glScaled(.1, .1, .1);
//	glutSolidTetrahedron();
//	glPopMatrix();
//
//	//X = blue
//	glPushMatrix();
//	glRotated(-90, 0, 1, 0);
//	glColor4d(0, 0, 1, 0.5);
//	//	glutSolidCylinder(0.05,1,15,20);
//	glBegin(GL_LINES);
//	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
//	glEnd();
//	glTranslated(0, 0, 1);
//	glScaled(.1, .1, .1);
//	glutSolidTetrahedron();
//	glPopMatrix();
//}
//
//void display(void)
//{
//	// draw the image in the back
//	int vPort[4]; glGetIntegerv(GL_VIEWPORT, vPort);
//	glEnable2D();
//	//drawOpenCVImageInGL(imgTex);
//	//glTranslated(vPort[2] / 2.0, 0, 1);
//	//drawOpenCVImageInGL(imgWithDrawing);
//	glDisable2D();
//
//	glClear(GL_DEPTH_BUFFER_BIT); // we want to draw stuff over the image
//
//								  // draw only on left part
//	glViewport(0, 0, vPort[2], vPort[3]);
//
//	glPushMatrix();
//
//	gluLookAt(0, 0, 0, 0, 0, 1, 0, -1, 0);
//
//
//	vector<cv::Point3f > modelPoints;
//	vector<int> vertexPoint;
//	vertexPoint.push_back(8638);	 // l eye (v 8083)
//	vertexPoint.push_back(8091);
//	vertexPoint.push_back(8467);
//	vertexPoint.push_back(8068);
//	vertexPoint.push_back(8928);
//	vertexPoint.push_back(8874);
//
//	vertexPoint.push_back(2325);	 // r eye (v 2343)
//	vertexPoint.push_back(2538);
//	vertexPoint.push_back(2710);
//	vertexPoint.push_back(3350);
//	vertexPoint.push_back(3216);
//	vertexPoint.push_back(3158);
//
//	vertexPoint.push_back(8710);	 //nose (v 1264)
//	vertexPoint.push_back(11151);
//	vertexPoint.push_back(1235);
//	vertexPoint.push_back(1197);
//	vertexPoint.push_back(1549);
//	vertexPoint.push_back(1553);
//	vertexPoint.push_back(1263);
//	vertexPoint.push_back(7328);
//	vertexPoint.push_back(7321);
//
//	vertexPoint.push_back(7381);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6637);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6683);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6601);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6747);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6565);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(6523);	 // l mouth (v 6601) or line 6605
//
//	vertexPoint.push_back(6551);
//	vertexPoint.push_back(6632);
//	vertexPoint.push_back(6730);
//	vertexPoint.push_back(6531);
//	vertexPoint.push_back(6522);
//
//	vertexPoint.push_back(832);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(879);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(728);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(754);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(866);		 // r mouth (v 728) or line 732 
//
//	vertexPoint.push_back(825);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(893);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(876);		 // r mouth (v 728) or line 732 
//
//	vertexPoint.push_back(9374);	 // l contour (v 8699)
//	vertexPoint.push_back(11176);	 // l contour (v 8699)
//	vertexPoint.push_back(7560);	 // l contour (v 8699)
//	vertexPoint.push_back(7900);	 // l contour (v 8699)
//	vertexPoint.push_back(7646);	 // l contour (v 8699)
//	vertexPoint.push_back(7718);	 // l contour (v 8699)
//	vertexPoint.push_back(7836);	 // l contour (v 8699)
//	vertexPoint.push_back(7579);	 // l contour (v 8699)
//	vertexPoint.push_back(1690);	 // l contour (v 8699)
//
//	vertexPoint.push_back(3917);	 // r contour (v 3484)
//	vertexPoint.push_back(5490);	 // r contour (v 3484)
//	vertexPoint.push_back(2245);	 // r contour (v 3484)
//	vertexPoint.push_back(1776);	 // r contour (v 3484)
//	vertexPoint.push_back(1805);	 // r contour (v 3484)
//	vertexPoint.push_back(2079);	 // r contour (v 3484)
//	vertexPoint.push_back(1671);	 // r contour (v 3484)
//	vertexPoint.push_back(1824);	 // r contour (v 3484)
//
//	vertexPoint.push_back(2712);	 // r eyebrow(v 1138)
//	vertexPoint.push_back(2495);	 // r eyebrow(v 1138)
//	vertexPoint.push_back(2700);	 // r eyebrow(v 1138)
//	vertexPoint.push_back(2698);	 // r eyebrow(v 1138)
//	vertexPoint.push_back(2481);	 // r eyebrow(v 1138)
//
//	vertexPoint.push_back(8814);	 // l eyebrow(v 1138)
//	vertexPoint.push_back(8815);	 // l eyebrow(v 1138)
//	vertexPoint.push_back(8830);	 // l eyebrow(v 1138)
//	vertexPoint.push_back(8773);	 // l eyebrow(v 1138)
//	vertexPoint.push_back(8776);	 // l eyebrow(v 1138)
//
//									 //modelPoints.push_back(cv::Point3f(0.294959, 0.493100, 0.389258));	// l eye (v 8083)
//									 //modelPoints.push_back(cv::Point3f(-0.263214, 0.495387, 0.373535));	// r eye (v 2343)
//									 //modelPoints.push_back(cv::Point3f(0.029220, 0.188548, 0.673340));	//nose (v 1264)
//									 //modelPoints.push_back(cv::Point3f(0.291080, -0.127212, 0.454637));	// l mouth (v 6601) or line 6605
//									 //modelPoints.push_back(cv::Point3f(-0.171941, -0.121142, 0.518103));	// r mouth (v 728) or line 732 
//									 //modelPoints.push_back(cv::Point3f(0.857392, 0.459262, -0.547478));	// l ear (v 8699)
//									 //modelPoints.push_back(cv::Point3f(-0.853804, 0.425016, -0.562550));		// r ear (v 3484)
//	for (int i = 0; i < vertexPoint.size(); i++) {
//		double x = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i]]);
//		double y = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i] + 1]);
//		double z = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i] + 2]);
//		modelPoints.push_back(cv::Point3d(
//			x,
//			y,
//			z
//		));
//	}
//	op = cv::Mat(modelPoints);
//	loadNext();
//
//	// put the object in the right position in space - make sure that the 3D head will move with the web cam detected face
//	cv::Vec3d tvv(tv[0], tv[1], tv[2]);
//	glTranslated(tvv[0], tvv[1], tvv[2]);
//	//glTranslated(0, 0, 5);
//	// control the head rotation
//	double _d[16] = { rot[0],rot[1],rot[2],0,
//		rot[3],rot[4],rot[5],0,
//		rot[6],rot[7],rot[8],0,
//		0,	   0,	  0		,1 };
//
//	glMultMatrixd(_d);
//	GLMmodel* original_obj = new GLMmodel();
//	//Store the neutral face vertex for reset the neutral face after drawing
//	vector<GLfloat> originalVertex;
//	vector<double> wexp = findWExp(vertexPoint);
//	vector<double> wid = findWid(vertexPoint);
//	//double sum = 0;
//	//for (int i = 0; i < wid.size(); i++) {
//	//	sum += wid[i];
//	//}
//	//cout << sum << endl;
//	//cout << sum << endl;
//	//cout << opt_wid << endl;
//	for (int i = 1; i <= numOfVertex; i++) {
//		originalVertex.push_back(empty_face[0]->vertices[3 * i]);
//		originalVertex.push_back(empty_face[0]->vertices[3 * i + 1]);
//		originalVertex.push_back(empty_face[0]->vertices[3 * i + 2]);
//		empty_face[0]->vertices[3 * i] = 0;
//		empty_face[0]->vertices[3 * i + 1] = 0;
//		empty_face[0]->vertices[3 * i + 2] = 0;
//	}
//
//	if (wexp.size() > 0) {
//		for (int i = 0; i < numOfuser; i++) {
//			for (int j = 0; j < numOfBlendshape; j++) {
//				for (int k = 0; k < numOfVertex; k++) {
//					double faceX = rank3Tensor[i][j]->vertices[3 * k];// -rank3Tensor[i][0]->vertices[3 * k];
//					double faceY = rank3Tensor[i][j]->vertices[3 * k + 1];// -rank3Tensor[i][0]->vertices[3 * k + 1];
//					double faceZ = rank3Tensor[i][j]->vertices[3 * k + 2];// -rank3Tensor[i][0]->vertices[3 * k + 2];
//																		  /*double faceX = coreTensor(i, j, k);
//																		  double faceY = coreTensor(i, j, k + 1);
//																		  double faceZ = coreTensor(i, j, k + 2);
//																		  double vertexX = faceX   * wid[i] * wexp[j];
//																		  double vertexY = faceY   * wid[i] * wexp[j];
//																		  double vertexZ = faceZ   * wid[i] * wexp[j];*/
//																		  /*cout << vertexX << endl;
//																		  cout << vertexX << endl;
//																		  cout << vertexX << endl;*/
//					empty_face[0]->vertices[3 * k] += faceX * wid[i] * wexp[j];
//					empty_face[0]->vertices[3 * k + 1] += faceY * wid[i] * wexp[j];
//					empty_face[0]->vertices[3 * k + 2] += faceZ * wid[i] * wexp[j];
//				}
//			}
//		}
//	}
//	// draw the 3D head model
//	glColor4f(1, 1, 1, 0.75);
//	glmDraw(empty_face[0], GLM_SMOOTH);
//
//	//Reset the neutral face to the orignial
//	for (int i = 0; i < originalVertex.size(); i++) {
//		empty_face[0]->vertices[3 + i] = originalVertex[i];
//	}
//	//----------Axes
//	//glScaled(.5, .5,.5);
//	//drawAxes();
//	//----------End axes
//
//	glPopMatrix();
//	// restore to looking at complete viewport
//	glutSwapBuffers();
//}
//
//void init_opengl(int argc, char** argv) {
//	glutInitWindowSize(640, 640);
//	glutInitWindowPosition(40, 40);
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH); // | GLUT_MULTISAMPLE
//	glutCreateWindow("head pose");
//
//	myGLinit();
//
//	glutReshapeFunc(resize);
//	glutDisplayFunc(display);
//	glutKeyboardFunc(key);
//	//glutSpecialFunc(special);
//	glutIdleFunc(idle);
//}
//
//int start_opengl() {
//
//	glutMainLoop();
//
//	return 1;
//}
//
//
//void loadNext() {
//	static int counter = 1;
//
//	//printf("load %d\n", counter);
//
//	const char* workingDir = "./";
//	//Micky 
//
//	cv::Mat temp;
//	//temp = cv::imread("Stephen-640x480.jpg");
//	webcam >> temp;
//	cv::Mat dst;
//	cv::resize(temp, dst, cv::Size(352, 240));
//	dlib::cv_image <dlib::bgr_pixel > cimg(dst);
//	std::vector<dlib::rectangle> dets = detector(cimg);
//	if (dets.size() != 0) {
//		dlib::full_object_detection shape = sp(cimg, dets[0]);
//		//char buf[256] = { 0 };
//		//sprintf(buf, "%sAngelina_Jolie/Angelina_Jolie_%04d.txt", workingDir, counter);
//
//		vector<cv::Point2f> imagePoints;
//		/*ifstream inputfile(buf);*/
//		/*if (inputfile.fail()) {
//		cerr << "can't read " << buf << endl; return;
//		}*/
//		vector<cv::Point2d> tempPt;
//		tempPt.push_back(cv::Point(shape.part(45).x(), shape.part(45).y()));
//		tempPt.push_back(cv::Point(shape.part(37).x(), shape.part(37).y()));
//		tempPt.push_back(cv::Point(shape.part(33).x(), shape.part(33).y()));
//		tempPt.push_back(cv::Point(shape.part(54).x(), shape.part(54).y()));
//		tempPt.push_back(cv::Point(shape.part(48).x(), shape.part(48).y()));
//		tempPt.push_back(cv::Point(shape.part(16).x(), shape.part(16).y()));
//		tempPt.push_back(cv::Point(shape.part(0).x(), shape.part(0).y()));
//
//		imagePoints.push_back(cv::Point(shape.part(44).x(), shape.part(44).y()));	 // l eye (v 314)
//		imagePoints.push_back(cv::Point(shape.part(45).x(), shape.part(45).y()));	// l eye (v 314)
//		imagePoints.push_back(cv::Point(shape.part(46).x(), shape.part(46).y()));	// l eye (v 314)
//		imagePoints.push_back(cv::Point(shape.part(47).x(), shape.part(47).y()));	// l eye (v 314)
//		imagePoints.push_back(cv::Point(shape.part(48).x(), shape.part(48).y()));	// l eye (v 314)
//		imagePoints.push_back(cv::Point(shape.part(49).x(), shape.part(49).y()));	// l eye (v 314)
//
//
//		imagePoints.push_back(cv::Point(shape.part(36).x(), shape.part(36).y()));	 // r eye (v 0)
//		imagePoints.push_back(cv::Point(shape.part(37).x(), shape.part(37).y()));	 // r eye (v 0)
//		imagePoints.push_back(cv::Point(shape.part(38).x(), shape.part(38).y()));	 // r eye (v 0)
//		imagePoints.push_back(cv::Point(shape.part(39).x(), shape.part(39).y()));	 // r eye (v 0)
//		imagePoints.push_back(cv::Point(shape.part(40).x(), shape.part(40).y()));	 // r eye (v 0)
//		imagePoints.push_back(cv::Point(shape.part(41).x(), shape.part(41).y()));	 // r eye (v 0)
//
//
//		imagePoints.push_back(cv::Point(shape.part(27).x(), shape.part(27).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(28).x(), shape.part(28).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(29).x(), shape.part(29).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(30).x(), shape.part(30).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(31).x(), shape.part(31).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(32).x(), shape.part(32).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(33).x(), shape.part(33).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(34).x(), shape.part(34).y()));	 //nose (v 1879)
//		imagePoints.push_back(cv::Point(shape.part(35).x(), shape.part(35).y()));	 //nose (v 1879)
//
//
//		imagePoints.push_back(cv::Point(shape.part(51).x(), shape.part(51).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(52).x(), shape.part(52).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(53).x(), shape.part(53).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(54).x(), shape.part(54).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(55).x(), shape.part(55).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(56).x(), shape.part(56).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(57).x(), shape.part(57).y()));	 // l mouth (v 1502)
//
//		imagePoints.push_back(cv::Point(shape.part(62).x(), shape.part(62).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(63).x(), shape.part(63).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(64).x(), shape.part(64).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(65).x(), shape.part(65).y()));	 // l mouth (v 1502)
//		imagePoints.push_back(cv::Point(shape.part(66).x(), shape.part(66).y()));	 // l mouth (v 1502)
//
//		imagePoints.push_back(cv::Point(shape.part(50).x(), shape.part(50).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(49).x(), shape.part(49).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(48).x(), shape.part(48).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(59).x(), shape.part(59).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(58).x(), shape.part(58).y()));	 // r mouth (v 695)
//
//		imagePoints.push_back(cv::Point(shape.part(61).x(), shape.part(61).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(60).x(), shape.part(60).y()));	 // r mouth (v 695)
//		imagePoints.push_back(cv::Point(shape.part(67).x(), shape.part(67).y()));	 // r mouth (v 695)
//
//
//		imagePoints.push_back(cv::Point(shape.part(16).x(), shape.part(16).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(15).x(), shape.part(15).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(14).x(), shape.part(14).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(13).x(), shape.part(13).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(12).x(), shape.part(12).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(11).x(), shape.part(11).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(10).x(), shape.part(10).y()));	 	// l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(9).x(), shape.part(9).y()));	 	    // l coutour (v 2011)
//		imagePoints.push_back(cv::Point(shape.part(8).x(), shape.part(8).y()));	 	    // l coutour (v 2011)
//
//		imagePoints.push_back(cv::Point(shape.part(0).x(), shape.part(0).y()));	 	// r coutour (v 1138
//		imagePoints.push_back(cv::Point(shape.part(1).x(), shape.part(1).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(2).x(), shape.part(2).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(3).x(), shape.part(3).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(4).x(), shape.part(4).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(5).x(), shape.part(5).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(6).x(), shape.part(6).y()));	 	// r coutour (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(7).x(), shape.part(7).y()));	 	// r coutour (v 1138)
//
//		imagePoints.push_back(cv::Point(shape.part(17).x(), shape.part(17).y()));	 	// r eyebrow(v 1138)
//		imagePoints.push_back(cv::Point(shape.part(18).x(), shape.part(18).y()));	 	// r eyebrow(v 1138)
//		imagePoints.push_back(cv::Point(shape.part(19).x(), shape.part(19).y()));	 	// r eyebrow(v 1138)
//		imagePoints.push_back(cv::Point(shape.part(20).x(), shape.part(20).y()));	 	// r eyebrow(v 1138)
//		imagePoints.push_back(cv::Point(shape.part(21).x(), shape.part(21).y()));	 	// r eyebrow(v 1138)
//
//		imagePoints.push_back(cv::Point(shape.part(22).x(), shape.part(22).y()));	 	// l eyebrow (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(23).x(), shape.part(23).y()));	 	// l eyebrow (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(24).x(), shape.part(24).y()));	 	// l eyebrow (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(25).x(), shape.part(25).y()));	 	// l eyebrow (v 1138)
//		imagePoints.push_back(cv::Point(shape.part(26).x(), shape.part(26).y()));	 	// l eyebrow (v 1138)
//
//
//		cv::Mat ip = cv::Mat(imagePoints);
//		imgPoint = imagePoints;
//		/*sprintf(buf, "%sAngelina_Jolie/Angelina_Jolie_%04d.jpg", workingDir, counter);*/
//		cv::Mat ip_temp = cv::Mat(tempPt);
//		/*cv::Mat img = cv::imread(buf);*/
//
//		//imgTex.set(dst); //TODO: what if different size??
//
//		// paint 2D feature points
//		//for (unsigned int i = 0; i < imagePoints.size(); i++) circle(dst, imagePoints[i], 1, cv::Scalar(255, 0, 255), CV_FILLED);
//		for (unsigned int i = 0; i < shape.num_parts(); i++) {
//			circle(dst, cv::Point(shape.part(i).x(), shape.part(i).y()), 1, cv::Scalar(255, 0, 255), CV_FILLED);
//		}
//		loadWithPoints(ip_temp, dst, imagePoints);
//		cv::imshow("micky", dst);
//		//imgWithDrawing.set(dst);
//	}
//	counter = (counter + 1);
//}
//
//void loadWithPoints(cv::Mat& ip, cv::Mat& img, vector<cv::Point2f> imagePoints) {
//	int max_d = MAX(img.rows, img.cols);
//	camMatrix = (cv::Mat_<double>(3, 3) <<
//		max_d, 0, img.cols / 2.0,
//		0, max_d, img.rows / 2.0,
//		0, 0, 1.0);
//	//cout << "using cam matrix " << endl << camMatrix << endl;
//
//	double _dc[] = { 0,0,0,0 };
//
//	//input 3D point, 2D landmark to get the head pose estimation with:
//	//1. Head position (tvec)
//	//2. Head rotation (face left, face right, face top , face bottom) (rvec)
//	solvePnP(op_min, ip, camMatrix, cv::Mat(1, 4, CV_64FC1, _dc), rvec, tvec, false, CV_EPNP);
//	cv::Mat rotM(3, 3, CV_64FC1, rot);
//	// convert the vector to matrix
//	Rodrigues(rvec, rotM);
//	double* _r = rotM.ptr<double>();
//	//printf("rotation mat: \n %.3f %.3f %.3f\n%.3f %.3f %.3f\n%.3f %.3f %.3f\n",
//	//	_r[0], _r[1], _r[2], _r[3], _r[4], _r[5], _r[6], _r[7], _r[8]);
//
//	//printf("trans vec: \n %.3f %.3f %.3f\n", tv[0], tv[1], tv[2]);
//
//	double _pm[12] = { _r[0],_r[1],_r[2],tv[0],
//		_r[3],_r[4],_r[5],tv[1],
//		_r[6],_r[7],_r[8],tv[2] };
//
//	cv::Matx34d P(_pm);
//	cv::Mat KP = camMatrix * cv::Mat(P);
//	////cout << "KP " << endl << KP << endl;
//
//	//reproject object points - check validity of found projection matrix
//	//op - vertex point in 3D mesh, ip - 2D landmark of the images
//	//for (int i = 0; i < op.rows; i++) {
//	//	
//	//	cv::Mat_<double> X = (cv::Mat_<double>(4, 1) << op.at<float>(i, 0), op.at<float>(i, 1), op.at<float>(i, 2), 1.0);
//	//	//cout << "object point " << X << endl;
//	//	cv::Mat_<double> opt_p = KP * X;
//	//	//cout << "object point " << opt_p << endl;
//	//	double TwoDopX = opt_p(0) / opt_p(2);
//	//	double TwoDopY = opt_p(1) / opt_p(2);
//	//	/*double imageX = imagePoints[i].x;
//	//	double imageY = imagePoints[i].y;*/
//	//	//cout << " opX :" << TwoDopX << " opY :" << TwoDopY << endl;
//	//	/*	double eimX = findEimX(opt_p(0), imageX, opt_p(2),max_d, img.cols / 2.0);
//	//		double eimY = findEimX(opt_p(1), imageY, opt_p(2), max_d, img.rows / 2.0);
//	//		cout << "eimX :" << eimX << " eimY :" << eimY << endl;
//	//		cout << "ein " << eimX + eimY << endl;
//	//		cout << " opX :" << TwoDopX << " opY :" << TwoDopY << endl;
//	//		cout << " imgX :" << imageX << " imgY :" << imageY << endl;*/
//	//		//double diffX = imageX - TwoDopX;
//	//		//double diffY = imageY - TwoDopX;
//	//		/*double dist = sqrt(pow((imageX - TwoDopX), 2) + pow((imageY - TwoDopY), 2));
//	//		cout << "distance :" << dist << " X:" << TwoDopX << " Y:" << TwoDopY  << endl;*/
//	//	//cv::Point2f opt_p_img(TwoDopX, TwoDopY);
//	//	//cout << "object point reproj " << opt_p_img << endl;
//	//	//circle(img, opt_p_img, 4, cv::Scalar(0, 0, 255), 1);
//	//}
//	//transpose for following the head rotation is same as the web cam face rotation otherwise it will move with opposite direction
//	rotM = rotM.t();// transpose to conform with majorness of opengl matrix
//}
//
//
//int main(int argc, char** argv)
//{
//	cout << "read blendshape" << endl;
//	for (int j = 0; j < numOfuser; j++) {
//		vector<GLMmodel*> model;
//		for (int i = 0; i < numOfBlendshape; i++) {
//			char *theString = "FWH/";
//			stringstream ss;
//			ss << theString << j << "/" << "shape_" << i << ".obj";
//			char *cstr = new char[ss.str().length() + 1];
//			strcpy(cstr, ss.str().c_str());
//			model.push_back(glmReadOBJ(cstr));
//			if (j == 0) {
//				char *theString2 = "AverageFace/";
//				stringstream sss;
//				sss << theString2 << "shape_" << i << ".obj";
//				char *cstrr = new char[sss.str().length() + 1];
//				strcpy(cstrr, sss.str().c_str());
//				empty_face.push_back(glmReadOBJ(cstrr));
//			}
//		}
//		rank3Tensor.push_back(model);
//	}
//	vector<cv::Point3d> modelPoint;
//	vector<int> vertexPoint;
//	vertexPoint.push_back(8091);
//	vertexPoint.push_back(2538);
//	vertexPoint.push_back(1263);
//	vertexPoint.push_back(6601);	 // l mouth (v 6601) or line 6605
//	vertexPoint.push_back(728);		 // r mouth (v 728) or line 732 
//	vertexPoint.push_back(9374);	 // l contour (v 8699)
//	vertexPoint.push_back(3917);	 // r contour (v 3484)
//	for (int i = 0; i < vertexPoint.size(); i++) {
//		double x = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i]]);
//		double y = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i] + 1]);
//		double z = roundFloat(empty_face[0]->vertices[3 * vertexPoint[i] + 2]);
//		modelPoint.push_back(cv::Point3d(
//			x,
//			y,
//			z
//		));
//	}
//	dlib::deserialize("shape_predictor_68_face_landmarks.dat") >> sp;
//	op_min = cv::Mat(modelPoint);
//	cout << "finished blendshape :" << endl;
//
//
//	detector = dlib::get_frontal_face_detector();
//	//op = cv::Mat(modelPoints);
//	//cv::Scalar m = mean(cv::Mat(modelPoints));
//	//op = op - m;
//
//	//assert(norm(mean(op)) < 1e-05); //make sure is centered
//
//	//op = op + Scalar(0,0,25);
//
//	//cout << "model points " << op << endl;
//	rvec = cv::Mat(rv);
//	double _d[9] = { 1,0,0,
//		0,-1,0,
//		0,0,-1 };
//	Rodrigues(cv::Mat(3, 3, CV_64FC1, _d), rvec);
//	tv[0] = 0; tv[1] = 0; tv[2] = 1;
//	tvec = cv::Mat(tv);
//
//	camMatrix = cv::Mat(3, 3, CV_64FC1);
//
//	init_opengl(argc, argv); // get GL context, for loading textures
//
//							 // prepare OpenCV-OpenGL images
//	imgTex = MakeOpenCVGLTexture(cv::Mat());
//	imgWithDrawing = MakeOpenCVGLTexture(cv::Mat());
//
//	//loadNext();
//
//	start_opengl();
//
//	return 0;
//}
//
//
