//#include <stdio.h>
//#include <math.h>
//#include <vector>
//#include "opencv\cv.h"
//#include "opencv\highgui.h"
//#include "glut.h"
//
//#include "glm.h"
//#define T(x) (glm_model->triangles[(x)])
//#define	G_PI 3.14159265358979323846f
//
//void prepare_lighting();
//void display();
//void keyboard(unsigned char key, int x, int y);
//GLMmodel *glm_model, *glm_model2;
//float theta, phi;
//float change;
//
//GLuint list_id;
//
//void main(int argc, char** argv)
//{
//	theta = G_PI / 2;
//	phi = -G_PI / 2;
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
//	glutInitWindowSize(640, 640);
//	glutCreateWindow("glutTest08");
//	std::vector<cv::Point3f > modelPoints;
//	modelPoints.push_back(cv::Point3f(0.294959, 0.493100, 0.389258));	// l eye (v 314)
//	modelPoints.push_back(cv::Point3f(-0.263214, 0.495387, 0.373535));	// r eye (v 0)
//	modelPoints.push_back(cv::Point3f(0.029220, 0.188548, 0.673340));	//nose (v 1879)
//	modelPoints.push_back(cv::Point3f(0.266336, -0.128219, 0.464839));	// l mouth (v 1502)
//	modelPoints.push_back(cv::Point3f(-0.171941, -0.121142, 0.518103));	// r mouth (v 695)
//	modelPoints.push_back(cv::Point3f(0.857392, 0.459262, -0.547478));	// l ear (v 2011)
//	modelPoints.push_back(cv::Point3f(-0.853804, 0.425016, -0.562550));		// r ear (v 1138)
//	glutDisplayFunc(display);
//	glutKeyboardFunc(keyboard);
//	{
//		glm_model = glmReadOBJ("shape_1.obj");
//		glm_model2 = glmReadOBJ("shape_0_blender_2.obj");
//		int numtri = glm_model->numvertices;
//		//int nunmtri2 = glm_model2->groups->numtriangles;
//		for (int i = 1; i < numtri; i++) {
//			GLMtriangle* triangle = &T(glm_model->groups->triangles[i]);
//			GLMtriangle* triangle2 = &T(glm_model2->groups->triangles[i]);
//			//printf("vindices %f %f \n", glm_model->vertices[3 * i],glm_model2->vertices[3 * i]);
//			//if (glm_model->vertices[3 * i] == glm_model2->vertices[3 * i]) {
//			//	//printf("same x");
//			//	continue;
//			//}
//			glm_model->vertices[3 * i] = (glm_model->vertices[3 * i] + glm_model2->vertices[3 * i]);
//			glm_model->vertices[3 * i + 1] = (glm_model->vertices[3 * i + 1] + glm_model2->vertices[3 * i + 1]);
//			glm_model->vertices[3 * i + 2] = (glm_model->vertices[3 * i + 2] + glm_model2->vertices[3 * i + 2]);
//			/*printf("vindices %f %f %f\n", 3 * triangle->vindices[0], 3 * triangle->vindices[1], 3 * triangle->vindices[2]);
//				printf("vindex :%f %f %f\n", glm_model->vertices[3 * i], glm_model->vertices[3 * i + 1], glm_model->vertices[3 * i + 2]);
//				printf("tindex : %f %f %f\n", &glm_model->texcoords[2 * triangle->tindices[0]], &glm_model->texcoords[2 * triangle->tindices[1]], &glm_model->texcoords[2 * triangle->tindices[1]]);
//				printf("nindex : %f %f %f\n", triangle->nindices[0], triangle->nindices[1], triangle->nindices[2]);
//				printf("findex : %f\n", triangle->findex);
//
//				printf("---------------------\n");*/
//		}
//		//printf("nindex : %f %f %f\n", v1, v2, v3);
//		//for (int i = 0; i < numtri; i++) {
//		//	GLMtriangle* triangle = &T(glm_model->groups->triangles[i]);
//		//	GLMtriangle* triangle2 = &T(glm_model2->groups->triangles[i]);
//		//	glm_model->vertices[3 * triangle->vindices[0]] = glm_model->vertices[3 * triangle->vindices[0]] * v1;
//		//	glm_model->vertices[3 * triangle->vindices[1]] = glm_model->vertices[3 * triangle->vindices[1]] * v2;
//		//	glm_model->vertices[3 * triangle->vindices[2]] = glm_model->vertices[3 * triangle->vindices[2]] * v3;
//
//		//	/*	printf("vindex :%f %f %f\n", glm_model->vertices[3 * triangle->vindices[0]], glm_model->vertices[3 * triangle->vindices[1]], glm_model->vertices[3 * triangle->vindices[2]]);
//		//	printf("tindex : %f %f %f\n", &glm_model->texcoords[2 * triangle->tindices[0]], &glm_model->texcoords[2 * triangle->tindices[1]], &glm_model->texcoords[2 * triangle->tindices[1]]);
//		//	printf("nindex : %f %f %f\n", triangle->nindices[0], triangle->nindices[1], triangle->nindices[2]);
//		//	printf("findex : %f\n", triangle->findex);
//
//		//	printf("---------------------\n");*/
//		//}
//
//	/*	glmUnitize(glm_model);
//		glmScale(glm_model, .1);
//		glmFacetNormals(glm_model);
//		glmVertexNormals(glm_model, 90);
//
//		list_id = glmList(glm_model, GLM_MATERIAL | GLM_SMOOTH);
//
//		glmDelete(glm_model);*/
//	}
//
//	prepare_lighting();
//
//	glutMainLoop();
//}
//
//void keyboard(unsigned char key, int x, int y)
//{
//	switch (key)
//	{
//	case 'w':
//		theta -= .05;
//		prepare_lighting();
//		glutPostRedisplay();
//		break;
//
//	case 's':
//		theta += .05;
//		prepare_lighting();
//		glutPostRedisplay();
//		break;
//
//	case 'a':
//		phi -= .05;
//		prepare_lighting();
//		glutPostRedisplay();
//		break;
//
//	case 'd':
//		phi += .05;
//		prepare_lighting();
//		glutPostRedisplay();
//		break;
//	};
//}
//
//
//void display()
//{
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	gluPerspective(20, 1, 0.1, 10);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//	gluLookAt(
//		0, 0, 1,
//		0, 0, 0,
//		0, 1, 0);
//
//	glEnable(GL_LIGHTING);
//	glEnable(GL_DEPTH_TEST);
//
//	glCallList(list_id);
//
//	glutSwapBuffers();
//}
//
//void prepare_lighting()
//{
//	theta = fmodf(theta, 2 * G_PI);
//	phi = fmodf(phi, 2 * G_PI);
//
//	float light_diffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
//	float mat_diffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
//	float light_position[4] = { sinf(theta) * cosf(phi), cosf(theta), -sinf(theta) * sinf(phi), 0 };
//
//	//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
//	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
//	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//	glEnable(GL_LIGHT0);
//	//glm_model->vertices = glm_model->vertices + glm_model2->vertices;
//	//int numtri = glm_model->groups->numtriangles;
//	//int nunmtri2 = glm_model2->groups->numtriangles;
//	//for (int i = 0; i < numtri; i++) {
//	//	GLMtriangle* triangle = &T(glm_model->groups->triangles[i]);
//	//	GLMtriangle* triangle2 = &T(glm_model2->groups->triangles[i]);
//
//	//	/*glm_model->vertices[3 * triangle->vindices[0]] = (glm_model->vertices[3 * triangle->vindices[0]]);
//	//	glm_model->vertices[3 * triangle->vindices[1]] = (glm_model->vertices[3 * triangle->vindices[1]]);
//	//	glm_model->vertices[3 * triangle->vindices[2]] = (glm_model->vertices[3 * triangle->vindices[2]]);*/
//	//	printf("vindex :%f %f %f\n", glm_model->vertices[3 * triangle->vindices[0]], glm_model->vertices[3 * triangle->vindices[1]], glm_model->vertices[3 * triangle->vindices[2]]);
//	//	printf("tindex : %f %f %f\n", &glm_model->texcoords[2 * triangle->tindices[0]], &glm_model->texcoords[2 * triangle->tindices[1]], &glm_model->texcoords[2 * triangle->tindices[1]]);
//	//	printf("nindex : %f %f %f\n", triangle->nindices[0], triangle->nindices[1], triangle->nindices[2]);
//	//	printf("findex : %f\n", triangle->findex);
//
//	//	printf("---------------------\n");
//	//}
//	//printf("nindex : %f %f %f\n", v1, v2, v3);
//	//for (int i = 0; i < numtri; i++) {
//	//	GLMtriangle* triangle = &T(glm_model->groups->triangles[i]);
//	//	GLMtriangle* triangle2 = &T(glm_model2->groups->triangles[i]);
//	//	glm_model->vertices[3 * triangle->vindices[0]] = glm_model->vertices[3 * triangle->vindices[0]] * v1;
//	//	glm_model->vertices[3 * triangle->vindices[1]] = glm_model->vertices[3 * triangle->vindices[1]] * v2;
//	//	glm_model->vertices[3 * triangle->vindices[2]] = glm_model->vertices[3 * triangle->vindices[2]] * v3;
//
//	//	/*	printf("vindex :%f %f %f\n", glm_model->vertices[3 * triangle->vindices[0]], glm_model->vertices[3 * triangle->vindices[1]], glm_model->vertices[3 * triangle->vindices[2]]);
//	//	printf("tindex : %f %f %f\n", &glm_model->texcoords[2 * triangle->tindices[0]], &glm_model->texcoords[2 * triangle->tindices[1]], &glm_model->texcoords[2 * triangle->tindices[1]]);
//	//	printf("nindex : %f %f %f\n", triangle->nindices[0], triangle->nindices[1], triangle->nindices[2]);
//	//	printf("findex : %f\n", triangle->findex);
//
//	//	printf("---------------------\n");*/
//	//}
//
//	glmUnitize(glm_model);
//	glmScale(glm_model, .1);
//	glmFacetNormals(glm_model);
//	glmVertexNormals(glm_model, 90);
//
//	list_id = glmList(glm_model, GLM_MATERIAL | GLM_SMOOTH);
//
//	//glmDelete(glm_model);
//}