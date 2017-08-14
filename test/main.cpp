#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "face.h"
#define M_PI 3.14159265358979323846
#define FILE_NUM 11


using namespace boost::math;
//#define TRANSFORM

GLfloat angle_y = 0.0f;
GLfloat angle_x = 0.0f;
GLfloat xDirection = 2.0f;
GLfloat yDirection = 2.0f;
GLfloat zDirection = 6.0f;


//wrong assumption
Position Face::default_left_ear = Position(-95, -5, -30);
Position Face::default_right_ear = Position(95, -5, -30);
Position Face::default_left_eye = Position(-20, 20, 80);
Position Face::default_right_eye = Position(20, 20, 80);
Position Face::default_nose = Position(0, -20, 100);


double threshold1 = 0.06;
double eps = 8;
int minPts = 5;

int flag[FILE_NUM];
void display_mesh();

float * vb[FILE_NUM];

float * fvb[FILE_NUM];


int * ib[FILE_NUM];
int nv[FILE_NUM], nf[FILE_NUM], fnv[FILE_NUM];

MatrixXd f_centroid[FILE_NUM];
MatrixXd f_nor[FILE_NUM];

Matrix<std::complex<double>, 3, Dynamic> ff[FILE_NUM];
Face face_model[FILE_NUM];

Position face_vertex_mean[FILE_NUM];
Position face_vertex_var[FILE_NUM];

Position face_nor_mean[FILE_NUM];
Position face_nor_var[FILE_NUM];

int file_index = 0;
int display_index = 0;
int display_mode = 1;
bool isFilter = false;

//Vertex * v[FILE_NUM];

double max_pos = -9999999;
double min_pos = 9999999;//对点云进行在[0,1]之间的归一化

void cal_spheric_harmonic(int face_index)
{
	face_model[face_index].feature_vec.clear();

	char file_path[50];
	FILE * fp;
	sprintf_s(file_path, "..\\2d_feature\\%d\\spheric_harmonic_feature.txt", face_index);
	fopen_s(&fp, file_path, "r");
	if (fp)
	{
		printf("spheric_harmonic has already done\n");
		for (int l = 0; l < 4; l++)
		{
				for (int i = 0; i < 3; i++)
				{
					double feature;
					fscanf_s(fp, "%lf", &feature);
					face_model[face_index].read_feature(feature);
				}
		}
	}
	else
	{
		fopen_s(&fp, file_path, "w");
		for (int l = 0; l < 4; l++)
		{
			Matrix<std::complex<double>, 3, 1> fl(std::complex<double>::complex(0, 0), std::complex<double>::complex(0, 0), std::complex<double>::complex(0, 0));
			std::cout << fl;
			for (int m = -l; m <= l; m++)
			{
				Matrix<std::complex<double>, 3, 1> Alm(std::complex<double>::complex(0, 0), std::complex<double>::complex(0, 0), std::complex<double>::complex(0, 0));
				//std::complex<Vector3d> Alm(Vector3d::Zero(), Vector3d::Zero());
				std::complex<double> Ylm_mul_area(0, 0);
				for (int i = 0; i < nf[face_index]; ++i)
				{
					double theta, phi, r;
					double area;
					Vector3d centroid = f_centroid[face_index].col(i);
					Matrix<std::complex<double>, 3, 1> nor_i = f_nor[face_index].col(i);

					r = centroid.norm();
					theta = acos(centroid.z() / r);
					if (fabs(centroid.x() - 0) < 1e-15)
					{
						if (fabs(centroid.y() - 0) < 1e-15)
						{
							phi = 0;
						}
						else
						{
							if (centroid.y() < 0)
								phi = M_PI;
							else
								phi = 0;
						}
					}
					else
					{
						double temp_phi = atan(centroid.y() / centroid.x());
						if (centroid.x() > 0 && centroid.y() > 0)
							phi = temp_phi;
						else if (centroid.x() < 0 && centroid.y() < 0)
							phi = temp_phi + M_PI;
						else if (centroid.x() < 0 && centroid.y() < 0)
							phi = temp_phi + M_PI;
						else if (centroid.x() > 0 && centroid.y() < 0)
							phi = temp_phi + 2 * M_PI;
					}
					if (i > 0)
						area = face_model[face_index].acc_face_area[i] - face_model[face_index].acc_face_area[i - 1];
					else
						area = face_model[face_index].acc_face_area[i];
					std::complex<double> Ylm_star = spherical_harmonic(l, -m, theta, phi);
					std::complex<double> Ylm = spherical_harmonic(l, m, theta, phi);

					Ylm_mul_area += area * Ylm;
					if (m % 2)//pow(-1, m)
						Ylm_star = -Ylm_star;
					//printf("Ylm%lf %lfi\n", Ylm_star.real(), Ylm_star.imag());

					//printf("area %lf\n", area);
					nor_i *= area;
					//printf("nor_i %lf %lf %lf\n", nor_i[0], nor_i[1], nor_i[2]);
					//Vector3d real = nor_i * Ylm_star.real();
					//Vector3d image = nor_i * Ylm_star.imag();
					//printf("real %lf %lf %lf\n", real[0], real[1], real[2]);
					//printf("image %lf %lf %lf\n", image[0], image[1], image[2]);

					Alm += nor_i * Ylm_star;


				}
				//fl += Alm;
				//Vector3d real = Alm.real();
				//Vector3d image = Alm.imag();

				fl += Ylm_mul_area * Alm;

				printf("real %lf %lf %lf\n", Alm[0].real(), Alm[1].real(), Alm[2].real());
				printf("image %lf %lf %lf\n", Alm[0].imag(), Alm[1].imag(), Alm[2].imag());

				/*for (int ii = 0; ii < 3; ++ii)
				{
					double feature = sqrt(Alm[ii].real() * Alm[ii].real() + Alm[ii].imag() * Alm[ii].imag());
					face_model[face_index].read_feature(feature);
				}*/
			}

			for (int ii = 0; ii < 3; ++ii)
			{
				double feature = sqrt(fl[ii].real() * fl[ii].real() + fl[ii].imag() * fl[ii].imag());
				fprintf_s(fp, "%lf ", feature);
				face_model[face_index].read_feature(feature);
			}

			/*Vector3d real = fl.real();
			Vector3d image = fl.imag();
			Vector3d feature(sqrt(real.x() * real.x() + image.x() * image.x()),
				sqrt(real.y() * real.y() + image.y() * image.y()), sqrt(real.z() * real.z() + image.z() * image.z()));
			for (int ii = 0; ii < 3; ++ii)
			{
				face_model[face_index].read_feature(feature[ii]);
			}*/
		}
	}
	fclose(fp);
}

void cal_facial_feature()
{
	printf("input a threshold, eps, minPts:\n");
	//scanf_s("%lf%lf%d", &threshold1, &eps, &minPts);
	for (size_t k = 0; k < FILE_NUM; k++)
	{
		printf("starting read the curve point for %d face\n", k);
		face_model[k].feature_point.clear();
		face_model[k].feature_vec.clear();
		char file_path[50];
		FILE * fp;
		sprintf_s(file_path, "..\\2d_feature\\%d\\sharp_point.txt", k);
		fopen_s(&fp, file_path, "r");
		if (fp)
		{
			printf("sharp_point is already calculated\n");
			int num_sharp_point;
			fscanf_s(fp, "%d", &num_sharp_point);
			printf("read num_sharp_point\n");
			for (int i = 0; i < num_sharp_point; ++i)
			{
				double point[3];
				for (int j = 0; j < 3; ++j)
				{
					fscanf_s(fp, "%lf", &point[j]);
				}
				Position p(point[0], point[1], point[2]);
				face_model[k].read_right_point(p);
			}
			printf("read done\n");
		}
		else
		{
			fopen_s(&fp, file_path, "w");
			printf("total size if %d\n", face_model[k].v.size());
			int cnt = 0;
			for (int i = 0; i < face_model[k].v.size(); i++)
			{
				double mean_cur = (face_model[k].c[i].max_curvature + face_model[k].c[i].min_curvature) / 2;
				double s = 0.5 - atan((face_model[k].c[i].max_curvature + face_model[k].c[i].min_curvature)
					/ (face_model[k].c[i].max_curvature - face_model[k].c[i].min_curvature)) / M_PI;
				if (mean_cur > threshold1)
				{
					Position p = Position((double)vb[k][i * 3], (double)vb[k][i * 3 + 1], (double)vb[k][i * 3 + 2]);
					face_model[k].read_right_point(p);
					++cnt;
				}
			}
			fprintf_s(fp, "%d", (int)(face_model[k].feature_point.size()));
			for (int i = 0; i < face_model[k].feature_point.size(); ++i)
			{
				Position p = face_model[k].feature_point[i];
				fprintf_s(fp, " %lf %lf %lf\n", p.getX(), p.getY(), p.getZ());
			}
			printf("end reading, read %d points\nstarting DBScan for %d face\n", cnt, k);
		}
		fclose(fp);
		printf("end reading, read %d points\nstarting DBScan for %d face\n", face_model[k].feature_point.size(), k);
		face_model[k].DBScan(eps, minPts, k);//半径为5， 点数目最少为6
		printf("end scaning\n");
		face_model[k].feature_extract(vb[k], nv[k], 10, k);
		
	}
}

void feature_normalize()
{
	std::vector<double> max, min;
	max.clear(), min.clear();
	for (size_t i = 0; i < face_model[0].feature_vec.size(); ++i)//ith feature
	{
		for (size_t j = 0; j < FILE_NUM; ++j)//jth face
		{
			if (max.size() == i)//1st feature
			{
				max.push_back(face_model[j].feature_vec[i]);
				min.push_back(face_model[j].feature_vec[i]);
			}
			else
			{
				if (max[i] < face_model[j].feature_vec[i])
					max[i] = face_model[j].feature_vec[i];
				if (min[i] > face_model[j].feature_vec[i])
					min[i] = face_model[j].feature_vec[i];
			}
		}
	}
	for (int i = 0; i < FILE_NUM; ++i)
	{
		face_model[i].feature_normalize(max, min);
	}
}

int find_nearst_eurler_dis(int index)
{
	double min_dis = 9999999;
	int next_display_index = index;
	for (int i = 0; i < FILE_NUM; i++)
	{
		if (i == index) continue;
		if (!fvb[i]) continue;
		double dis = Face::euler_dis(face_model[index], face_model[i]);
		if (dis < min_dis)
		{
			next_display_index = i;
			min_dis = dis;
		}
	}
	printf("the nearest face  of %d face is %d face, the eurler dis is %lf\n", index, next_display_index, min_dis);
	return next_display_index;
}

int find_nearst_L1_dis(int index)
{
	double min_dis = 9999999;
	int next_display_index;
	for (int i = 0; i < FILE_NUM; i++)
	{
		if (i == index) continue;
		double dis = Face::L1_dis(face_model[index], face_model[i]);
		if (dis < min_dis)
		{
			next_display_index = i;
			min_dis = dis;
		}
	}
	printf("the nearest face  of %d face is %d face, the L1 dis is %lf\n", index, next_display_index, min_dis);
	return next_display_index;
}

int find_nearst_cos_dis(int index)
{
	double max_dis = 0;//余弦距离越大越接近
	int next_display_index;
	for (int i = 0; i < FILE_NUM; i++)
	{
		if (i == index) continue;
		double cos_dis = Face::cos_dis(face_model[index], face_model[i]);
		if (cos_dis > max_dis)
		{
			next_display_index = i;
			max_dis = cos_dis;
		}
	}
	printf("the nearest face  of %d face is %d face, the cos dis is %lf\n", index, next_display_index, max_dis);
	return next_display_index;
}

int find_nearst_ICP(int index, int nor)
{
	double min_dis = 99999;//余弦距离越大越接近
	int next_display_index;
	for (int i = 0; i < FILE_NUM; i++)
	{
		if (i == index) continue;
		double ICP_dis = face_model[i].ICP_rotate(face_model[index], nor);
		if (ICP_dis < min_dis)
		{
			next_display_index = i;
			min_dis = ICP_dis;
		}
	}
	printf("the nearest face  of %d face is %d face, the ICP dis is %lf\n", index, next_display_index, min_dis);
	return next_display_index;
}

void KeyBoards(unsigned char key, int x, int y)
{
	GLfloat pre_xDirection = xDirection;
	GLfloat pre_yDirection = yDirection;
	printf("\nenter keyBoard %c %lf, %lf, %lf, %lf, %lf\n", key, angle_y, angle_x, xDirection, yDirection, zDirection);
	if ('0' <= key && key <= '0' + FILE_NUM)
	{
		display_index = key - '0';
		face_model[display_index].print_feature();
	}
	switch (key)
	{
	case '+':
		display_index++;
		display_index %= FILE_NUM;
		glutPostRedisplay();
		break;
	case '-':
		display_index--;
		if (display_index < 0) display_index = FILE_NUM - 1;
		glutPostRedisplay();
		break;
	case 'z':
		cal_facial_feature();
		break;
	case 'x':
		for (int i = 0; i < FILE_NUM; i++)
		{
			cal_spheric_harmonic(i);
		}
		break;
	case 'c':
		for (int i = 0; i < FILE_NUM; i++)
		{
			if (fvb[i])
			{
				face_model[i].feature_vec.clear();
				face_model[i].cal_2d_feature_2_3d_feature(vb[i], nv[i], i);
			}
		}
		glutPostRedisplay();
		break;
	case 's':
		angle_x += 0.5;
		glutPostRedisplay();
		break;
	case 'e':
		display_mode = 0;
		glutPostRedisplay();
		break;
	case 'w':
		angle_x -= 0.5;
		glutPostRedisplay();
		break;
	case 'q':
		display_mode = 1;
		glutPostRedisplay();
		break;
	case 'f':		
		display_index = find_nearst_eurler_dis(display_index);
		glutPostRedisplay();
		break;
	case 'g':
		display_index = find_nearst_cos_dis(display_index);
		glutPostRedisplay();
		break;
	case 'h':
		display_index = find_nearst_L1_dis(display_index);
		glutPostRedisplay();
		break;
	case 'j':
		display_index = find_nearst_ICP(display_index, 0);
		glutPostRedisplay();
		break;
	case 'k':
		flag[display_index] = 0;
		face_model[display_index].feature_point.clear();
		glutPostRedisplay();
		break;
	case '0':
		glutPostRedisplay();
		break;
	case '4':
		glutPostRedisplay();
		break;
	case '5':
		glutPostRedisplay();
		break;
	case '6':
		glutPostRedisplay();
		break;
	case '7':
		glutPostRedisplay();
		break;
	case '8':
		glutPostRedisplay();
		break;
	case '9':
		glutPostRedisplay();
		break;
	case '1':
		glutPostRedisplay();
		break;
	case '2':
		glutPostRedisplay();
		break;
	case '3':
		glutPostRedisplay();
		break;
	case 27:
		exit(0);
		break;
	case 'd':
		angle_y += 5;
		glutPostRedisplay();
		break;
	/*case 'w':
		xDirection -= 0.05;
		glutPostRedisplay();
		break;*/
	case 'a':
		angle_y -= 5;
		glutPostRedisplay();
		break;
	case 'r':
		isFilter = !isFilter;
		glutPostRedisplay();
		break;
	case 't':
		zDirection += 0.05;
		glutPostRedisplay();
		break;
	case 'y':
		zDirection -= 0.05;
		glutPostRedisplay();
		break;
	case 'v':
		face_model[display_index].rotate_by_2d_feature(vb[display_index], nv[display_index], fvb[display_index]);
		glutPostRedisplay();
		break;
	}
}

void display_mesh(void)
{
	printf("display_index = %d\n", display_index);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  //清空颜色缓冲区  
	glColor3f(0.0, 0.0, 1.0);   //重置颜色  
	glLoadIdentity();   //清空矩阵
	gluLookAt(0, 0, zDirection, 0, 0, 0, 0, 1, 0);
	glTranslatef(0, 0, 0); //将场景中的物体沿z轴负方向移动5个单位长  
							 //glRotatef(40, 0, 1, 0);
							 //gluLookAt(0,0,5,0,0,0,0,2,0); //视点变换  
							 //模型变换  
							 //glutWireCube(1.2); //绘制实心立方体和线框立方体
	glRotatef(angle_y, 0, 1, 0);
	glRotatef(angle_x, 1, 0, 0);
	int cnt = 0;
	if (!display_mode)
	{
		

		glColor3b(0, 0, 1);
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1000);
		glEnd();

		glColor3b(0, 1, 0);
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1000, 0);
		glEnd();

		glColor3b(1, 0, 0);
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(1000, 0, 0);
		glEnd();

		glScalef(0.02, 0.02, 0.02);
		glPointSize(1.0f);
		glColor3b(1, 0, 1);
		glBegin(GL_POINTS);
		if (isFilter)
		{
			for (int i = 0; i < face_model[display_index].feature_point.size(); i++)
			{
				Position p = face_model[display_index].feature_point[i];
				glVertex3f(p.getX(), p.getY(), p.getZ());
			}
			glEnd();
			glPointSize(3.0f);
			int cnt = 0;
			for (int i = 0; i < face_model[display_index].feature_point.size(); i++)
			{
				if (face_model[display_index].class_no[i])
				{
					Position p = face_model[display_index].feature_point[i];
					glColor3f((face_model[display_index].cluster_set[i] % 255) / 255.0, 0.0, 1.0);
					glBegin(GL_POINTS);
					glVertex3f(p.getX(), p.getY(), p.getZ());
					glEnd();
					cnt++;
				}
			}
			cnt = 0;
			glPointSize(10.0f);
			glColor3f(0.0, 0.0, 0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < face_model[display_index].class_mean.size(); i++)
			{
				Position p = face_model[display_index].class_mean[i];
				if (p.getZ() > -40)
				{
					glVertex3f(p.getX(), p.getY(), p.getZ());
					//printf("point at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
					cnt++;
				}
			}
			glEnd();
			if (face_model[display_index].left_eye.getPos() != Position())
			{
				glPointSize(10.0f);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POINTS);
				Position p = face_model[display_index].left_eye.getPos();
				glVertex3f(p.getX(), p.getY(), p.getZ());
				printf("left_eye at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
				glEnd();
			}
			if (face_model[display_index].right_eye.getPos() != Position())
			{
				glPointSize(10.0f);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POINTS);
				Position p = face_model[display_index].right_eye.getPos();
				glVertex3f(p.getX(), p.getY(), p.getZ());
				printf("right_eye at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
				glEnd();
			}
			if (face_model[display_index].right_ear.getPos() != Position())
			{
				glPointSize(10.0f);
				glColor3f(0.0, 1.0, 0.0);
				glBegin(GL_POINTS);
				Position p = face_model[display_index].right_ear.getPos();
				glVertex3f(p.getX(), p.getY(), p.getZ());
				printf("right_ear at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
				glEnd();
			}
			if (face_model[display_index].left_ear.getPos() != Position())
			{
				glPointSize(10.0f);
				glColor3f(0.0, 1.0, 0.0);
				glBegin(GL_POINTS);
				Position p = face_model[display_index].left_ear.getPos();
				glVertex3f(p.getX(), p.getY(), p.getZ());
				printf("left_ear at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
				glEnd();
			}
			if (face_model[display_index].nose.getPos() != Position())
			{
				glPointSize(10.0f);
				glColor3f(0.0, 0.0, 1.0);
				glBegin(GL_POINTS);
				Position p = face_model[display_index].nose.getPos();
				glVertex3f(p.getX(), p.getY(), p.getZ());
				printf("nose at %lf, %lf, %lf\n", p.getX(), p.getY(), p.getZ());
				glEnd();
			}
			//printf("画了有%d重心\n", cnt);
		}
		else
		{
			glPointSize(1.0f);
			glBegin(GL_POINTS);
			for (int i = 0; i < nv[display_index]; ++i)
			{
				glVertex3f(vb[display_index][3 * i], vb[display_index][3 * i + 1], vb[display_index][3 * i + 2]);
			}
			glEnd();
		}
		glPointSize(9.0f);
		glColor3f(0.92, 0.89, 0.41);
		
		for (int i = 0; i < fnv[display_index]; i++)
		{
			glBegin(GL_POINTS);
			glVertex3f(fvb[display_index][3 * i], fvb[display_index][3 * i + 1], fvb[display_index][3 * i + 2]);
			glEnd();
			if (!flag[display_index])
			{
				printf("%d : %f, %f, %f\n", i, fvb[display_index][3 * i], fvb[display_index][3 * i + 1], fvb[display_index][3 * i + 2]);
				//if(i >= 35 && i < 46)
					//system("PAUSE");
				if (i == fnv[display_index] - 1)
					flag[display_index] = 1;
			}
				
			glFlush();
		}
		
	}
	else
	{
		glLineWidth(0.1f);
		glScalef(0.02, 0.02, 0.02);
		for (int i = 0; i < nf[display_index]; i++)
		{
			//printf("i = %d\n", i);
			int index = i * 3;
			float triangle[3][3];
			for (int j = 0; j < 3; j++)
			{
				int point_index = ib[display_index][index + j];
				//printf("point_index = %d\n", point_index);
				if (ib[display_index][index] == ib[display_index][index + 1] || ib[display_index][index + 0] == ib[display_index][index + 2] || ib[display_index][index + 1] == ib[display_index][index + 2])
				{
					//printf("??? at %d\n", i);
				}
				for (int k = 0; k < 3; k++)
				{
					triangle[j][k] = vb[display_index][point_index * 3 + k];
					//printf("k = %d, pos = %f\n", k, vb[point_index + k]);
				}
				//printf("the %d triangle %d vertex at (%f, %f, %f)\n", i, j, triangle[j][0], triangle[j][1], triangle[j][0]);
				//system("PAUSE");
			}

			glBegin(GL_LINE_LOOP);
			glVertex3f(triangle[0][0], triangle[0][1], triangle[0][2]);
			glVertex3f(triangle[1][0], triangle[1][1], triangle[1][2]);
			glVertex3f(triangle[2][0], triangle[2][1], triangle[2][2]);
			glEnd();


			/*glBegin(GL_LINES);
			glVertex3f(triangle[0][0], triangle[0][1], triangle[0][2]);
			glVertex3f(triangle[1][0], triangle[1][1], triangle[1][2]);
			glEnd();

			glBegin(GL_LINES);
			glVertex3f(triangle[0][0], triangle[0][1], triangle[0][2]);
			glVertex3f(triangle[2][0], triangle[2][1], triangle[2][2]);
			glEnd();

			glBegin(GL_LINES);
			glVertex3f(triangle[2][0], triangle[2][1], triangle[2][2]);
			glVertex3f(triangle[1][0], triangle[1][1], triangle[1][2]);
			glEnd();*/
		}
	}
	glFlush();   //刷新窗口以显示当前绘制图形  
}

void init(void)
{
	glEnable(GL_DEPTH);
	glClearColor(1.0, 1.0, 1.0, 0);
	glShadeModel(GL_FLAT); //选择平面明暗模式或光滑明暗模式  */
	
	/*GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };  //镜面反射参数
	GLfloat mat_shininess[] = { 50.0 };               //高光指数
	GLfloat light_position[] = { 0.0, 0.0, 6.0, 0.0 };
	GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };   //灯位置(1,1,1), 最后1-开关
	GLfloat Light_Model_Ambient[] = { 0.2, 0.2, 0.2, 1.0 }; //环境光参数

	glClearColor(0.0, 0.0, 0.0, 0.0);  //背景色
	glShadeModel(GL_SMOOTH);           //多变性填充模式

									   //材质属性
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	//灯光设置
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);   //散射光属性
	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);  //镜面反射光
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Light_Model_Ambient);  //环境光参数

	glEnable(GL_LIGHTING);   //开关:使用光
	glEnable(GL_LIGHT0);     //打开0#灯
	glEnable(GL_DEPTH_TEST); //打开深度测试*/
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);   //设置机口  
	glMatrixMode(GL_PROJECTION);  //指定哪一个矩阵是当前矩阵  
	glLoadIdentity();
	gluPerspective(60, 1, 0, 20);   //创建透视投影矩阵(fovy,aspect,zNear,zFar);  
									  //glFrustum(-1,1,-1,1,1.5,20.0);  //用透视矩阵乘以当前矩阵(left,Right,bottom,top,near,far);  
	glMatrixMode(GL_MODELVIEW);
}

void read_a_face()
{
	max_pos = -9999999;
	min_pos = 9999999;
	FILE * fp;
	char file_path[50];
	sprintf_s(file_path, "..\\2d_feature\\%d\\mesh.bin", file_index);
	fopen_s(&fp, file_path, "rb");


	FILE * fp1;
	char file_path1[50];
	sprintf_s(file_path1, "..\\2d_feature\\%d\\marks.txt", file_index);
	fopen_s(&fp1, file_path1, "r");

	if (fp1)
	{
		fscanf_s(fp1, "%d", &fnv[file_index]);
		fvb[file_index] = (float *)malloc(sizeof(float) * 3 * fnv[file_index]);
		for (int i = 0; i < fnv[file_index]; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				fscanf_s(fp1, "%f", &fvb[file_index][3 * i + j]);
			}
		}
	}
	else
	{
		fnv[file_index] = 0;
		fvb[file_index] = NULL;
		printf("%s :\n", file_path1);
		printf("2Dfeature not exist\n");
	}



	if (!fp)
	{
		printf("file_not_exist\n");
		system("PAUSE");
	}
	if (!feof(fp))
		fread(&nv[file_index], sizeof(int), 1, fp);
	
	if (!feof(fp))
		fread(&nf[file_index], sizeof(int), 1, fp);
	
	vb[file_index] = (float *)malloc(nv[file_index] * 3 * sizeof(float));
	fread(vb[file_index], sizeof(float), nv[file_index] * 3, fp);

	ib[file_index] = (int *)malloc(3 * nf[file_index] * sizeof(int));
	fread(ib[file_index], sizeof(float), nf[file_index] * 3, fp);
	printf("\n\nf0 = (%d, %d, %d)", ib[file_index][0], ib[file_index][1], ib[file_index][2]);
	face_model[file_index] = Face(nv[file_index]);

	f_nor[file_index] = MatrixXd(3, nf[file_index]);
	f_centroid[file_index] = MatrixXd(3, nf[file_index]);

	double sum_x = 0;
	double sum_y = 0;
	double max_y = 0;
	double min_y = 0;
	double sum_z = 0;
	for (int i = 0; i < nv[file_index]; i++)
	{
		double x = vb[file_index][i * 3];
		double y = vb[file_index][i * 3 + 1];
		double z = vb[file_index][i * 3 + 2];
		if (max_y < y)
			max_y = y;
		if (min_y > y)
			min_y = y;
		sum_x += x;
		sum_y += y;
		sum_z += z;
		face_model[file_index].read_vertex(Vertex(Position(x, y, z), Position()));
		for (int j = 0; j < 3; j++)
		{
			if (max_pos < vb[file_index][i * 3 + j])
				max_pos = vb[file_index][i * 3 + j];
			if (min_pos > vb[file_index][i * 3 + j])
				min_pos = vb[file_index][i * 3 + j];
		}
	}
	sum_x /= nv[file_index];
	sum_y = (max_y + min_y) / 2;
	sum_z /= nv[file_index];
	printf("%d model nv = %d, nv[file_index] = %d\n", file_index, face_model[file_index].vertex_num, nv[file_index]);
	MatrixXd m(3, nv[file_index]);

	//0, 1标准化
	for (int i = 0; i < nv[file_index]; i++)
	{
		vb[file_index][i * 3] -= sum_x;
		vb[file_index][i * 3 + 1] -= sum_y;
		vb[file_index][i * 3 + 2] -= sum_z;
		m.block<3, 1>(0, i) = Vector3d(vb[file_index][i * 3], vb[file_index][i * 3 + 1], vb[file_index][i * 3 + 2]);
		//face_model[file_index].v[i].mean(max_pos - min_pos);
		face_model[file_index].v[i].mean(max_pos - min_pos);
		//face_vertex_mean[file_index] += face_model[file_index].v[i].getPos();
		face_vertex_mean[file_index] += face_model[file_index].v[i].getPos();
	}

	//record the 2d feature
	for (int i = 0; i < fnv[file_index]; ++i)
	{
		fvb[file_index][i * 3] -= sum_x;
		fvb[file_index][i * 3 + 1] -= sum_y;
		fvb[file_index][i * 3 + 2] -= sum_z;
		if (i < 15)
		{
			face_model[file_index].face_curve.col(i) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 27)
		{
			face_model[file_index].eyebrow.col(i - 15) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 35)
		{
			face_model[file_index].eye_curve.col(i - 27) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 46)
		{
			face_model[file_index].nose_curve.col(i - 35) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 64)
		{
			face_model[file_index].lips.col(i - 46) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 65)
		{
			face_model[file_index].nose_tip.col(i - 64) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 73)
		{
			face_model[file_index].eyeball.col(i - 65) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else if (i < 75)
		{
			face_model[file_index].eyes.col(i - 73) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
		else
		{
			face_model[file_index].ear_tip.col(i - 75) = Vector3d(fvb[file_index][i * 3], fvb[file_index][i * 3 + 1], fvb[file_index][i * 3 + 2]);
		}
	}


	face_vertex_mean[file_index].mean(nv[file_index]);
	for (int i = 0; i < nf[file_index]; i++)
	{
		int point_index[3];
		Matrix3d triangle;
		for (int j = 0; j < 3; j++)
		{
			point_index[j] = ib[file_index][i * 3 + j];
			Position p =  face_model[file_index].v[point_index[j]].getPos();
			triangle.col(j) = Vector3d(p.getX(), p.getY(), p.getZ());
			if (point_index[j] >= face_model[file_index].v.size() || point_index[j] < 0)
			{
				printf("wrong at f%d,  j%d, v(%d,%d,%d)\n", i, j, point_index[0], point_index[1], point_index[2]);
			}
		}
		index3 index = { point_index[0], point_index[1], point_index[2] };
		f_centroid[file_index].col(i) = Vector3d(triangle.row(0).mean(), triangle.row(0).mean(), triangle.row(0).mean());
		face_model[file_index].read_face(index);
		if (point_index[0] == point_index[1] || point_index[0] == point_index[2] || point_index[2] == point_index[1])
		{
			continue;
		}
		Position v1 = face_model[file_index].v[point_index[0]].getPos() - face_model[file_index].v[point_index[1]].getPos();
		Position v2 = face_model[file_index].v[point_index[0]].getPos() - face_model[file_index].v[point_index[2]].getPos();
		Position nor = Position::cal_nor(v1, v2);//计算得到法向量
		for (int j = 0; j < 3; j++)//向三个顶点加上法向量
		{
			double cos_dis = Position::cos(nor, face_vertex_var[file_index] - face_model[file_index].v[point_index[j]].getPos());
			int flag = 0;
			if (cos_dis < 0)
			{
				if (flag) printf("wrong here %d\n", i);
				flag++;
				//永远朝外方向(背离质心)*/
				nor = -nor;
			}
			double cos_angle = Position::cos(face_model[file_index].v[point_index[(j + 1) % 3]].getPos() - face_model[file_index].v[point_index[j]].getPos(), face_model[file_index].v[point_index[(j + 2) % 3]].getPos() - face_model[file_index].v[point_index[j]].getPos());
			double degree = acos(cos_angle);
			f_nor[file_index].col(i) = Vector3d(nor.getX(), nor.getY(), nor.getZ());
			face_model[file_index].v[point_index[j]].addNor(nor, degree);
		}
	}
	for (int i = 0; i < nv[file_index]; i++)
	{
		face_model[file_index].v[i].normalize();
		face_model[file_index].readVertex2Matrix(face_model[file_index].v[i]);
		face_nor_mean[file_index] += face_model[file_index].v[i].getNor();
		
	}
#ifdef TRANSFORM
	printf("start to cal M\n");
	FILE * fo;
	char out_path[50];
	sprintf_s(out_path, "..\\temp_data\\%d_rotate.bin", file_index);
	fopen_s(&fo, out_path, "wb");

	Matrix3d M = (m * m.transpose() / face_model[file_index].acc_face_area.back());
	EigenSolver<Matrix3d> es(M);
	Matrix3d D = es.pseudoEigenvalueMatrix();
	Matrix3d V = es.pseudoEigenvectors();
	Matrix3d RM;
	int idx[3];
	memset(idx, 0, sizeof(int) * 3);
	double minValue, maxValue;
	minValue = maxValue = fabs(D(0, 0));
	for (int i = 1; i < 3; i++)
	{
		if (minValue > fabs(D(i, i)))
		{
			minValue = fabs(D(i, i));
			idx[2] = i;
		}
		if (maxValue < fabs(D(i, i)))
		{
			maxValue = fabs(D(i, i));
			idx[0] = i;
		}
	}
	idx[1] = 3 - (idx[0] + idx[2]);
	for (int i = 0; i < 3; i++)
	{
		RM.block<3, 1>(0, i) = V.col(idx[i]);
	}
	printf("end  cal M\n");
	fwrite(&nv[file_index], sizeof(int), 1, fo);
	fwrite(&nf[file_index], sizeof(int), 1, fo);
	printf("start to cal new point\n");
	for (int i = 0; i < nv[file_index]; i++)
	{
		Vector3d point(vb[file_index][3 * i], vb[file_index][3 * i + 1], vb[file_index][3 * i + 2]);
		point = RM * point;
		float x[3];
		x[0] = point[0];
		x[1] = point[1];
		x[2] = point[2];
		fwrite(&x, sizeof(float), 3, fo);
	}
	printf("end  cal new point\n");
	fwrite(ib[file_index], sizeof(float), nf[file_index] * 3, fo);
	fclose(fo);
#else	

	face_model[file_index].cal_curvature();
	face_model[file_index].cal_varOfVertexNNor();
	//face_model[file_index].add_D2_feature(100);
	printf("v.size %d vertex.cols = %d\n", face_model[file_index].v.size(), face_model[file_index].vertex.cols());
	printf(" %d face\n", file_index);
	face_model[file_index].print_feature();
	fclose(fp);
	if(fp1)
		fclose(fp1);
#endif // TRANSFORM
}

int main(int argc, char** argv) {
	for (; file_index < FILE_NUM; file_index++)
	{
		read_a_face();
		//cal_spheric_harmonic(file_index);
		if (fvb[file_index])
		{
			//face_model[file_index].rotate_by_2d_feature(vb[file_index], nv[file_index], fvb[file_index]);
			//face_model[file_index].cal_2d_feature_2_3d_feature(vb[file_index], nv[file_index]);
		}
	}
		
#ifndef TRANSFORM
	//cal_facial_feature();
	//feature_normalize();

	//display = display_mesh;
	glutInit(&argc, argv); //固定格式  
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);   //缓存模式  
	glutInitWindowSize(800, 800);    //显示框的大小  
	glutInitWindowPosition(400, 400); //确定显示框左上角的位置  
	glutCreateWindow("绘制立方体");
	init();
	glutKeyboardFunc(&KeyBoards);  //注册键盘事件 
	glutDisplayFunc(display_mesh);
	glutReshapeFunc(reshape);
	glutMainLoop(); //进人GLUT事件处理循环  
#endif // !TRANSFORM
	return 0;

}