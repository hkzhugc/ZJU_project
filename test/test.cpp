#include "stdlib.h"
#include"time.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "windows.h"
#include "iostream"
#include"fstream"
using namespace std;
#define PI 3.1415926535  

float thetaX = 0.0, thetaY = 0.0, scaleFactor = 1.0;
static float dx = 0, dy = 0, oldy = -1, oldx = -1;
int width = 300, height = 300;       //设定窗口大小
struct Point
{
	float x;
	float y;
	float z;
};
vector data;
vector ptsmm;
Point ptsCen = {};             //点云中心
void drawBunny()
{
	glPointSize(1.0f);
	glColor3f(0.0, 1.0, 0.0);    //点云是绿色的
	glBegin(GL_POINTS);
	for (int i = 0; i
	{
		glVertex3f(data[i].x,data[i].y,data[i].z);
	}
	glEnd();
		glPointSize(5.0f);
		glColor3f(1.0, 0.0, 0.0);    //原点是红色的
		glBegin(GL_POINTS);
		glVertex3f(0, 0, 0);
		glEnd();
}
void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);     //设置填屏颜色
	glShadeModel(GL_FLAT);
}
void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	if (thetaY<0)
	{
		thetaY = thetaY + 360;
	}
	if (thetaY>360)
	{
		thetaY = thetaY - 360;
	}
	if (thetaX<0)
	{
		thetaX = thetaX + 360;
	}
	if (thetaX>360)
	{
		thetaX = thetaX - 360;
	}
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 2, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);      //设置照相机位置
														   //glTranslatef(ptsCen.x,ptsCen.y,ptsCen.z);
	glRotatef(thetaX, 1, 0, 0);
	glRotatef(thetaY, 0, 1, 0);
	glScalef(scaleFactor, scaleFactor, scaleFactor);
	glTranslatef(-ptsCen.x, -ptsCen.y, -ptsCen.z);
	drawBunny();
	glutSwapBuffers();
}
void reshape(int width, int height)
{
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);         //投影变换，确定显示空间的大小
	glLoadIdentity();
	glOrtho(-1.5, 1.5, -1.5, 1.5, -5, 5);
}
void keyBoard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'A':              //放大和缩小
	case 'a':
		scaleFactor *= 0.9;
		glutPostRedisplay();
		break;
	case 'D':
	case 'd':
		scaleFactor *= 1.1;
		glutPostRedisplay();
		break;
	case 'R':                     //恢复原状
	case 'r':
		thetaX = 0;
		thetaY = 0;
		scaleFactor = 1.0;
		glutPostRedisplay();
		break;
	case 'Q':                    //退出程序
	case 'q':
		exit(0);
		break;
	default:
		break;
	}
}
void myMouse(int button, int state, int x, int y)        //处理鼠标点击 
{
	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)                  //左键按下时,记录初始坐标 
		oldx = x, oldy = y;
	if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)             //点击右键时还原
	{
		thetaX = 0; thetaY = 0; scaleFactor = 1;
		glutPostRedisplay();
	}
	if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
	{
	}
}
void onMouseMove(int x, int y)     //处理鼠标拖动 
{
	dx += x - oldx;
	dy += y - oldy;
	thetaX = dy / width * 90;
	thetaY = dx / width * 90;
	oldx = x, oldy = y;               //用新值覆盖旧值
	glutPostRedisplay();
}


//////////////////////////////////

void calminmax(vector & ptsmm, vector pts);

int _tmain(int argc, char ** argv)
{
	ifstream f("H:\\Dragon.txt");
	if (!f)
	{
		cout << "cannot read file!!" << endl;
		getchar();
	}
	Point temp = {};
	while (!f.eof())
	{
		f >> temp.x;
		f >> temp.y;
		f >> temp.z;
		data.push_back(temp);
	}
	calminmax(ptsmm, data);
	float  width = ptsmm[0] - ptsmm[1];
	float height = ptsmm[2] - ptsmm[3];
	float depth = ptsmm[4] - ptsmm[5];
	int dataLength = data.size();
	ptsCen.x = (ptsmm[0] + ptsmm[1]) / 2; ptsCen.x = (ptsmm[2] + ptsmm[3]) / 2; ptsCen.x = (ptsmm[4] + ptsmm[5]) / 2;
	cout << "Points Number: " << dataLength << endl << "depth: " << depth << endl << "height: " << height << endl << "width: " << width << endl;
	cout << "Space scope: " << endl << "x: " << ptsmm[1] << "  " << ptsmm[0] << endl << "y: " << ptsmm[3] << "  " << ptsmm[2] << endl << "z: " << ptsmm[5] << "  " << ptsmm[4] << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(600, 600);
	glutCreateWindow("Bunny Display");
	glutReshapeFunc(reshape);
	glutDisplayFunc(&display);
	glutKeyboardFunc(keyBoard);
	init();
	glutMouseFunc(myMouse);
	glutMotionFunc(onMouseMove);
	glutMainLoop();
	return 0;
}

void calminmax(vector & ptsmm, vector data)     //分别计算x,y,z的最大值最小值
{
	int i, j;
	float tmax, tmin;
	float ** a;
	a = new float *[3];
	for (i = 0; i<3; i++)
	{
		a[i] = new float[data.size()];
	}
	for (i = 0; i
	{
		a[0][i] = data[i].x;
		a[1][i] = data[i].y;
		a[2][i] = data[i].z;
	}
	for (i = 0; i<3; i++)
	{
		tmax = a[i][0];
			tmin = a[i][0];
			for (j = 0; j
			{
				if (a[i][j]>tmax)
				{
					tmax = a[i][j];
				}
		if (a[i][j]
		{
			tmin = a[i][j];
		}

			}
		ptsmm.push_back(tmax);
		ptsmm.push_back(tmin);
	}

	for (int i = 0; i<3; i++)    //由里向外释放内存
	{
		delete[]a[i];
		a[i] = NULL;
	}
	delete[]a;
	a = NULL;
}