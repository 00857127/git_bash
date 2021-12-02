#include<stdio.h>
#include<math.h>
#include<GL/glut.h>
typedef struct point {
	long double x;
	double y;
};
point grid_points[60];
int count = 0, countC = 0;;
long double C[500];

void reshape(int w, int h) {
	glViewport(0, 0, w, h);	//調整視角
}

long double f(long double x, long double y) {
	return (((x - 3) * (x - 3) / 64) + ((y - 4) * (y - 4) / 64) - 1);
}

long double g(long double x, long double y) {
	return ((x - 3) * (y - 4) - 1);
}

long double fx(long double x, long double y) {
	return ((x - 3) / 32);
}

long double fy(long double x, long double y) {
	return ((y - 4) / 32);
}

long double gx(long double x, long double y) {
	return (y - 4);
}

long double gy(long double x, long double y) {
	return (x - 3);
}

void eigenvalues(long double fx, long double fy, long double gx, long double gy) {
	long double b, c, D, lamda_max, lamda_min;
	b = -(gy + fx);
	c = (fx * gy - fy * gx);
	D = sqrt(b * b - 4 * c);
	lamda_max = (-b + D) / 2;
	lamda_min = (-b - D) / 2;
	if (lamda_min > lamda_max) {
		long double temp = lamda_min;
		lamda_min = lamda_max;
		lamda_max = temp;
	}
	C[countC] = fabs(lamda_max) / fabs(lamda_min);
	countC++;
}

void newton2D(long double x, long double y, int type) {
	int i = 0, lock;
	long double xold, yold, xnew, ynew, err;
	long double fold, gold, fdx, fdy, gdx, gdy, delta;
	if (type == 3) {
		xnew = x + sqrt(fabs(x));
		ynew = y + 0.01;
	}
	else {
		xnew = x;
		ynew = y;
	}
	lock = 1;
	err = 1;
	if (type==1) {
		printf("  i     xn          yn          error\n");
		printf("------------------------------------------------------------\n");
		printf("%3d  %.8lf  %.8lf  ---\n", i, xnew, ynew);
	}
	while (err >= 0.000001) {
		//設定初值
		xold = xnew;
		yold = ynew;
		//計算各函式
		fold = f(xold, yold);
		gold = g(xold, yold);
		fdx = fx(xold, yold);
		fdy = fy(xold, yold);
		gdx = gx(xold, yold);
		gdy = gy(xold, yold);
		if (type == 0 && lock == 1) {
			eigenvalues(fdx, fdy, gdx, gdy);
			lock = 0;
		}
		delta = fdx * gdy - fdy * gdx;
		//計算根
		xnew = xold - (gdy * fold - fdy * gold) / delta;
		ynew = yold - (-gdx * fold + fdx * gold) / delta;
		i++;	
		//計算2-norm err
		err = sqrt((xnew - xold) * (xnew - xold) + (ynew - yold) * (ynew - yold));
		//output
		if (type == 1) {
			printf("%3d  %.8lf  %.8lf  %.8lf\n", i, xnew, ynew, err);
		}
	}
	
	//判斷是否為grid_points
	if (type == 0) {
		if (isnan(xnew) && isnan(ynew)) {		
			grid_points[count].x = x;
			grid_points[count].y = y;
			count++;
		}
	}
	if (type == 9) {
		if (isnan(xnew) && isnan(ynew)) {
			glColor4f(1,0, 0, 0.1);
			glVertex2f(x / 16, y / 16);
		}
		else {
			glColor4f(0, 1, 1, 0.1);
			glVertex2f(x / 16, y / 16);
		}
	}
	//非grid_points的迭代次數
	if (type == 2) {
		printf("(%lf,%lf) iterations = %d\n", x, y, i);
	}
	//grid_points加入perturbation後的狀態
	if (type == 3) {
		if (isnan(xnew) && isnan(ynew)) {
			printf("(%lf,%lf) is still disconverge\n", x, y);
		}
		else {
			printf("(%lf,%lf) becomes converge\n", x, y);
		}
	}
}

void display_func(void) {
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	//畫g函式
	float x1, y1;
	glBegin(GL_LINE_STRIP);
	glColor3f(1, 1, 1);
	for (x1 = -6.0; x1 < 3; x1 = x1 + 0.05) {
		y1 = (4 * x1 - 11) / (x1 - 3);
		glVertex2f(x1 / 16, y1 / 16);
	}
	glEnd();
	glBegin(GL_LINE_STRIP);
	glColor3f(1, 1, 1);
	for (x1 = 3.05; x1 <= 12.0; x1 = x1 + 0.05) {
		y1 = (4 * x1 - 11) / (x1 - 3);
		glVertex2f(x1 / 16, y1 / 16);
	}
	glEnd();
	//畫f函式	
	glBegin(GL_LINE_STRIP);
	glColor3f(1, 1, 1);
	float R = 8;
	for (int i = 0; i < 3600; i++) {
		glVertex2f((R * cos(i / (float)360) + 3) / 16, (R * sin(i / (float)360) + 4) / 16);
	}
	glEnd();
	//畫點
	glPointSize(3);
	glBegin(GL_POINTS);
	glColor3f(1, 1, 1);
	for (float i = -50; i < 50; i++) {
		glVertex2f(0.0, i / 16);
		glVertex2f(i / 16, 0.0);
	}
	glEnd();
	//畫格線
	glBegin(GL_LINES);
	glColor3f(0.3, 0.3, 0.3);
	glVertex2f(-1.0, 0.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(0.0, 1.0);
	glVertex2f(0.0, -1.0);
	for (int i = -50; i < 50; i++) {
		glVertex2f(float(i) / 16, 1.0);
		glVertex2f(float(i) / 16, -1.0);
		glVertex2f(1.0, float(i) / 16);
		glVertex2f(-1.0, float(i) / 16);
	}
	glEnd();
	//3
	/*glBegin(GL_POLYGON);
	glColor4f(0, 0.5, 0.5, 0.5);
	glVertex2f(-6.0/16, -5.0/16);
	glVertex2f(12.0/16, -5.0/16);	
	glVertex2f(12.0/16, 13.0/16);
	glVertex2f(-6.0/16, 13.0/16);
	glEnd();
	glBegin(GL_LINES);
	for (int p = 0; p < count; p++) {
			glColor3f(1, 0, 0);
			glVertex2f((float)grid_points[p].x / 16, (float)grid_points[p].y / 16);
	}
	glEnd();*/

	glBegin(GL_POINTS);
	for (long double i = -6; i < 12; i = i + 0.5) {
		for (long double j = -5; j < 13; j = j + 0.5) {
			newton2D(i, j, 9);
		}
	}
	glEnd();

	glFlush();
}

int main(int argc, char** argv) {
	//1-(a)
	printf("1-(a): compute initial solution(2.0,-4.0)\n");
	newton2D(2.0, -4.0, 1);
	//1-(b)
	printf("-----------------------------------------------------------------------------\n");
	printf("1-(b): list all grid-points which converge\n");
	for (int i = -6; i <= 12; i++) {
		for (int j = -5; j <= 13; j++) {
			newton2D(i, j, 0);
		}
	}
	for (int i = -6; i <= 12; i++) {
		for (int j = -5; j <= 13; j++) {
			int not_print = 0;
			for (int p = 0; p < count; p++) {
				if ((long double)i == grid_points[p].x && (long double)j == grid_points[p].y) {
					not_print = 1;
					break;
				}
			}
			if (!not_print) {
				printf("(%lf,%lf)\n", (long double)i, (long double)j);
			}
		}
	}
	
	//1-(c)
	printf("-----------------------------------------------------------------------------\n");
	printf("1-(c): iterations of converge points\n");
	for (int i = -6; i <= 12; i++) {
		for (int j = -5; j <= 13; j++) {
			newton2D(i, j, 2);
		}
	}
	//1-(d)
	printf("-----------------------------------------------------------------------------\n");
	printf("1-(d): grid-points which diverge add small perturbarion\n");
	for (int i = 0; i < count; i++) {
		newton2D(grid_points[i].x, grid_points[i].y, 3);
	}
	
	//4-(a)
	printf("-----------------------------------------------------------------------------\n");
	printf("4-(a): 若初始解沒有使delta等於零，則會converge\n");
	//4-(b)
	printf("-----------------------------------------------------------------------------\n");
	printf("4-(b): print all points and eigenvalues\n");
	int k = 0;
	for (int i = -6; i <= 12; i++) {
		for (int j = -5; j <= 13; j++) {
			int exp = 1;
			printf("(%lf,%lf) eigenvalues is %lf ", (long double)i, (long double)j, C[k]);
			for (int p = 0; p < count; p++) {
				if ((long double)i == grid_points[p].x && (long double)j == grid_points[p].y) {
					printf("which diverge\n");
					exp = 0;
				}
			}
			if (exp) {
				printf("which converge\n");
			}
			k++;
		}
	}
	printf("根據以上資料，若 eigenvalues 等於零、無限、NAN則為 disconverge\n");
	//2-(a) and 3
	glutInit(&argc, argv);
	glutInitWindowPosition(500, 250);	//窗口初始位置
	glutInitWindowSize(700, 700);	//窗口初始大小
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);	//設定模式
	glutCreateWindow("2-(a) and 3");
	glutReshapeFunc(reshape);
	glutDisplayFunc(display_func);
	glutMainLoop();
	return 0;
}
