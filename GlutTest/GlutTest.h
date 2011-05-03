
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "imageloader.h"


#if defined(_WIN32) || defined(__CYGWIN__)
# include <windows.h>
# include "GL\glut.h"
#elif defined(__APPLE__)
# include <GLUT/glut.h>
#else 
# include <GL/glut.h>
#endif

void GLUTResize(int w, int h);
//static void drawBrickBox(GLfloat size, GLenum type);
//void glutSolidBrickCube(GLdouble size);
//void drawBuilding();
GLuint loadTexture(Image* image);
//void computeDir(float deltaAngle);
//void computePos(float deltaMoveX, float deltaMoveY);
//void moveForward();
void GLUTRedraw(void);
void GLUTKeyboard(unsigned char key, int xx, int yy);
void GLUTSpecial(int key, int xx, int yy);
void releaseKey(int key, int x, int y);





