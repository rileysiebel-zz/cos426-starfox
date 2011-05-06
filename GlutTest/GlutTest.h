
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

//ship movement and view functions
void computeDir(double deltaAngle);
void computePos(double deltaMoveX, double deltaMoveY);
void moveForward(void);
void updateShip(void);
void peakRight(void);
void peakLeft(void);
void peakUp(void);
void peakDown(void);
void lookStraightLR(void);
void lookStraightUD(void);
void computeRotation(void);


//enemy and projectile functions
void updateEnemies(void);
void updateProjectiles(void);


void GLUTRedraw(void);
void GLUTKeyboard(unsigned char key, int xx, int yy);
void GLUTSpecial(int key, int xx, int yy);
void releaseKey(int key, int x, int y);