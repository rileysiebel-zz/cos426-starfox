#include "GlutTest.h"

using namespace std;

// texture variables
GLuint texBrick;
GLuint texStone;

// camera direction
float lx=0.0f,ly=0.0f,lz=-1.0f;

// XYZ position of the camera
float x=0.0f, y=1.0f, z=-5.0f;

// the key states. These variables will be zero when no key is being presses
float deltaMoveX = 0, deltaMoveY = 0;

// angles
// angle of rotation for the camera direction
float angle = 0.0f;
float deltaAngle = 0.0f;

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio =  w * 1.0 / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

// code adapted from http://www.opengl.org/resources/libraries/glut/ 
// changed the make box function to include the texure
static void drawBrickBox(GLfloat size, GLenum type)
{
  static GLfloat n[6][3] =
  {
    {-1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0},
    {0.0, 0.0, 1.0},
    {0.0, 0.0, -1.0}
  };
  static GLint faces[6][4] =
  {
    {0, 1, 2, 3},
    {3, 2, 6, 7},
    {7, 6, 5, 4},
    {4, 5, 1, 0},
    {5, 6, 2, 1},
    {7, 4, 0, 3}
  };
  GLfloat v[8][3];
  GLint i;

  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -size / 2;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] = size / 2;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] = -size / 2;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] = size / 2;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] = -size / 2;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] = size / 2;

	glBindTexture(GL_TEXTURE_2D, texBrick);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glColor3f(1.0f, 1.0f, 1.0f);

  for (i = 5; i >= 0; i--) {
    glBegin(type);
    glNormal3fv(&n[i][0]);
	glTexCoord2f(0.0f, 0.0f);
    glVertex3fv(&v[faces[i][0]][0]);
	glTexCoord2f(0.0f, 1.0f);
    glVertex3fv(&v[faces[i][1]][0]);
	glTexCoord2f(1.0f, 1.0f);
    glVertex3fv(&v[faces[i][2]][0]);
	glTexCoord2f(1.0f, 0.0f);
    glVertex3fv(&v[faces[i][3]][0]);
    glEnd();
  }
}

// complementing the above function
void glutSolidBrickCube(GLdouble size)
{
  drawBrickBox(size, GL_QUADS);
}

void drawBuilding() {

	glTranslatef(0.0f ,0.0f, 0.0f);
	glutSolidBrickCube(1.0f);

}

//Makes the image into a texture, and returns the id of the texture
GLuint loadTexture(Image* image) {
	GLuint textureId;
	glGenTextures(1, &textureId); //Make room for our texture
	glBindTexture(GL_TEXTURE_2D, textureId); //Tell OpenGL which texture to edit
	//Map the image to the texture
	glTexImage2D(GL_TEXTURE_2D,                //Always GL_TEXTURE_2D
				 0,                            //0 for now
				 GL_RGB,                       //Format OpenGL uses for image
				 image->width, image->height,  //Width and height
				 0,                            //The border of the image
				 GL_RGB, //GL_RGB, because pixels are stored in RGB format
				 GL_UNSIGNED_BYTE, //GL_UNSIGNED_BYTE, because pixels are stored
				                   //as unsigned numbers
				 image->pixels);               //The actual pixel data
	return textureId; //Returns the id of the texture
}

void computeDir(float deltaAngle) {

	angle += deltaAngle;
	if (angle > .02f)
		angle = .02f;
	if (angle < -.02f)
		angle = -.02f;

	lx = sin(angle);
	lz = -cos(angle);
}

void computePos(float deltaMoveX, float deltaMoveY) {

	x += deltaMoveX;
	y += deltaMoveY;

}

// move forward with a constant speed
void moveForward() {
	z -= 0.01f;
}

void renderScene(void) {

	moveForward();

	if (deltaMoveX || deltaMoveY)
		computePos(deltaMoveX, deltaMoveY);

	if (deltaAngle)
		computeDir(deltaAngle);

	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();
	// Set the camera
	gluLookAt(	x, y, z,
				x+lx, y+ly,  z+lz,
				0.0f, 1.0f,  0.0f);

	// set background
	glClearColor(0.53f, 0.8f, 0.98f, 1.0f);

	// Draw ground
	glColor3f(0.0f, 0.4f, 0.0f);
	glBegin(GL_QUADS);
		glVertex3f(-100.0f, 0.0f, -1000.0f);
		glVertex3f(-100.0f, 0.0f,  100.0f);
		glVertex3f( 100.0f, 0.0f,  100.0f);
		glVertex3f( 100.0f, 0.0f, -1000.0f);
	glEnd();

	// Draw Buildings
	for(int i=-30; i < 0; i++) {
		glPushMatrix();
		glTranslatef(-5.0,0.5,i * 10.0);
		drawBuilding();
		glPopMatrix();
		glPushMatrix();
		glTranslatef(-5.0,1.5,i * 10.0);
		drawBuilding();
		glPopMatrix();
	}
	for(int i=-30; i < 0; i++) {
		glPushMatrix();
		glTranslatef(5.0,0.5,i * 10.0);
		drawBuilding();
		glPopMatrix();
		glPushMatrix();
		glTranslatef(5.0,1.5,i * 10.0);
		drawBuilding();
		glPopMatrix();
	}

	glutSwapBuffers();
}

void processNormalKeys(unsigned char key, int xx, int yy) {

	// escape
	if (key == 27)
		exit(0);
}

void pressKey(int key, int xx, int yy) {

	switch (key) {
		case GLUT_KEY_LEFT : deltaMoveX = -0.01f; deltaAngle = -0.001f; break;
		case GLUT_KEY_RIGHT : deltaMoveX = 0.01f; deltaAngle = 0.001f; break;
		case GLUT_KEY_UP : deltaMoveY = 0.01f; break;
		case GLUT_KEY_DOWN : deltaMoveY = -0.01f; break;
	}
}

void releaseKey(int key, int x, int y) {

	switch (key) {
		case GLUT_KEY_LEFT :
		case GLUT_KEY_RIGHT : 
			deltaMoveX = 0.0f; 
			if (angle > 0)
				deltaAngle = -0.001f; 
			if (angle < 0)
				deltaAngle = 0.001f;
			break;
		case GLUT_KEY_UP :
		case GLUT_KEY_DOWN : deltaMoveY = 0; deltaAngle = 0.0f; break;
	}
}

int main(int argc, char **argv) {

	// init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(320,320);
	glutCreateWindow("StarFox");

	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);

	// handlers
	glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(pressKey);

	// here are the new entries
	glutIgnoreKeyRepeat(1);
	glutSpecialUpFunc(releaseKey);

	// OpenGL init
	glEnable(GL_DEPTH_TEST);

	// texturing
	Image* image = loadBMP("images\\brick.bmp");
	texBrick = loadTexture(image);
	delete image;

	image = loadBMP("images\\stone.bmp");
	texStone = loadTexture(image);
	delete image;

	glEnable(GL_TEXTURE_2D);

	// enter GLUT event processing cycle
	glutMainLoop();

	return 0;
}





