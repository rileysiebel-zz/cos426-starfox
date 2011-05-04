#include "GlutTest.h"
#include "R3/R3.h"
#include "R3Scene.h"

using namespace std;

// Arguments
static char *input_scene_name = "../art/level1.scn";

// Display Variables
// use these, throw the stuff below away
static R3Scene *scene = NULL;
static R3Camera camera;

// GLUT variables
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;

// START: Varaibles by Awais

// this is Arwing
R3Mesh *ship;

// texture variables
GLuint texBrick;
GLuint texStone;

// camera direction
float lx,ly,lz;
// XYZ position of the camera
float x,y,z;

// the key states. These variables will be zero when no key is being presses
float deltaMoveX = 0, deltaMoveZ = 0;
float moveStep = 0.1; // left-right move speed

// angles
// angle of rotation for the camera direction
float angleLR = 0.0f;
float deltaAngleLR = 0.0f;
float angleUD = 0.0f;
float deltaAngleUD = 0.0f;

bool rightAngle=false, leftAngle=false; 
bool upAngle=false, downAngle=false;

float angleStep = 0.01;
float angleCutOff = 0.2;

// rotation
float rotationAngle = 0.0;

// speed varibales
float cameraSpeed = 0.1;
float shipSpeed = 0.1;
// END: Variables by Awais

void DrawShape(R3Shape *shape)
{
  // Check shape type
  if (shape->type == R3_BOX_SHAPE) shape->box->Draw();
  else if (shape->type == R3_SPHERE_SHAPE) shape->sphere->Draw();
  else if (shape->type == R3_CYLINDER_SHAPE) shape->cylinder->Draw();
  else if (shape->type == R3_CONE_SHAPE) shape->cone->Draw();
  else if (shape->type == R3_MESH_SHAPE) shape->mesh->Draw();
  else if (shape->type == R3_SEGMENT_SHAPE) shape->segment->Draw();
  else fprintf(stderr, "Unrecognized shape type: %d\n", shape->type);
}



void LoadMatrix(R3Matrix *matrix)
{
  // Multiply matrix by top of stack
  // Take transpose of matrix because OpenGL represents vectors with 
  // column-vectors and R3 represents them with row-vectors
  R3Matrix m = matrix->Transpose();
  glMultMatrixd((double *) &m);
}



void LoadMaterial(R3Material *material) 
{
  GLfloat c[4];

  // Check if same as current
  static R3Material *current_material = NULL;
  if (material == current_material) return;
  current_material = material;

  // Compute "opacity"
  double opacity = 1 - material->kt.Luminance();

  // Load ambient
  c[0] = material->ka[0];
  c[1] = material->ka[1];
  c[2] = material->ka[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);

  // Load diffuse
  c[0] = material->kd[0];
  c[1] = material->kd[1];
  c[2] = material->kd[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c);

  // Load specular
  c[0] = material->ks[0];
  c[1] = material->ks[1];
  c[2] = material->ks[2];
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c);

  // Load emission
  c[0] = material->emission.Red();
  c[1] = material->emission.Green();
  c[2] = material->emission.Blue();
  c[3] = opacity;
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, c);

  // Load shininess
  c[0] = material->shininess;
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, c[0]);

  // Load texture
  if (material->texture) {
    if (material->texture_index <= 0) {
      // Create texture in OpenGL
      GLuint texture_index;
      glGenTextures(1, &texture_index);
      material->texture_index = (int) texture_index;
      glBindTexture(GL_TEXTURE_2D, material->texture_index); 
      R2Image *image = material->texture;
      int npixels = image->NPixels();
      R2Pixel *pixels = image->Pixels();
      GLfloat *buffer = new GLfloat [ 4 * npixels ];
      R2Pixel *pixelsp = pixels;
      GLfloat *bufferp = buffer;
      for (int j = 0; j < npixels; j++) { 
        *(bufferp++) = pixelsp->Red();
        *(bufferp++) = pixelsp->Green();
        *(bufferp++) = pixelsp->Blue();
        *(bufferp++) = pixelsp->Alpha();
        pixelsp++;
      }
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glTexImage2D(GL_TEXTURE_2D, 0, 4, image->Width(), image->Height(), 0, GL_RGBA, GL_FLOAT, buffer);
      delete [] buffer;
    }

    // Select texture
    glBindTexture(GL_TEXTURE_2D, material->texture_index); 
    glEnable(GL_TEXTURE_2D);
  }
  else {
    glDisable(GL_TEXTURE_2D);
  }

  // Enable blending for transparent surfaces
  if (opacity < 1) {
    glDepthMask(false);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
  }
  else {
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDepthMask(true);
  }
}


void LoadCamera(R3Camera *camera)
{
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(2*180.0*camera->yfov/M_PI, (GLdouble) GLUTwindow_width /(GLdouble) GLUTwindow_height, 0.01, 10000);

  // Set camera transformation
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(x, y, z,
			x+lx, y+ly,  z+lz,
			0, 1,  0); 
}

void LoadLights(R3Scene *scene)
{
  GLfloat buffer[4];

  // Load ambient light
  static GLfloat ambient[4];
  ambient[0] = scene->ambient[0];
  ambient[1] = scene->ambient[1];
  ambient[2] = scene->ambient[2];
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Load scene lights
  for (int i = 0; i < (int) scene->lights.size(); i++) {
    R3Light *light = scene->lights[i];
    int index = GL_LIGHT0 + i;

    // Temporarily disable light
    glDisable(index);

    // Load color
    buffer[0] = light->color[0];
    buffer[1] = light->color[1];
    buffer[2] = light->color[2];
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);

    // Load attenuation with distance
    buffer[0] = light->constant_attenuation;
    buffer[1] = light->linear_attenuation;
    buffer[2] = light->quadratic_attenuation;
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);

    // Load spot light behavior
    buffer[0] = 180.0 * light->angle_cutoff / M_PI;
    buffer[1] = light->angle_attenuation;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glLightf(index, GL_SPOT_EXPONENT, buffer[1]);

    // Load positions/directions
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Load direction
      buffer[0] = -(light->direction.X());
      buffer[1] = -(light->direction.Y());
      buffer[2] = -(light->direction.Z());
      buffer[3] = 0.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_POINT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);
    }
    else if (light->type == R3_SPOT_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
    else if (light->type == R3_AREA_LIGHT) {
      // Load position
      buffer[0] = light->position.X();
      buffer[1] = light->position.Y();
      buffer[2] = light->position.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_POSITION, buffer);

      // Load direction
      buffer[0] = light->direction.X();
      buffer[1] = light->direction.Y();
      buffer[2] = light->direction.Z();
      buffer[3] = 1.0;
      glLightfv(index, GL_SPOT_DIRECTION, buffer);
    }
     else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->type);
      return;
    }

    // Enable light
    glEnable(index);
  }
}



void DrawNode(R3Scene *scene, R3Node *node)
{
  // Push transformation onto stack
  glPushMatrix();
  LoadMatrix(&node->transformation);

  // Load material
  if (node->material) LoadMaterial(node->material);

  // Draw shape
  if (node->shape) DrawShape(node->shape);

  // Draw children nodes
  for (int i = 0; i < (int) node->children.size(); i++) 
    DrawNode(scene, node->children[i]);

  // Restore previous transformation
  glPopMatrix();
}



void DrawLights(R3Scene *scene)
{
  // Check if should draw lights

  // Setup
  GLboolean lighting = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);

  // Draw all lights
  double radius = scene->bbox.DiagonalRadius();
  for (int i = 0; i < scene->NLights(); i++) {
    R3Light *light = scene->Light(i);
    glColor3d(light->color[0], light->color[1], light->color[2]);
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Draw direction vector
      glLineWidth(5);
      glBegin(GL_LINES);
      R3Point centroid = scene->bbox.Centroid();
      R3Vector vector = radius * light->direction;
      glVertex3d(centroid[0] - vector[0], centroid[1] - vector[1], centroid[2] - vector[2]);
      glVertex3d(centroid[0] - 1.25*vector[0], centroid[1] - 1.25*vector[1], centroid[2] - 1.25*vector[2]);
      glEnd();
      glLineWidth(1);
    }
    else if (light->type == R3_POINT_LIGHT) {
      // Draw sphere at point light position
      R3Sphere(light->position, 0.1 * radius).Draw();
    }
    else if (light->type == R3_SPOT_LIGHT) {
      // Draw sphere at point light position and line indicating direction
      R3Sphere(light->position, 0.1 * radius).Draw();
  
      // Draw direction vector
      glLineWidth(5);
      glBegin(GL_LINES);
      R3Vector vector = radius * light->direction;
      glVertex3d(light->position[0], light->position[1], light->position[2]);
      glVertex3d(light->position[0] + 0.25*vector[0], light->position[1] + 0.25*vector[1], light->position[2] + 0.25*vector[2]);
      glEnd();
      glLineWidth(1);
    }
    else if (light->type == R3_AREA_LIGHT) {
      // Draw circular area
      R3Vector v1, v2;
      double r = light->radius;
      R3Point p = light->position;
      int dim = light->direction.MinDimension();
      if (dim == 0) { v1 = light->direction % R3posx_vector; v1.Normalize(); 
	v2 = light->direction % v1; }
      else if (dim == 1) { v1 = light->direction % R3posy_vector; v1.Normalize(); 
	v2 = light->direction % v1; }
      else { v1 = light->direction % R3posz_vector; v1.Normalize(); v2 = light->direction % v1; }
      glBegin(GL_POLYGON);
      glVertex3d(p[0] +  1.00*r*v1[0] +  0.00*r*v2[0], p[1] +  1.00*r*v1[1] +  0.00*r*v2[1], p[2] +  1.00*r*v1[2] +  0.00*r*v2[2]);
      glVertex3d(p[0] +  0.71*r*v1[0] +  0.71*r*v2[0], p[1] +  0.71*r*v1[1] +  0.71*r*v2[1], p[2] +  0.71*r*v1[2] +  0.71*r*v2[2]);
      glVertex3d(p[0] +  0.00*r*v1[0] +  1.00*r*v2[0], p[1] +  0.00*r*v1[1] +  1.00*r*v2[1], p[2] +  0.00*r*v1[2] +  1.00*r*v2[2]);
      glVertex3d(p[0] + -0.71*r*v1[0] +  0.71*r*v2[0], p[1] + -0.71*r*v1[1] +  0.71*r*v2[1], p[2] + -0.71*r*v1[2] +  0.71*r*v2[2]);
      glVertex3d(p[0] + -1.00*r*v1[0] +  0.00*r*v2[0], p[1] + -1.00*r*v1[1] +  0.00*r*v2[1], p[2] + -1.00*r*v1[2] +  0.00*r*v2[2]);
      glVertex3d(p[0] + -0.71*r*v1[0] + -0.71*r*v2[0], p[1] + -0.71*r*v1[1] + -0.71*r*v2[1], p[2] + -0.71*r*v1[2] + -0.71*r*v2[2]);
      glVertex3d(p[0] +  0.00*r*v1[0] + -1.00*r*v2[0], p[1] +  0.00*r*v1[1] + -1.00*r*v2[1], p[2] +  0.00*r*v1[2] + -1.00*r*v2[2]);
      glVertex3d(p[0] +  0.71*r*v1[0] + -0.71*r*v2[0], p[1] +  0.71*r*v1[1] + -0.71*r*v2[1], p[2] +  0.71*r*v1[2] + -0.71*r*v2[2]);
      glEnd();
    }
    else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->type);
      return;
    }
  }

  // Clean up
  if (lighting) glEnable(GL_LIGHTING);
}


void DrawScene(R3Scene *scene) 
{
  // Draw nodes recursively
  DrawNode(scene, scene->root);
}

void GLUTResize(int w, int h) {

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

void GLUTRedraw(void) {

  // Awais
  // move camera and the ship forward
  moveForward();
  updateShip();
  // left-right-up-down movement in viewplane
  if (deltaMoveX || deltaMoveZ)
	computePos(deltaMoveX, deltaMoveZ);
  // angle movement in LR direction
  if (rightAngle)
	  peakRight();
  if (leftAngle)
	  peakLeft();
  if (!rightAngle && !leftAngle)
	  lookStraightLR();
  // angle movement in UD direction
  if (upAngle)
	  peakUp();
  if (downAngle)
	  peakDown();
  if (!upAngle && !downAngle)
	  lookStraightUD();

  if (rotationAngle) 
	  computeRotation();

  // Initialize OpenGL drawing modes
  glEnable(GL_LIGHTING);
  glDisable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDepthMask(true);

  // Clear window 
  R3Rgb background = scene->background;
  glClearColor(background[0], background[1], background[2], background[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Load camera
  LoadCamera(&camera);

  // Load scene lights
  LoadLights(scene);

  // Draw scene camera
  //DrawCamera(scene);

  // Draw scene lights
  //DrawLights(scene);

  // Draw scene surfaces
  glEnable(GL_LIGHTING);
  DrawScene(scene);
  

  // Draw scene edges
  glDisable(GL_LIGHTING);
  glColor3d(1 - background[0], 1 - background[1], 1 - background[2]);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  DrawScene(scene);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  glutSwapBuffers();
 
}

// Awais
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

// Awais
// add any normal keyboard keys over here. Normal Keys are defined on glut pages
void GLUTKeyboard(unsigned char key, int xx, int yy) {
	// escape
	if (key == 27)
		exit(0);

}

// Awais
// add any special keys over here. Special keys are defined on glut pages
void GLUTSpecial(int key, int xx, int yy) {

	switch (key) {
		case GLUT_KEY_LEFT : 
			deltaMoveX = -moveStep;
			leftAngle = true;
			break;
		case GLUT_KEY_RIGHT : 
			deltaMoveX = moveStep; 
			rightAngle = true;
			break;
		case GLUT_KEY_UP : 
			deltaMoveZ = moveStep; 
			upAngle = true;
			break;
		case GLUT_KEY_DOWN : 
			deltaMoveZ = -moveStep;
			downAngle = true;
			break;
		case GLUT_KEY_F1 :
			rotationAngle = 0.01;
			break;
		case GLUT_KEY_F2 :
			rotationAngle = -0.01;
			break;
	}
}

// Awais
// release key tells what to do upon relese of a key
void releaseKey(int key, int x, int y) {

	switch (key) {
		case GLUT_KEY_LEFT :
			deltaMoveX = 0.0f; 
			leftAngle = false;
			break;
		case GLUT_KEY_RIGHT : 
			deltaMoveX = 0.0f;
			rightAngle = false;
			break;
		case GLUT_KEY_UP :
			deltaMoveZ = 0;
			upAngle = false;
			break;
		case GLUT_KEY_DOWN : 
			deltaMoveZ = 0; 
			downAngle = false;
			break;
		case GLUT_KEY_F1 :
			rotationAngle = 0.0;
			break;
		case GLUT_KEY_F2 :
			rotationAngle = 0.0;
			break;
	}
}

// Awais
void computePos(float dx, float dz) {
	x += deltaMoveX;
	z += deltaMoveZ;	
}
// Below: fncitons for angle movement
void lookStraightLR() {

	if (angleLR > angleCutOff)
		angleLR = angleCutOff;
	if (angleLR < -angleCutOff)
		angleLR = -angleCutOff;

	if (angleLR > 0) {
		deltaAngleLR = -angleStep;
		angleLR += deltaAngleLR;
		ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
	}
	else if (angleLR < 0) {
		deltaAngleLR = angleStep;
		angleLR += deltaAngleLR;
		ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
	}
	else
		deltaAngleLR = 0.0;
}
void lookStraightUD() {

	if (angleUD > angleCutOff)
		angleUD = angleCutOff;
	if (angleUD < -angleCutOff)
		angleUD = -angleCutOff;

	if (angleUD > 0) {
		deltaAngleUD = -angleStep;
		angleUD += deltaAngleUD;
		ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
	}
	else if (angleUD < 0) {
		deltaAngleUD = angleStep;
		angleUD += deltaAngleUD;
		ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
	}
	else
		deltaAngleUD = 0.0;
}

void peakRight() {
	if (angleLR >= -angleCutOff) {
		deltaAngleLR = -angleStep;
		angleLR += deltaAngleLR;
		ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
	}
}
void peakLeft() {
	if (angleLR <= angleCutOff) {
		deltaAngleLR = angleStep;
		angleLR += deltaAngleLR;
		ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
	}
}

void peakUp() {
	if (angleUD <= angleCutOff) {
		deltaAngleUD = angleStep;
		angleUD += deltaAngleUD;
		ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
	}
}
void peakDown() {
	if (angleUD >= -angleCutOff) {
		deltaAngleUD = -angleStep;
		angleUD += deltaAngleUD;
		ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
	}
}

//Awais
// move forward with a constant speed
void moveForward() {
	y += cameraSpeed;
}

void computeRotation(void) {
	ship->Rotate(rotationAngle,R3Line(ship->Center(), camera.up));
}

//Awais
// movement for ship as camera moves
void updateShip() {
	ship->Translate(deltaMoveX,deltaMoveZ,-shipSpeed);
}

// Riley added this, just took some stuff out of main
// for organizational purposes
void GLUTInit(int *argc, char **argv) {
   // init GLUT and create window
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowPosition(100,100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutCreateWindow("StarFox");
  
  // register callbacks
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutIdleFunc(GLUTRedraw);
  
  // handlers
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  
  // here are the new entries
  glutIgnoreKeyRepeat(1);
  glutSpecialUpFunc(releaseKey);
  
  // OpenGL init
  glEnable(GL_DEPTH_TEST);
  // texturing
  /*	Image* image = loadBMP("images\\brick.bmp");
	texBrick = loadTexture(image);
	delete image;
	
	image = loadBMP("images\\stone.bmp");
	texStone = loadTexture(image);
	delete image; */

  glEnable(GL_TEXTURE_2D);
}

// Riley
// More organizational stuff
// Just in case we want it to be more complex later
void GLUTMainLoop() {
  glutMainLoop();
}

// Riley
// Reading the scene
R3Scene *
ReadScene(const char *filename) {
  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Read file
  if (!scene->Read(filename)) {
    fprintf(stderr, "Unable to read scene from %s\n", filename);
    return NULL;
  }

  // Remember initial camera
  camera = scene->camera;

  // Awais
  // set the global variables from the camera
  x = camera.eye.X();
  y = camera.eye.Y();
  z = camera.eye.Z();
  lx = camera.towards.X();
  ly = camera.towards.Y();
  lz = camera.towards.Z();
  // get the ship
  ship = scene->Root()->children.at(0)->children.at(0)->shape->mesh;

  // Return scene
  return scene;
}

int main(int argc, char **argv) {

  GLUTInit(&argc, argv);

  // Allocate mesh
  scene = ReadScene(input_scene_name);
  if(!scene) exit(-1);
   
  GLUTMainLoop();

  return 0;
}





