#include "GlutTest.h"
#include "R3/R3.h"
#include "R3Scene.h"

using namespace std;

// Arguments
static char *input_scene_name = "../art/level1.scn";

//Network things..
#include <stdlib.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/select.h>    

#define BUFLEN 512
#define NPACK 10
#define SRV_PORT 9930
#define CLI_PORT 9999
#define SRV_IP "140.180.29.64"
#define CLI_IP "140.180.29.64"

pthread_t thr;

static struct sockaddr_in si_me, si_other;
static int sock_in, sock_out;
static socklen_t slen=sizeof(si_other);
static bool two_player = false;
static bool is_server = false;
static bool is_client = false;

struct info_to_send
{
    double xp;
    double yp;
    double zp;
    double health;
};

static struct info_to_send net_info;
static struct info_to_send my_info;

struct ProjectileInter 
{
    int hit;
    R3Point position;
    R3Node *node;
};

// Display Variables
// use these, throw the stuff below away
static R3Scene *scene = NULL;
static R3Camera camera;

// GLUT variables
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;

// START: Varaibles by Awais

double epsilon = 10e-5;
double collision_epsilon = .1;
double cull_depth = 100;
double cull_behind_cutoff = 5;

// this is Arwing
R3Mesh *ship;
R3Node *other_ship;
R3Matrix trans;
double shipTipX,shipTipY,shipTipZ;

// texture variables
GLuint texBrick;
GLuint texStone;

// camera direction
double lx,ly,lz;
// XYZ position of the camera
double x,y,z;

// the key states. These variables will be zero when no key is being presses
double deltaMoveX = 0, deltaMoveZ = 0;
double moveStep = 0.1; // left-right move speed

// angles
// angle of rotation for the camera direction
double deltaAngleLR = 0.0;
double deltaAngleUD = 0.0;

bool rightAngle=false, leftAngle=false; 
bool upAngle=false, downAngle=false;

double angleStep = 0.01;
double angleCutoffR = 0.2;
double angleCutoffL = 0.2;
double angleCutoffU = 0.2;
double angleCutoffD = 0.1;

// camera angles
double cameraAngleLR = 0.0;
double cameraAngleUD = 0.0;

// intersction
bool front_intersection = false;
bool back_intersection = false;
bool left_intersection = false;
bool right_intersection = false;
bool top_intersection = false;
bool bottom_intersection = false;

ProjectileInter projIntersect(R3Node *node, SFProjectile *proj);

R3Point ship_pos = R3Point(0,0,0);

// rotation
double rotationAngle = 0.0;
double rotationStep = 0.01;

// speed variables
double cameraSpeed = 0.05;
double shipSpeed = 0.05;

// mutilple views
enum view {INSIDE, OUTSIDE};
enum view currView = OUTSIDE;

// Health
double health = 100;

// END: Variables by Awais


static double GetTime(void);
double RandomNumber(void);
void setCumulativeTransformations(R3Scene *scene, R3Node *node, R3Matrix transformation);


void diep(char *s)
{
    perror(s);
    exit(1);
}


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
    if (material == current_material) 
        return;
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
              0, 0,  1); 
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


void compute_intersections(R3Point p, R3Point c, double xmid, 
                           double ymid, double zmid)
{
    double front_diff = 0;
    double top_diff = 0;
    double bottom_diff = 0;
    double left_diff = 0;
    double right_diff = 0;
    double back_diff = 0;
    
    front_diff = abs(p.Y() - c.Y() + ymid +  0);
    back_diff = abs(p.Y() - c.Y() - ymid -  0);
    left_diff = abs(p.X() - c.X() + xmid +  0);
    right_diff = abs(p.X() - c.X() - xmid -  0);
    top_diff = abs(p.Z() - c.Z() - zmid -  0 );
    bottom_diff = abs(p.Z() - c.Z() + zmid +  0 );
    
    if (front_diff <= collision_epsilon)
        front_intersection = true;
    if (back_diff <= collision_epsilon)
        back_intersection = true;
    if (left_diff <= collision_epsilon)
        left_intersection = true;
    if (right_diff <= collision_epsilon)
        right_intersection = true;
    if (top_diff <= collision_epsilon)
        top_intersection = true;
    if (bottom_diff <= collision_epsilon)
        bottom_intersection = true;	
}


//Updated this to find when the ship is inside object bboxes
//Also made method adjust node coordinates for transformations
void DrawNode(R3Scene *scene, R3Node *node, R3Matrix transformation)
{
    // Push transformation onto stack
    glPushMatrix();
    LoadMatrix(&node->transformation);
    //Update the transformation
    transformation *= node->transformation;	
    
    //This shows how you would get the *proper* coordinates. 
    //YOU MUST APPLY THE TRANSFORMATION TO THE POINT FIRST.	
    //cout << (transformation * node->bbox.Centroid()).X() << endl;
    
    //Extract the ship's location (from inside a node)
    if (node->shape != NULL &&
        node->shape->mesh == ship)
    {
        R3Point t = transformation * node->shape->mesh->Center();
        ship_pos = t;
        trans = transformation;
        DrawShape(node->shape);
    }
    R3Point box_center = transformation * node->bbox.Centroid();
    double edge_length = (node->bbox.YMax() - node->bbox.YMin())/2;
    //cout << box_center.Y() << " " << y << endl;
   	
    //cout << ship_pos.Y() << endl;
   	
    //Check to see if the ship is inside anything
    if (node->shape != NULL
        &&  node->shape->mesh != ship
        && ship_pos != R3Point(0,0,0))
    {   
        R3Point p = ship_pos;
        R3Point c = R3Point(0,0,0);
        double xmid = 0;
        double ymid = 0;
        double zmid = 0;
        //These are the bbox edge lengths divided by 2
        
        //BBoxes are not explicit-- we have to compute them.
        if (node->shape->mesh != NULL)
        {
            c = transformation * node->shape->mesh->Center();
            xmid = node->shape->mesh->Radius();
            ymid = node->shape->mesh->Radius();
            zmid = node->shape->mesh->Radius();
        }
        else if (node->shape->box != NULL)
        {
            c = transformation * node->shape->box->Centroid();
            xmid = (node->shape->box->XMax() - node->shape->box->XMin())/2;
            ymid = (node->shape->box->YMax() - node->shape->box->YMin())/2;
            zmid = (node->shape->box->ZMax() - node->shape->box->ZMin())/2;
        }
        else if (node->shape->cylinder != NULL)
        {
            c = transformation * node->shape->cylinder->Center();
            xmid = node->shape->cylinder->Radius();
            ymid = node->shape->cylinder->Radius();
            zmid = node->shape->cylinder->Height()/2;
        }
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        
        p.SetY(p.Y() + ship->Radius());
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetY(p.Y() - ship->Radius());
        
        p.SetX(p.X() + ship->Radius());
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetX(p.X() - ship->Radius());
        
        p.SetX(p.X() - ship->Radius());
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetX(p.X() + ship->Radius());
        
        p.SetX(p.X() - ship->Radius()/2);
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetX(p.X() + ship->Radius());
        
        p.SetX(p.X() + ship->Radius()/2);
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetX(p.X() - ship->Radius()/2);
        
        p.SetZ(p.Z() - ship->Radius()/5);
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetZ(p.Z() + ship->Radius()/5);
        
        p.SetZ(p.Z() + ship->Radius()/5);
        //Check if we have an intersection and set the right booleans
        if (p.Z() < (c.Z() + zmid + 0) && (p.Z() > c.Z() - zmid -  0 )
            &&  p.Y() < (c.Y() + ymid +  0) && p.Y() > (c.Y() - ymid -  0)
            &&  p.X() < (c.X() + xmid +  0) && p.X() > (c.X() - xmid -  0))
        {
            compute_intersections(p, c, xmid, ymid, zmid);
            //cout << front_intersection << endl;
        }
        p.SetZ(p.Z() - ship->Radius()/5);
    }
    
    // Load material
    if (node->material) LoadMaterial(node->material);
    
    // Draw shape
    if (node->shape) 
    {
        if (!(box_center.Y() + edge_length < y)
            &&  !(box_center.Y() - edge_length > y + cull_depth))
        {
            DrawShape(node->shape); 
        }
    }
    // Draw children nodes
    for (int i = 0; i < (int) node->children.size(); i++) 
    {
        DrawNode(scene, node->children[i], transformation);
    }
    
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
    //Set the ship position to some value (to be replaced with actual val in DrawNode)
    ship_pos = R3Point(0,0,0);
    //Draw the node-- note we call it with the identity.
    DrawNode(scene, scene->root, R3identity_matrix);
}


void DrawProjectiles(R3Scene *scene)
{
    for (int i = 0; i < scene->NProjectiles(); i++)
    {
        SFProjectile *proj = scene->Projectile(i);
        proj->segment.Draw();
    }
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

static void* receive_data(void *threadid)
{
    while(1)
    {
      	//cout <<"receiving" << endl;
        if (recvfrom(sock_in, &net_info, sizeof(info_to_send), 0, 
                     (struct sockaddr*)&si_other, &slen)==-1)
            diep("recvfrom()");
    }
    pthread_exit(NULL);
}

void GLUTRedraw(void) {
    //fprintf(stderr, "%f\n", ship_pos.Z());
    //fprintf(stderr, "Y\n%f\n%f\n", ship->Center().Y(), ship_pos.Y());
    //fprintf(stderr, "Z\n%f\n%f\n", ship->Center().Z(), ship_pos.Z());
    //double diff = ship->Center().Y() - ship->Face(200)->vertices.at(0)->position.Y();
    //fprintf(stderr, "%f\n", diff-shipTipY);
    //printf("%f:%f:%f\n", ship->Center().X(),ship->Center().Y(),ship->Center().Z());
    
    R3Point old_ship_pos = ship_pos;
    // Awais
    // move camera and the ship forward
    if (deltaMoveX || deltaMoveZ)
        computePos(deltaMoveX, deltaMoveZ);
    
    moveForward();
    updateShip();
    //Reset the intersection
    
    // left-right-up-down movement in viewplane
    // angle movement in LR direction
    if (!left_intersection && !right_intersection)
    {
        if (rightAngle)
            peakRight();
        if (leftAngle)
            peakLeft();
        if (!rightAngle && !leftAngle)
            lookStraightLR();
    }
    // angle movement in UD direction
    if (!top_intersection && !bottom_intersection)
    {
        if (upAngle)
            peakUp();
        if (downAngle)
            peakDown();
        if (!upAngle && !downAngle)
            lookStraightUD();
    }
    if (rotationAngle) 
        computeRotation();
    
    //Update enemies
    updateEnemies();
    //Update projectiles
    updateProjectiles();
    
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
    
    front_intersection = false;
    back_intersection = false;
    left_intersection = false;
    right_intersection = false;
    top_intersection = false;
    bottom_intersection = false;
    
    DrawScene(scene);
    
    
    DrawProjectiles(scene);
    
    // Draw scene edges
    glDisable(GL_LIGHTING);
    glColor3d(1 - background[0], 1 - background[1], 1 - background[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    DrawScene(scene);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    // Write The altitude
    writeAltitude();
    
    //Nader
    //Receive data from companion
    if (two_player)
    {
        if (is_server)//checked
        {
         	//cout << "here" << endl;
         	
            if (sendto(sock_out, &my_info, sizeof(struct info_to_send), 0, 
                       (struct sockaddr*)&si_other, slen)==-1)
                diep("sendto()"); 
         	
         	my_info.xp = (-old_ship_pos.X() + ship_pos.X());
            my_info.zp = (old_ship_pos.Y() - ship_pos.Y());
            my_info.yp = (-old_ship_pos.Z() + ship_pos.Z());
            
         	//cout << (other_ship->transformation * other_ship->shape->mesh->Center()).X() << endl; 
            R3Vector vec = R3Vector(net_info.yp
                                    ,net_info.xp
                                    , net_info.zp
                                    );
            other_ship->transformation.Translate(vec);
            net_info.xp = 0;
         	net_info.yp = 0;
         	net_info.zp = 0;
        }
        if (is_client)//checked
        {
            if (sendto(sock_out, &my_info, sizeof(struct info_to_send), 0, 
                       (struct sockaddr*)&si_other, slen)==-1)
                diep("sendto()");
            
            
            my_info.xp = (ship_pos.X()-old_ship_pos.X());
            my_info.zp = (old_ship_pos.Y() - ship_pos.Y());
            my_info.yp = (-old_ship_pos.Z() + ship_pos.Z());
            
            
            //cout << (other_ship->transformation * other_ship->shape->mesh->Center()).X() << endl; 
         	R3Vector vec = R3Vector(net_info.xp
                                    ,net_info.yp
                                    , net_info.zp
                                    );
            other_ship->transformation.Translate(vec);
            net_info.xp = 0;
         	net_info.yp = 0;
         	net_info.zp = 0;
        }
    }
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
    // spacebar
    else if (key == 32)
        arwingShoot();
    
}

// Awais
// add any special keys over here. Special keys are defined on glut pages
void GLUTSpecial(int key, int xx, int yy) {
    
    switch (key) {
        case GLUT_KEY_LEFT : 
            if(!right_intersection)
            {
                deltaMoveX = -moveStep;
                leftAngle = true;
            }
            break;
        case GLUT_KEY_RIGHT : 
        {
            if (!left_intersection)
            {
                deltaMoveX = moveStep; 
                rightAngle = true;
            }
            break;
        }
        case GLUT_KEY_UP : 
            if (!bottom_intersection)
            {
                deltaMoveZ = moveStep; 
                upAngle = true;
            }
            break;
        case GLUT_KEY_DOWN : 
            if (!top_intersection)
            {
                deltaMoveZ = -moveStep;
                downAngle = true;
            }
            break;
        case GLUT_KEY_F1 :
            rotationAngle = rotationStep;
            break;
        case GLUT_KEY_F2 :
            rotationAngle = -rotationStep;
            break;
        case GLUT_KEY_F3 :
            if (currView == INSIDE) {
                y += -6;
                z += -1;
                currView = OUTSIDE;
            }
            else {
                y += 6;
                z += 1;
                currView = INSIDE;
            }
            break;
    }
}

// Awais
// release key tells what to do upon relese of a key
void releaseKey(int key, int x, int y) {
    
    switch (key) {
        case GLUT_KEY_LEFT :
            deltaMoveX = 0.0; 
            leftAngle = false;
            break;
        case GLUT_KEY_RIGHT : 
            deltaMoveX = 0.0;
            rightAngle = false;
            break;
        case GLUT_KEY_UP :
            deltaMoveZ = 0.0;
            upAngle = false;
            break;
        case GLUT_KEY_DOWN : 
            deltaMoveZ = 0.0; 
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
void computePos(double dx, double dz) {
    double x_move = deltaMoveX;
    double y_move = deltaMoveZ;
    
    if (left_intersection && x_move > 0)
        x_move = 0*x_move;
    if (right_intersection && x_move < 0)
        x_move = 0*x_move;
    if (top_intersection && y_move < 0)
        y_move = 0*y_move;
    if (bottom_intersection && y_move > 0)
        y_move = 0*y_move;
    
    x += x_move;
    z += y_move;	
}

// Below: fncitons for angle movement
void lookStraightLR() {
    double diff = ship->Center().X() - ship->Face(200)->vertices.at(0)->position.X();
    if ((diff - shipTipX) > epsilon) {
        deltaAngleLR = -angleStep;
        ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
        if (currView == INSIDE) {
            cameraAngleLR += -deltaAngleLR;
            lx = sin(cameraAngleLR);
            ly = cos(cameraAngleLR);
        }
    }       
    else if ((diff - shipTipX) < -epsilon) {
        deltaAngleLR = angleStep;
        ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
        if (currView == INSIDE) {
            cameraAngleLR += -deltaAngleLR;
            lx = sin(cameraAngleLR);
            ly = cos(cameraAngleLR);
        }
    }
}

void lookStraightUD() {
    double diff = ship->Center().Y() - ship->Face(200)->vertices.at(0)->position.Y();
    if ((diff - shipTipY) < -epsilon) {
        deltaAngleUD = -angleStep;
        ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
        if (currView == INSIDE) {
            cameraAngleUD += deltaAngleUD;
            lz = sin(cameraAngleUD);
            ly = cos(cameraAngleUD);
        }
    }       
    else if ((diff - shipTipY) > epsilon) {
        deltaAngleUD = angleStep;
        ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
        if (currView == INSIDE) {
            cameraAngleUD += deltaAngleUD;
            lz = sin(cameraAngleUD);
            ly = cos(cameraAngleUD);
        }
    }
}

void peakRight() {
    double diff = ship->Center().X() - ship->Face(200)->vertices.at(0)->position.X();
    if ((diff - shipTipX) > -angleCutoffR) {
        deltaAngleLR = -angleStep;
        ship->Rotate(deltaAngleLR,R3Line(ship->Center(), camera.towards));
        if (currView == INSIDE) {
            cameraAngleLR += -deltaAngleLR;
            lx = sin(cameraAngleLR);
            ly = cos(cameraAngleLR);
        }
    }
}

void peakLeft() {
    double diff = ship->Center().X() - ship->Face(200)->vertices.at(0)->position.X();
    if ((diff - shipTipX) < angleCutoffL) {
        deltaAngleLR = angleStep;
        ship->Rotate(angleStep,R3Line(ship->Center(), camera.towards));
        if (currView == INSIDE) {
            cameraAngleLR += -deltaAngleLR;
            lx = sin(cameraAngleLR);
            ly = cos(cameraAngleLR);
        }
    }
}

void peakUp() {
    double diff = ship->Center().Y() - ship->Face(200)->vertices.at(0)->position.Y();
    if ((diff - shipTipY) > -angleCutoffU) {
        deltaAngleUD = angleStep;
        ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
        if (currView == INSIDE) {
            cameraAngleUD += deltaAngleUD;
            lz = sin(cameraAngleUD);
            ly = cos(cameraAngleUD);
        }
    }
}

void peakDown() {
    double diff = ship->Center().Y() - ship->Face(200)->vertices.at(0)->position.Y();
    if ((diff - shipTipY) < angleCutoffD) {
        deltaAngleUD = -angleStep;
        ship->Rotate(deltaAngleUD,R3Line(ship->Center(), camera.right));
        if (currView == INSIDE) {
            cameraAngleUD += deltaAngleUD;
            lz = sin(cameraAngleUD);
            ly = cos(cameraAngleUD);
        }
    }
}

//Awais
// move forward with a constant speed
void moveForward() {
    //Check if the ship is hitting anything; don't move camera forward if it is.
    if (!front_intersection)
        y += cameraSpeed;
}

void computeRotation(void) {
    ship->Rotate(rotationAngle,R3Line(ship->Center(), camera.up));
}

//Awais
// movement for ship as camera moves
void updateShip() {
    //Check if the ship is hitting anything. Don't move it if it is.
    double x_move = deltaMoveX;
    double y_move = deltaMoveZ;
    double z_move = -shipSpeed;
    
    if (left_intersection &&  x_move > 0)
        x_move = 0* x_move;
    if (right_intersection && x_move < 0)
        x_move = 0*x_move;
    if (top_intersection && y_move < 0)
        y_move = 0*y_move;
    if (bottom_intersection && y_move > 0)
        y_move = 0*y_move;
    if (front_intersection)
        z_move = 0;
    ship->Translate(x_move, y_move, z_move);
}

//Awais
// code adapted from http://www.lighthouse3d.com
// for displaying text on the screen
void renderBitmapString(double x, double y, void *font, char *string) {
    char *c;
    glRasterPos2f(x, y);
    for (c=string; *c != '\0'; c++) {
        glutBitmapCharacter(font, *c);
    }
}

// Awais
// code adapted from http://www.lighthouse3d.com
// used for text display
void setOrthographicProjection() {
    
    // switch to projection mode
    glMatrixMode(GL_PROJECTION);
    
    // save previous matrix which contains the
    //settings for the perspective projection
    glPushMatrix();
    
    // reset matrix
    glLoadIdentity();
    
    // set a 2D orthographic projection
    gluOrtho2D(0, 800, 0, 800);
    
    // invert the y axis, down is positive
    glScalef(1, -1, 1);
    
    // mover the origin from the bottom left corner
    // to the upper left corner
    glTranslatef(0, -800, 0);
    
    // switch back to modelview mode
    glMatrixMode(GL_MODELVIEW);
}

//Awais
// code adapted from http://www.lighthouse3d.com
// used for text display
void restorePerspectiveProjection() {
    
    glMatrixMode(GL_PROJECTION);
    // restore previous projection matrix
    glPopMatrix();
    
    // get back to modelview mode
    glMatrixMode(GL_MODELVIEW);
}

// Awais
// display altitude
void writeAltitude(void) {
    setOrthographicProjection();
    glPushMatrix();
    glLoadIdentity();
    char number1[256];
    sprintf(number1,"%.0f",ship_pos.Z()*10-14);
    char *quant1 = "Altitude: ";
    char *unit1 = " units";
    char *result1 = new char[strlen(quant1)+strlen(number1)+strlen(unit1)];
    sprintf(result1,"%s%s%s",quant1,number1,unit1);
    renderBitmapString(10,20,GLUT_BITMAP_HELVETICA_18,result1);
    glPopMatrix();
    restorePerspectiveProjection();
}

// Awais
// display health
void writeHealth(void) {
	setOrthographicProjection();
	glPushMatrix();
	glLoadIdentity();
	char number1[256];
	sprintf(number1,"%.0f",health);
	char *quant1 = "Health: ";
	char *unit1 = "%";
	char *result1 = new char[strlen(quant1)+strlen(number1)+strlen(unit1)];
	sprintf(result1,"%s%s%s",quant1,number1,unit1);
	renderBitmapString(10,40,GLUT_BITMAP_HELVETICA_18,result1);
	glPopMatrix();
	restorePerspectiveProjection();
}

//Kevin
//have enemies move, shoot
void updateEnemies(void)
{
    for (int i = 0; i < scene->NEnemies(); i++)
    {
        R3Vector y = *(new R3Vector(0,1,0));
        SFEnemy *enemy = scene->Enemy(i);
        
        if (!enemy->fixed)
        {
            //move
            R3Vector *v = &enemy->movementPath;
            enemy->node->shape->mesh->Translate(v->X(), v->Y(), v->Z());
            enemy->position = enemy->node->shape->mesh->Center();
        }
        
        //shoot (will change to a static rate in the future)
        if ((int)GetTime() % 4 == 0 && (int)RandomNumber() % 4 == 0)
        {
            SFProjectile *proj = new SFProjectile(.1);
            R3Point arwingPos = ship_pos + shipSpeed * y;//ship->Center();
            R3Point enemyPos = enemy->position;//enemy->projectileSource;
            
            enemyPos.Transform(scene->Enemy(i)->node->cumulativeTransformation);
            
            //      printf("arwing pos %f:%f:%f\n", arwingPos.X(), arwingPos.Y(), arwingPos.Z());
            
            //      printf("enemy pos %f:%f:%f\n", enemyPos.X(), enemyPos.Y(), enemyPos.Z());
            
            R3Vector projDir = arwingPos - enemyPos;
            
            //       printf("difference vect %f:%f:%f\n", projDir.X(), projDir.Y(), projDir.Z());
            
            //temporary code for simplification. projectile length = 1
            projDir.Normalize();
            
            //    projDir *= 4;
            
            proj->segment = *(new R3Segment(enemyPos, projDir));
            
            scene->projectiles.push_back(proj);
            
            
            // Create scene graph node
            R3Shape *shape = new R3Shape();
            shape->type = R3_SEGMENT_SHAPE;  //will change to mesh when we change laser shape to mesh
            shape->box = NULL;
            shape->sphere = NULL;
            shape->cylinder = NULL;
            shape->cone = NULL;
            shape->mesh = NULL;
            shape->segment = &proj->segment;
            
            R3Node *projNode = new R3Node();
            //note, this material should be the laser material later
            projNode->material = enemy->node->material;
            projNode->shape = shape;
            projNode->bbox = proj->segment.BBox();
            
            //Note: specifically choosing NOT to merge bboxes with the generating enemy
            enemy->node->children.push_back(projNode);
            projNode->parent = enemy->node;
        }
    }
}


//move projectiles, check for intersections
void updateProjectiles(void)
{
    vector<int> deletionIndices;
    for (int i = 0; i < scene->NProjectiles(); i++)
    {
        SFProjectile *proj = scene->Projectile(i);
        
        if (proj->segment.Start().Y() > ship_pos.Y() + cull_depth
            || proj->segment.Start().Y() < ship_pos.Y() - cull_behind_cutoff)
        {
            //projectile beyond active area, so delete
            deletionIndices.push_back(i);
        }
        else
        {
            //calcluate for intersection with object on next move
            
            //move
            //    printf("projectile endpoint %f:%f:%f\n", proj->segment.End().X(), proj->segment.End().Y(), proj->segment.End().Z());
            
            proj->segment.Reset(proj->segment.Start() + proj->speed * proj->segment.Vector(), 
                                proj->segment.End() + proj->speed * proj->segment.Vector());
            
            //    printf(" to %f:%f:%f\n",proj->segment.Start().X(), proj->segment.Start().Y(), proj->segment.Start().Z());
        }
        sort (scene->projectiles.begin(), scene->projectiles.end());
        for (int i = deletionIndices.size() - 1; i >= 0; i--)
        {
            //   delete scene->Projectile(deletionIndices[i]);
            scene->projectiles.erase(scene->projectiles.begin() + deletionIndices[i]);
        } 
    }
}

ProjectileInter projIntersect(R3Node *node, SFProjectile *proj)
{
    
}

void arwingShoot(void)
{
    R3Vector toLeftWing = *(new R3Vector(-1.5,-.5,-.5));
    R3Vector toRightWing = *(new R3Vector(1.5,-.5,-.5));
    
    SFProjectile *left = new SFProjectile(.3);
    SFProjectile *right = new SFProjectile(.3);
    
    //Kevin - I would like to get these shooting in the ship towards
    //rather than the camera, but this will do
    left->segment = *(new R3Segment(ship_pos + toLeftWing, 2 * camera.towards));
    right->segment = *(new R3Segment(ship_pos + toRightWing, 2 * camera.towards));
    
    scene->projectiles.push_back(left);
    scene->projectiles.push_back(right);
    
    
    // Create scene graph node for left laser
    R3Shape *shape = new R3Shape();
    shape->type = R3_SEGMENT_SHAPE;  //will change to mesh when we change laser shape to mesh
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = &left->segment;
    
    R3Node *leftNode = new R3Node();
    //note, this material should be the laser material later
    leftNode->material = scene->arwingNode->material;
    leftNode->shape = shape;
    leftNode->bbox = left->segment.BBox();
    
    // Create scene graph node for right laser
    shape = new R3Shape();
    shape->type = R3_SEGMENT_SHAPE;  //will change to mesh when we change laser shape to mesh
    shape->box = NULL;
    shape->sphere = NULL;
    shape->cylinder = NULL;
    shape->cone = NULL;
    shape->mesh = NULL;
    shape->segment = &right->segment;
    
    R3Node *rightNode = new R3Node();
    //note, this material should be the laser material later
    rightNode->material = scene->arwingNode->material;
    rightNode->shape = shape;
    rightNode->bbox = right->segment.BBox();
    
    //Note: specifically choosing NOT to merge bboxes with the arwing
    scene->arwingNode->children.push_back(leftNode);
    scene->arwingNode->children.push_back(rightNode);
    leftNode->parent = scene->arwingNode;
    rightNode->parent = scene->arwingNode;
}

static void set_up_socket()
{
    if (is_server) //checked
    {
        if ((sock_in=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
            diep("socket");
        
        memset((char *) &si_me, 0, sizeof(si_me));
        si_me.sin_family = AF_INET;
        si_me.sin_port = htons(SRV_PORT);
        si_me.sin_addr.s_addr = htonl(INADDR_ANY);
        //cout << "before bind" << endl;
        if (bind(sock_in, (struct sockaddr*)&si_me, sizeof(si_me))==-1)
            diep("bind");
        
        
        if ((sock_out=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
            diep("socket");
        
        bzero(&si_other, sizeof(si_other));
        //memset((char *) &si_other, 0, sizeof(si_other));
        si_other.sin_family = AF_INET;
        si_other.sin_port = htons(CLI_PORT);
        if (inet_aton(CLI_IP, &si_other.sin_addr)==0) 
        {
            fprintf(stderr, "inet_aton() failed\n");
            exit(1);
        }
    }
   	
    if (is_client) //checked
    {
        if ((sock_out=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
            diep("socket");
        
        bzero(&si_other, sizeof(si_other));
        //memset((char *) &si_other, 0, sizeof(si_other));
        si_other.sin_family = AF_INET;
        si_other.sin_port = htons(SRV_PORT);
        if (inet_aton(SRV_IP, &si_other.sin_addr)==0) 
        {
            fprintf(stderr, "inet_aton() failed\n");
            exit(1);
        }
        
        
        if ((sock_in=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
            diep("socket");
        
        memset((char *) &si_me, 0, sizeof(si_me));
        si_me.sin_family = AF_INET;
        si_me.sin_port = htons(CLI_PORT);
        si_me.sin_addr.s_addr = htonl(INADDR_ANY);
        //cout << "before bind" << endl;
        if (bind(sock_in, (struct sockaddr*)&si_me, sizeof(si_me))==-1)
            diep("bind");
        
    }
    
   	if (two_player)
   	{
        int rc = pthread_create(&thr, NULL, receive_data, NULL);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
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
    
    if (*argc == 2)
    {	
        two_player = true;
        if (strcmp(argv[1], "-s") == 0)
            is_server = true;
        if (strcmp(argv[1], "-c") == 0)
            is_client = true;
    }
    
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
    
   	//Set up a listening socket
    set_up_socket();
   	
   	
   	
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
    setCumulativeTransformations(scene, scene->root, scene->root->transformation);
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
    
    if (two_player)
    {
        if (is_server) //checked
        {
            ship = scene->Root()->children.at(0)->children.at(0)->shape->mesh;
            DrawScene(scene);   
            
            other_ship = new R3Node();
            other_ship->shape = new R3Shape();
            other_ship->shape->type = R3_MESH_SHAPE;
            other_ship->shape->mesh = new R3Mesh(*ship);
            other_ship->children = vector<R3Node*>();  
            other_ship->parent = scene->Root()->children[0];
            trans.Translate(R3Vector(0,2,0));
            other_ship->transformation = trans;
            trans.Translate(R3Vector(0,2,0));
            other_ship->bbox = scene->Root()->children[0]->children[0]->bbox;
            other_ship->material = NULL;
            other_ship->enemy = new SFEnemy();
            
            scene->Root()->children.push_back(other_ship);
        }
        if (is_client) //checked
        {
            ship = scene->Root()->children.at(0)->children.at(0)->shape->mesh;
            DrawScene(scene);   
            
            other_ship = new R3Node();
            other_ship->shape = new R3Shape();
            other_ship->shape->type = R3_MESH_SHAPE;
            other_ship->shape->mesh = new R3Mesh(*ship);
            other_ship->children = vector<R3Node*>();  
            other_ship->parent = scene->Root()->children[0];
            other_ship->transformation = trans;
            scene->Root()->children[0]->children[0]->transformation.Translate(R3Vector(0,2,0));
            trans = scene->Root()->children[0]->children[0]->transformation;
            other_ship->bbox = scene->Root()->children[0]->children[0]->bbox;
            other_ship->material = NULL;
            other_ship->enemy = new SFEnemy();
            
            scene->Root()->children.push_back(other_ship);
        }
        
      	net_info.xp = 0;
      	net_info.yp = 0;
      	net_info.zp = 0;
      	
      	my_info.xp = 0;
      	my_info.yp = 0;
        my_info.zp = 0;
      	
    }
    
    shipTipX = ship->Center().X() - ship->Face(200)->vertices.at(0)->position.X();
    shipTipY = ship->Center().Y() - ship->Face(200)->vertices.at(0)->position.Y();
    shipTipZ = ship->Center().Z() - ship->Face(200)->vertices.at(0)->position.Z();
    
    
    // Return scene
    return scene;
}

//
void setCumulativeTransformations(R3Scene *scene, R3Node *node, R3Matrix transformation)
{
    // Push transformation onto stack
    glPushMatrix();
    LoadMatrix(&node->transformation);
    //Update the transformation
    transformation *= node->transformation;	
    
    //This shows how you would get the *proper* coordinates. 
    //YOU MUST APPLY THE TRANSFORMATION TO THE POINT FIRST.	
    //cout << (transformation * node->bbox.Centroid()).X() << endl;
    
    if (node->enemy != NULL)
    {
        node->cumulativeTransformation *= transformation;
    }
    
    for (int i = 0; i < (int) node->children.size(); i++) 
    {
        setCumulativeTransformations(scene, node->children[i], transformation);
    }
    
    glPopMatrix();
}

////////////////////////////////////////////////////////////
// TIMER CODE
////////////////////////////////////////////////////////////

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>
#endif

static double GetTime(void)
{
#ifdef _WIN32
    // Return number of seconds since start of execution
    static int first = 1;
    static LARGE_INTEGER timefreq;
    static LARGE_INTEGER start_timevalue;
    
    // Check if this is the first time
    if (first) {
        // Initialize first time
        QueryPerformanceFrequency(&timefreq);
        QueryPerformanceCounter(&start_timevalue);
        first = 0;
        return 0;
    }
    else {
        // Return time since start
        LARGE_INTEGER current_timevalue;
        QueryPerformanceCounter(&current_timevalue);
        return ((double) current_timevalue.QuadPart - 
                (double) start_timevalue.QuadPart) / 
        (double) timefreq.QuadPart;
    }
#else
    // Return number of seconds since start of execution
    static int first = 1;
    static struct timeval start_timevalue;
    
    // Check if this is the first time
    if (first) {
        // Initialize first time
        gettimeofday(&start_timevalue, NULL);
        first = 0;
        return 0;
    }
    else {
        // Return time since start
        struct timeval current_timevalue;
        gettimeofday(&current_timevalue, NULL);
        int secs = current_timevalue.tv_sec - start_timevalue.tv_sec;
        int usecs = current_timevalue.tv_usec - start_timevalue.tv_usec;
        return (double) (secs + 1.0E-6F * usecs);
    }
#endif
}
////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double RandomNumber(void)
{
#ifdef _WIN32
    // Seed random number generator
    static int first = 1;
    if (first) {
        srand(GetTickCount());
        first = 0;
    }
    
    // Return random number
    int r1 = rand();
    double r2 = ((double) rand()) / ((double) RAND_MAX);
    return (r1 + r2) / ((double) RAND_MAX);
#else 
    // Seed random number generator
    static int first = 1;
    if (first) {
        struct timeval timevalue;
        gettimeofday(&timevalue, 0);
        srand48(timevalue.tv_usec);
        first = 0;
    }
    
    // Return random number
    return 100 * drand48();
#endif
}

int main(int argc, char **argv) {
    
    GLUTInit(&argc, argv);
    
    // Allocate mesh
    scene = ReadScene(input_scene_name);
    if(!scene) exit(-1);
    
    GLUTMainLoop();
    
    return 0;
}