#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "gluvi.h"
#include "gl/glu.h"

//Simple viewer for liquid simulator data
//Hold shift and use the mouse buttons to manipulate the camera

using namespace std;

string frame_number="frame 0";

string file_path;

unsigned int frame= 0;
unsigned int highest_frame = 0;

bool filming = false;
bool running = false;

char * ppmfileformat;

std::vector<Vec3f> particles;
float particle_radius;

bool read_frame(int newframe)
{
   if(newframe<0) newframe = highest_frame;

   std::ostringstream strout;
   
   strout << file_path << "particles_" << newframe << ".txt";
   printf("File path %s\n", strout.str().c_str());
   std::ifstream particles_in(strout.str().c_str());
   if(!particles_in.good()) {
      printf("Failed to open particles!\n");
      return false;
   }
   unsigned int p_count;
   particles_in >> p_count >> particle_radius;
   particles.resize(p_count);
   for(unsigned int p = 0; p < p_count; ++p) {
      Vec3f pos;
      particles_in >> pos[0] >> pos[1] >> pos[2];
      particles[p] = pos;
   }
   printf("Particle count: %d\n", p_count);
   printf("Particle radius: %f\n", particle_radius);

   if(newframe > (int)highest_frame)
      highest_frame = newframe;

   strout.str("");

   frame=newframe;

   strout.str("");
   strout << "frame " << frame;
   frame_number = strout.str();

   return true;
}

void set_view(Gluvi::Target3D &cam)
{
   cam.dist=0.5;
}

void set_lights_and_material(int object)
{
   glEnable(GL_LIGHTING);
   GLfloat global_ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
   glShadeModel(GL_SMOOTH);

   //Light #1
   GLfloat color[4] = {1.0f, 1.0f, 1.0f, 1.0f};
   GLfloat position[3] = {1.0f, 1.0f, 1.0f};
   glLightfv(GL_LIGHT0, GL_SPECULAR, color);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
   glLightfv(GL_LIGHT0, GL_POSITION, position);

   //Light #2
   GLfloat color2[4] = {1.0f, 1.0f, 1.0f, 1.0f};
   GLfloat position2[3] = {-1.0f, -1.0f, 1.0f};
   glLightfv(GL_LIGHT1, GL_SPECULAR, color2);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, color2);
   glLightfv(GL_LIGHT1, GL_POSITION, position2);

   GLfloat obj_color[4] = {.2, .3, .7};
   glMaterialfv (GL_FRONT, GL_AMBIENT, obj_color);
   glMaterialfv (GL_FRONT, GL_DIFFUSE, obj_color);

   GLfloat specular[4] = {.4, .2, .8};
   glMaterialf (GL_FRONT, GL_SHININESS, 32);
   glMaterialfv (GL_FRONT, GL_SPECULAR, specular);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);

}

void timer(int value)
{
   if(filming) {
      Gluvi::ppm_screenshot(ppmfileformat, frame);
      if(read_frame(frame+1)) {
         if(frame == 0) {
            filming = false;
         }
      }
      glutPostRedisplay();
      glutTimerFunc(200, timer, 0);
   }


   if(running) {
      read_frame(frame+1);
      glutTimerFunc(1, timer, 0);
      glutPostRedisplay();
   }

}


void display(void)
{
   glClearColor(0.6f, 0.7f, 0.9f, 1);

   //Coordinate system
   glDisable(GL_LIGHTING);
   glBegin(GL_LINES);
   glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(0.1,0,0);
   glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,0.1,0);
   glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,0.1);
   glEnd();

   glEnable(GL_LIGHTING);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

   set_lights_and_material(1); 

   //Draw the liquid particles as simple spheres for now.
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
   GLUquadric* particle_sphere;
   particle_sphere = gluNewQuadric();
   gluQuadricDrawStyle(particle_sphere, GLU_FILL );
   for(unsigned int p = 0; p < particles.size(); ++p) {
      glPushMatrix();
      Vec3f pos = particles[p].v;
      glTranslatef(pos[0], pos[1], pos[2]);
      gluSphere(particle_sphere, particle_radius, 20, 20);
      glPopMatrix();   
   }
 
   //Draw the bound box for good measure
   glDisable(GL_LIGHTING);
   glColor3f(0,0,0);
   glBegin(GL_LINES);
   glVertex3f(0,0,0);
   glVertex3f(0,0,1);

   glVertex3f(0,0,0);
   glVertex3f(0,1,0);

   glVertex3f(0,0,0);
   glVertex3f(1,0,0);

   glVertex3f(0,1,0);
   glVertex3f(1,1,0);

   glVertex3f(1,1,0);
   glVertex3f(1,0,0);

   glVertex3f(1,0,0);
   glVertex3f(1,0,1);

   glVertex3f(0,1,0);
   glVertex3f(0,1,1);

   glVertex3f(1,1,0);
   glVertex3f(1,1,1);

   glVertex3f(0,1,1);
   glVertex3f(1,1,1);

   glVertex3f(1,0,1);
   glVertex3f(1,1,1);

   glVertex3f(0,0,1);
   glVertex3f(1,0,1);

   glVertex3f(0,0,1);
   glVertex3f(0,1,1);

   glEnd();
   
   //Draw wireframe sphere geometry (specific to this scene).
   glColor3f(0,0,0);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
   GLUquadric* sphere;
   sphere = gluNewQuadric();
   gluQuadricDrawStyle(sphere, GLU_LINE );
   glPushMatrix();
   glTranslatef(0.5f, 0.5f,0.5f);
   gluSphere(sphere, 0.35, 20, 20);
   glPopMatrix();


}

struct ScreenShotButton : public Gluvi::Button{
   const char *filename_format;
   ScreenShotButton(const char *label, const char *filename_format_) : Gluvi::Button(label), filename_format(filename_format_) {}
   void action()
   { 
      Gluvi::ppm_screenshot(filename_format, frame); 
   }
};


struct MovieButton : public Gluvi::Button{
   const char *filename_format;
   MovieButton(const char *label, const char *filename_format_) : Gluvi::Button(label), filename_format(filename_format_) {}
   void action()
   { 
      if(!running) {
         if(!filming) {
            filming = true;
            glutTimerFunc(1000, timer, 0);
         }
         else {
            filming = false;
         }
      }
   }
};

struct RunButton : public Gluvi::Button{
   RunButton(const char *label) : Gluvi::Button(label){}
   void action()
   { 
      if(!filming) {
         if(!running) {
            running = true;
         }
         else {
            running = false;
         }
      }
   }
};

void keyPress(unsigned char key, int x, int y) {

   if(key == 'r') {
      if(!filming) {
         if(!running) {
            running = true;
            glutTimerFunc(1000, timer, 0);
         }
         else {
            running = false;
         }
      }
   }
   else if(key == 'f') {
      if(!running) {
         if(!filming) {
            filming = true;
         }
         else {
            filming = false;
         }
      }
   }
   glutPostRedisplay();



}

void special_key_handler(int key, int x, int y)
{
   int mods=glutGetModifiers();
   switch(key){
case GLUT_KEY_LEFT:
   if(mods == GLUT_ACTIVE_SHIFT) {
      if(read_frame(0))
         glutPostRedisplay();
   }
   else {
      if(read_frame(frame-1))
         glutPostRedisplay();
   }
   break;
case GLUT_KEY_RIGHT:
   if(mods == GLUT_ACTIVE_SHIFT) {
      if(read_frame(highest_frame))
         glutPostRedisplay();
   }
   else {
      if(read_frame(frame+1))
         glutPostRedisplay();
   }
   break;
default:
   ;
   }
}


Gluvi::Target3D* cam_local;
int main(int argc, char **argv)
{

   Gluvi::init("Liquid Data Viewer", &argc, argv);

   if(argc!=2){
      cerr<<"Expecting path to simulation data folder"<<endl;
      return 1;
   }


   file_path=argv[1];   

   //read grid dimensions
   char buffer[100];
   sprintf(buffer, "%s/simdata.txt", file_path.c_str());

   ostringstream strout;
   strout << file_path << "liquidmesh_" << 0 << ".obj";

   //read input mesh here

   if(!read_frame(0))
      return 1;

   glutSpecialFunc(special_key_handler);
   glutKeyboardFunc(keyPress);

   Gluvi::Target3D cam;
   set_view(cam);
   Gluvi::camera=&cam;
   cam_local = &cam;

   Gluvi::userDisplayFunc=display;

   Gluvi::StaticText frametext(frame_number.c_str());
   Gluvi::root.list.push_back(&frametext);

   ppmfileformat = new char[strlen(file_path.c_str())+100];
   sprintf(ppmfileformat, "%sscreenshot%%04d.ppm", file_path.c_str());
   printf("%s\n", ppmfileformat);

   ScreenShotButton screenshot("Screenshot", ppmfileformat);
   Gluvi::root.list.push_back(&screenshot);

   MovieButton movie("Movie", ppmfileformat);
   Gluvi::root.list.push_back(&movie);

   RunButton run("Run");
   Gluvi::root.list.push_back(&run);

   Gluvi::run();
   return 0;
}
