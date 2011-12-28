#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "array3_utils.h"

using namespace std;

int grid_resolution = 60;
float timestep = 0.01f;
int frame = 0;

float grid_width = 1;

FluidSim sim;

float sphere_phi(const Vec3f& position, const Vec3f& centre, float radius) {
   return (dist(position,centre) - radius);
}

Vec3f c0(0.5f,0.5f,0.5f);
float rad0 = 0.35f;

float boundary_phi(const Vec3f& position) {
   return -sphere_phi(position, c0, rad0);
}

float liquid_phi(const Vec3f& position) {
   return sphere_phi(position, Vec3f(0.55f, 0.55f, 0.4f), 0.23f);
}

void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius);

//Main testing code
//-------------
int main(int argc, char **argv)
{
   
   if(argc!=2){
      cerr << "The first parameter should be the folder to write the output liquid meshes into. (eg. c:\\output\\)" << endl;
      return 1;
   }

   string outpath(argv[1]);
   
   printf("Initializing data\n");
   sim.initialize(grid_width, grid_resolution, grid_resolution, grid_resolution);
   
   printf("Initializing boundary\n");
   sim.set_boundary(boundary_phi);
   
   printf("Initializing liquid\n");
   sim.set_liquid(liquid_phi);
      
   printf("Exporting initial data\n");
   export_particles(outpath, 0, sim.particles, sim.particle_radius);

   for(frame = 1; frame <1000; ++frame) {
      printf("--------------------\nFrame %d\n", frame);
      
      //Simulate
      printf("Simulating liquid\n");
      sim.advance(timestep);
      
      printf("Exporting particle data\n");
      export_particles(outpath, frame, sim.particles, sim.particle_radius);
   }

   return 0;
}


void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius) {
   //Write the output
   
   std::stringstream strout;
   strout << path << "particles_" << frame << ".txt";
   string filepath = strout.str();
   
   ofstream outfile(filepath.c_str());
   //write vertex count and particle radius
   outfile << particles.size() << " " << radius << std::endl;
   //write vertices
   for(unsigned int i = 0; i < particles.size(); ++i)
      outfile << particles[i][0] << " " << particles[i][1] << " " << particles[i][2] << std::endl;
   outfile.close();
}


