#include "molecule.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>
using namespace std;


Molecule::Molecule(const char *filename, int q)
{
  charge = q;

  // open filename
  std::ifstream is(filename);
  assert(is.good());

  // read the number of atoms from filename
  is >> natom;
  zvals = new int[natom];
  geom = new double* [natom];
  for(int i=0; i < natom; i++)
    geom[i] = new double[3];

  for(unsigned int i=0; i < natom; i++)
    is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  is.close();
}

//clear memory
Molecule::~Molecule()
{
  delete[] zvals;
  for(int i=0; i < natom; i++)
    delete[] geom[i];
  delete[] geom;
}
 
void Molecule::print_geom()
{
  for(int i=0; i < natom; i++)
    printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}
 
//shift coord 
void Molecule::translate(double x, double y, double z)
{
  for(int i=0; i < natom; i++) {
     geom[i][0] += x;
     geom[i][1] += y;
     geom[i][2] += z;
  }
}

//bond length
double Molecule::bond(int i, int j)
{
  return sqrt( pow((geom[i][0]-geom[j][0]), 2.0) + pow((geom[i][1]-geom[j][1]), 2.0) + pow((geom[i][2]-geom[j][2]), 2.0));
}

//unit vector
double Molecule::unit(int cart, int a, int b)
{
  return -(geom[a][cart]-geom[b][cart])/bond(a,b);
}

// angle(i,j,k) in rad
double Molecule::angle(int i, int j, int k)
{
  double ang = unit(0,j,i) * unit(0,j,k) + unit(1,j,i) * unit(1,j,k) + unit(2,j,i) * unit(2,j,k);
  return acos(ang);
}

//Out-of-Plane Angles
double Molecule::oop(int a, int b, int c, int d)
{
  double bcdx = (unit(1,c,b) * unit(2,c,d) - unit(2,c,b) * unit(1,c,d));
  double bcdy = (unit(2,c,b) * unit(0,c,d) - unit(0,c,b) * unit(2,c,d));
  double bcdz = (unit(0,c,b) * unit(1,c,d) - unit(1,c,b) * unit(0,c,d));
  
  double exx = bcdx * unit(0,c,a);
  double eyy = bcdy * unit(1,c,a);
  double ezz = bcdz * unit(2,c,a);

  double theta = (exx + eyy + ezz)/sin(angle(b,c,d));

  if(theta < -1.0) {
    theta = asin(-1.0);
  }
  else if(theta > 1.0) {
    theta = asin(1.0);
  }
  else {
    theta = asin(theta);
  }

  return theta;
}

double Molecule::torsion(int i, int j, int k, int l){
  double jkklx = (unit(1,j,k) * unit(2,k,l) - unit(2,j,k) * unit(1,k,l));
  double jkkly = (unit(2,j,k) * unit(0,k,l) - unit(0,j,k) * unit(2,k,l));
  double jkklz = (unit(0,j,k) * unit(1,k,l) - unit(1,j,k) * unit(0,k,l));

  double ijjkx = (unit(1,i,j) * unit(2,j,k) - unit(2,i,j) * unit(1,j,k));
  double ijjky = (unit(2,i,j) * unit(0,j,k) - unit(0,i,j) * unit(2,j,k));
  double ijjkz = (unit(0,i,j) * unit(1,j,k) - unit(1,i,j) * unit(0,j,k));
  
  double exx1 = jkklx * ijjkx;
  double eyy1 = jkkly * ijjky;
  double ezz1 = jkklz * ijjkz;
  
  double tor = (exx1 + eyy1 + ezz1) / (sin(angle(i,j,k)) * sin(angle(j,k,l)));

  if(tor < -1.0) {
    tor = acos(-1.0);
  }
  else if(tor > 1.0) {
    tor = acos(1.0);
  }
  else {
    tor = acos(tor);
  }

  // sign of the torsion
  double cross_x = ijjky * jkklz - ijjkz * jkkly;
  double cross_y = ijjkz * jkklx - ijjkx * jkklz;
  double cross_z = ijjkx * jkkly - ijjky * jkklx;
  double norm = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
  
  cross_x /= norm;
  cross_y /= norm;
  cross_z /= norm;

  double sign;
  double dot = cross_x * unit(0,j,k) + cross_y*unit(1,j,k) + cross_z*unit(2,j,k);
  
  if(dot < 0.0) {
    sign = -1.0;
  }
  else{
    sign = 1.0;
  }

  return tor*sign;
}

