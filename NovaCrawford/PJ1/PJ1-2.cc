#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#define _USE_MATH_DEFINES 

#include "molecule.h"
#include "mass.h"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

int main()
{  
    Molecule mol("geom.dat", 0);

    //Step1 Print Out 
    printf ("Number of atoms: ");
    printf ("%d\n", mol.natom);
    printf ("\nInput:\n");
    mol.print_geom();

    //Step2 Bond Length
    printf ("\nBond Length\n");
    for(int i=0; i < mol.natom; i++)
        for(int j=0; j < i; j++)
            printf("%d %d %10.6f\n", i, j, mol.bond(i,j));


    //Step3 Bond Angles
    printf ("\nBond Angles:\n");
    for(int i=0; i < mol.natom; i++) {
        for(int j=0; j < i; j++) {
            for(int k=0; k < j; k++) {
                if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)
                    printf("%2d %2d %2d %10.6f\n", i, j, k, mol.angle(i,j,k)*(180.0/M_PI));
            }
        }
    }

    //Step4 Out-of-Plane Angles
    printf ("\nOut-of-Plane Angles:\n");
    for(int i=0; i < mol.natom; i++) {
        for(int j=0; j < mol.natom; j++) {
            for(int k=0; k < mol.natom; k++) {
                for(int l=0; l < k; l++) {
                    if(i != j && i != k && i != l && j != k && j != l && k != l 
                        && mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(k,l) < 4.0)
                        printf("%2d %2d %2d %2d %10.6f\n", i, j, k, l, mol.oop(i,k,j,l)*(180.0/M_PI));
                        //oop() order matters
                }
            }
        }
    }

    //Step5 Torsion
    printf ("\nTorsion:\n");
    for(int i=0; i < mol.natom; i++) {
        for(int j=0; j < i; j++) {
            for(int k=0; k < j; k++) {
                for(int l=0; l < k; l++) {
                    if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(k,l) < 4.0)
                        printf("%2d %2d %2d %2d %10.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/M_PI));
                }
            }
        }
    }

    //Step6 Center of Mass
    printf ("\nCenter of Mass:\n");
    Mass m;
    double mtot;
    double xm;
    double ym;
    double zm;

    for(int i=0; i < mol.natom; i++) {
        mtot = mtot + m.printmass((int) mol.zvals[i]);
    }

    for(int i=0; i < mol.natom; i++) {
        xm = xm + (m.printmass((int) mol.zvals[i]) *  mol.geom[i][0]);
        ym = ym + (m.printmass((int) mol.zvals[i]) *  mol.geom[i][1]);
        zm = zm + (m.printmass((int) mol.zvals[i]) *  mol.geom[i][2]);
    }

    printf ("%12.8f %12.8f %12.8f\n", xm/mtot, ym/mtot, zm/mtot);

    xm = xm / mtot;
    ym = ym / mtot;
    zm = zm / mtot;

    mol.translate (-xm, -ym, -zm);

    //Step7 Moments of Inertia
    printf ("\nMoments of Inertia:\n");
    Matrix I(3,3);

    for(int i=0; i < 3; i++) {
        //Diagonal
        I(0,0) += m.printmass((int) mol.zvals[i]) * (mol.geom[i][1]*mol.geom[i][1] + mol.geom[i][2]*mol.geom[i][2]);
        I(1,1) += m.printmass((int) mol.zvals[i]) * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][2]*mol.geom[i][2]);
        I(2,2) += m.printmass((int) mol.zvals[i]) * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][1]*mol.geom[i][1]);
        //Off-Diagonal
        I(0,1) -= m.printmass((int) mol.zvals[i]) * (mol.geom[i][0]*mol.geom[i][1]);
        I(0,2) -= m.printmass((int) mol.zvals[i]) * (mol.geom[i][0]*mol.geom[i][2]);
        I(1,2) -= m.printmass((int) mol.zvals[i]) * (mol.geom[i][1]*mol.geom[i][2]);
    }

    // I(0,1) = I(1,0);
    // I(2,0) = I(0,2);
    // I(1,2) = I(2,1);
    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);

    printf ("\nthe moment of inertia tensor:\n");
    cout << I << endl;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    //unit conversions
    printf ("\nPrincipal moments of inertia (amu*bohr^2):\n");
    cout << evals;
    printf ("\n(unit: amu*AA^2):\n");
    cout << evals * (0.529177249 * 0.529177249);
    printf ("\n(unit: g*cm^2):\n");
    cout << evals * (1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8);

    printf ("\n\nClassification of molecular rotors:   ");
    if(mol.natom == 2) 
        printf("diatomic molecule\n");
    else if(evals(0) < 1e-4) 
        printf("linear\n");
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4)) 
    //Ia=Ib=Ic
        printf("spherical top\n");
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(2) - evals(1)) > 1e-4)) 
    // Ia=Ib and Ib<<Ic
        printf("oblate symmetric top\n");
    else if((fabs(evals(1) - evals(0)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
    // Ia<<Ib and Ib=Ic
        printf("prolate symmetric top\n");
    else 
        printf("asymmetric top\n");

    //Step8 Rotational Constants


    //rot constant

    return 0;
}