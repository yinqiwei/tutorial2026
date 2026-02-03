#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "molecule.h"
#include "Hessian.h"
#include "mass.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{
    Molecule mol("PJ2geom.dat", 0);
    Hessian hess("PJ2hess.dat");

    int natom  = mol.natom;
    int ncoord = 3 * natom;

    Mass mass;

    // Step 1: Mass-weighted Hessian
    // cout << "Atomic masses:\n";
    // for (int i = 0; i < natom; i++) {
    //     double m = mass.printmass(mol.zvals[i]);
    //     cout << "Atom " << i << " (Z=" << mol.zvals[i] << "): mass = " << m << "\n";
    // }
    // cout << "\n";

    MatrixXd H_mw(ncoord, ncoord);

    for (int i = 0; i < ncoord; i++) {
        int atom_i = i / 3;
        double mi = mass.printmass(mol.zvals[atom_i]);

        for (int j = 0; j < ncoord; j++) {
            int atom_j = j / 3;
            double mj = mass.printmass(mol.zvals[atom_j]);

            // if (mi <= 0 || mj <= 0) {
            //     cout << "ERROR: Invalid mass at indices " << i << ", " << j 
            //          << ": mi=" << mi << ", mj=" << mj << "\n";
            // }

            H_mw(i,j) = hess.H[i][j] / sqrt(mi * mj);
        }
    }

    // Step 2: Force symmetric
    MatrixXd H_sym = H_mw;
    for (int i = 0; i < ncoord; i++) {
        for (int j = i+1; j < ncoord; j++) {
            double avg = (H_sym(i,j) + H_sym(j,i)) / 2.0;
            H_sym(i,j) = avg;
            H_sym(j,i) = avg;
        }
    }

    // cout << "Mass-weighted Hessian:\n"; 
    // cout << H_sym << endl;
    // cout << endl;

    // int nan_count = 0, inf_count = 0;
    // for (int i = 0; i < ncoord; i++) {
    //     for (int j = 0; j < ncoord; j++) {
    //         if (isnan(H_sym(i,j))) {
    //             nan_count++;
    //             if (nan_count <= 5) cout << "NaN at [" << i << "][" << j << "]\n";
    //         }
    //         if (isinf(H_sym(i,j))) {
    //             inf_count++;
    //             if (inf_count <= 5) cout << "Inf at [" << i << "][" << j << "]\n";
    //         }
    //     }
    // }
    // cout << "Matrix contains " << nan_count << " NaN and " << inf_count << " Inf values\n";

    // Step 3: Diagonalize 
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(H_sym);
    
    if (solver.info() != Eigen::Success) {
        cerr << "Eigenvalue decomposition failed!\n";
        return 1;
    }

    VectorXd evals = solver.eigenvalues();

    cout << "Eigenvalues (a.u.):\n";
    for (int k = ncoord-1; k >= 0; k--)
        cout << setw(3) << k << "   " << setprecision(10) << evals[k] << endl;

    // Step 4: Convert to frequencies (cm^-1)
    cout << "\nHarmonic Vibrational Frequencies (cm^-1):\n";

    // sqrt(Hartree/(amu*bohr^2)) -> cm^-1
    double conv = 5140.484; 
    double tol  = 1e-15;

    for (int k = ncoord-1; k >= 0; k--) {
        double lam = evals[k];
        double freq;

        if (lam > tol)
            freq = sqrt(lam) * conv;
        else if (lam < -tol)
            freq = sqrt(-lam) * conv;
        else
            freq = 0.0;

        cout << setw(3) << k << "   " << fixed << setprecision(4) << freq;
        if (lam < -tol) cout << " i";
        cout << "\n";
    }

    return 0;
}

//g++ -std=c++17 PJ2-3.cc molecule.cc mass.cc Hessian.cc -I/opt/homebrew/include/eigen3 -o PJ2-3