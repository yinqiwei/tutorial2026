#include "Hessian.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
using namespace std;

Hessian::Hessian(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Cannot open Hessian file");
        exit(1);
    }

    // First number in file = natom
    if (fscanf(fp, "%d", &natom) != 1) {
        fprintf(stderr, "Error: could not read natom from %s\n", filename);
        fclose(fp);
        exit(1);
    }

    ncoord = 3 * natom;

    // Read all remaining doubles into a flat vector
    vector<double> vals;
    vals.reserve(ncoord * ncoord);

    double x;
    int count = 0;
    while (fscanf(fp, "%lf", &x) == 1) {
        vals.push_back(x);
        count++;
    }

    fclose(fp);

    int expected = ncoord * ncoord;
    if ((int)vals.size() >= expected) {
        vals.resize(expected);
        //cout << "DEBUG: Using first " << expected << " values\n";
    } else {
        cerr << "ERROR: Not enough values in file!\n";
        exit(1);
    }

    // Allocate H
    H = new double*[ncoord];
    for (int i = 0; i < ncoord; i++)
        H[i] = new double[ncoord];

    // Fill H row-major from vals
    int idx = 0;
    for (int i = 0; i < ncoord; i++) {
        for (int j = 0; j < ncoord; j++) {
            H[i][j] = vals[idx++];
        }
    }

    // Debug: print first 3x3 block
    // cout << "DEBUG: First 3x3 block of H:\n";
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         cout << H[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    // Debug: print last 3x3 block
    // cout << "DEBUG: Last 3x3 block of H:\n";
    // for (int i = ncoord-3; i < ncoord; i++) {
    //     for (int j = ncoord-3; j < ncoord; j++) {
    //         cout << H[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    // Check symmetry
    int asymmetric_count = 0;
    double max_asymmetry = 0.0;
    for (int i = 0; i < ncoord; i++) {
        for (int j = i+1; j < ncoord; j++) {
            double diff = fabs(H[i][j] - H[j][i]);
            if (diff > 1e-10) {
                asymmetric_count++;
                if (diff > max_asymmetry) {
                    max_asymmetry = diff;
                    // cout << "DEBUG: H[" << i << "][" << j << "] = " << H[i][j] 
                    //      << ", H[" << j << "][" << i << "] = " << H[j][i]
                    //      << ", diff = " << diff << "\n";
                }
            }
        }
    }
    //cout << "DEBUG: Found " << asymmetric_count << " asymmetric pairs, max diff = " << max_asymmetry << "\n";
}

Hessian::~Hessian()
{
    for (int i = 0; i < ncoord; i++)
        delete[] H[i];
    delete[] H;
}

void Hessian::print()
{
    for (int i = 0; i < ncoord; i++) {
        for (int j = 0; j < ncoord; j++)
            printf("%12.6f ", H[i][j]);
        printf("\n");
    }
}