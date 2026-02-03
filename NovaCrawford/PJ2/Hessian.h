#ifndef HESSIAN_H
#define HESSIAN_H
#include <stdio.h>

class Hessian {
public:
    int natom;          // number of atoms
    int ncoord;         // = 3*natom
    double **H;         // full Hessian matrix [3N][3N]

    Hessian(const char *filename);
    //Hessian(FILE *filename);
    ~Hessian();
    void print();
};

#endif