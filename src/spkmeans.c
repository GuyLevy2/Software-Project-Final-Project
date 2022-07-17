#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spkmeans.h"

int main(int argc, char *argv[]) {
    
}


/* Main Utility Functions */

int wam_func(double*** vectors_list, int N, int dim, double*** wamMat){

}

int ddg_func(int N, double*** wamMat, double*** ddgMat){

}

int lNorm_func(int N, double*** wamMat, double*** ddgMat, double*** lNormMat){
    
}

int jacobi_func(int N, double*** symMat, double*** eigeVectors, double** eigenValues){

}



/* Granular Utility Functions */

double vectorDist(double* v1, double* v2, int dim){
    double dist = 0;
    int i = 0;
    
    for(i = 0; i < dim; ++i){
        dist += pow(v1[i]-v2[i], 2.0);
    }

    return pow(dist, 0.5);
}

int minusSqrtD(int N, double *** D){

}

int diagMatMult(int N, int diagPosition, double*** mat1, double*** mat2, double*** outputMat){
    
}

int matMult(int N, double*** mat1, double*** mat2, double*** outputMat){

}

int matDup(double*** origMat, double*** newMat){

}

int matTranspose(double*** mat){

}

double*** initMat(int N){

}

int freeMat(int N){
    
}

int convergenceTest(double epsilon, double*** mat1, double*** mat2){

}

int offCalc(int N, double*** mat){

}