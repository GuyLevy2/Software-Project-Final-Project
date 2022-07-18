#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spkmeans.h"

int main(int argc, char *argv[]) {
    
}


/* Main Utility Functions */

int wam_func(double*** vectors_list, int N, int dim, double*** wamMat){
    int i;
    int j;
    
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            
            if (i == j){
                (*wamMat)[i][j] = 0;
            }
            
            else{
                double** vec1 = (*vectors_list)[i];
                double** vec2 = (*vectors_list)[j];
                double dist = vectorDist(vec1, vec2, dim);
                (*wamMat)[i][j] = exp(-(dist/2));
            }
        }
    }
}

int ddg_func(int N, double*** wamMat, double*** ddgMat){
    int i;
    int j;
    double sum;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            sum += (*wamMat)[i][j];
        }
        
        (*ddgMat)[i][i] = sum;
        sum = 0;
    }
}

int lNorm_func(int N, double*** wamMat, double*** ddgMat, double*** lNormMat){
    int i;
    int j;
    
    double*** minusSqrtMat = initMat(N);
    matDup(ddgMat, minusSqrtMat);

    minusSqrtD(N, minusSqrtMat);
    
    double*** resMat1 = initMat(N);
    matMult(N, minusSqrtMat, wamMat, resMat1);

    matMult(N, resMat1, minusSqrtMat, lNormMat);
    
    freeMat(N, resMat1);
    freeMat(N, minusSqrtMat);

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                (*lNormMat)[i][j] = 1 - (*lNormMat)[i][j];
            }
            else{
                (*lNormMat)[i][j] = - (*lNormMat)[i][j];
            }
        }
    }
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

void minusSqrtD(int N, double*** D){ // Guy - changed return type to void???
    int i;
    for (i = 0; i < N; i++){
        (*D)[i][i] = 1 / (sqrt((*D)[i][i]));
    }
}

int diagMatMult(int N, int diagPosition, double*** mat1, double*** mat2, double*** outputMat){
    
}

int matMult(int N, double*** mat1, double*** mat2, double*** outputMat){
    int i;
    int j;
    int m;
    int sum = 0;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            for (m = 0; m < N; m++){
                sum += (*mat1)[i][j + m] * (*mat2)[i + m][j];
            }
            (*outputMat)[i][j] = sum;
            sum = 0;
        }
    }
}

int matDup(double*** origMat, double*** newMat){
    int i;
    int j;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            (*newMat)[i][j] = (*origMat)[i][j];
        }
    }
}

int matTranspose(double*** mat){

}

double*** initMat(int N){
    int i;
    int j;
    double** newMat = malloc(N * sizeof(double*)); // Guy - should we check if succeeded???
    
    for (i = 0; i < N; i++){
        double* vec = malloc(N * sizeof(double));
        
        for (j = 0; j < N; j++){
            vec[j] = 0;
        }

        newMat[i] = vec;
    }

    return &newMat;
}

int freeMat(int N, double*** mat){
    int i;

    for (i = 0; i < N; i++){
        free((*mat)[i]);
    }
    
    free(*mat);
}

int convergenceTest(double epsilon, double*** mat1, double*** mat2){

}

int offCalc(int N, double*** mat){

}