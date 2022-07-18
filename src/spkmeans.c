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

/* 
 * Function: jacobi_func
 * ---------------------
 * computes the eigenvectors and eigenvalues of a given real, symmetric and full
 *      rank matrix using jacobi algorithm
 * 
 * N: the dimention of the matrix symMat (NxN)
 * symMat: a pointer to the given real, symmetric and full rank matrix that
 *      arranged as two dimensional array
 * eigenVectors: a pointer to the output eigenvectors that arranged as two
 *      dimensional array
 * eigenValues: a pointer to the output eigenvalues that arranged as one 
 *      dimensional array
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int jacobi_func(int N, double*** symMat, double*** eigenVectors, double** eigenValues){

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

/* 
 * Function: matMult
 * ---------------------
 * computes the multiplication of two given nxn matrices
 * 
 * N: the dimention of both given matrices (NxN)
 * mat1: a pointer to the first given matrix
 * mat2: a pointer to the second given matrix
 *      dimensional array
 * outputMat: a pointer to the output multipling matrix mat1*mat2 
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
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
    /* Guy
    int i;
    int j;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            (*newMat)[i][j] = (*origMat)[i][j];
        }
    }
    */

    int i,j,k;
    double*** mat2_T = NULL;
    double m_ij;

    mat2_T = initMat(N);
    if (mat2_T == NULL){
        return 1;
    }

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            m_ij = 0;
            for(k = 0; k < N; k++){
                m_ij += ((*mat1)[i][k]) * ((*mat2_T)[j][k]);    
            }
            (*outputMat)[i][j] = m_ij;
        }
    }

    if(freeMat(N,mat2_T)){
        return 1;
    }
    return 0;
}

/* 
 * Function: matDup
 * ---------------------
 * creates a duplication of a given matrix
 * 
 * N: the dimention of the given matrix (NxN)
 * origMat: a pointer to the given matrix
 * newMat: a pointer to the output dulicated matrix (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int matDup(int N, double*** origMat, double*** newMat){
    int i,j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            (*newMat)[i][j] = (*origMat)[i][j];
        }
    }
    return 0;
}

/* 
 * Function: matTranspose
 * ---------------------
 * computes a transpose matrix of a given nxn matrix
 * 
 * N: the dimention of the given matrix (NxN)
 * mat: a pointer to the given matrix
 * matT: a pointer to the output transpose matrix (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int matTranspose(int N, double*** mat, double*** matT){
    int i,j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            matT[j][i] = mat[i][j];
        }
    }

    return 0;
}

/* 
 * Function: computeA_tag
 * ---------------------
 * computes the A' matrix defind as P^T*A*P where A is a real symmetric
 *      matrix and P is a rotation matrix with c and s as rotation elemens:
 *          (1
 *             ...
 *                  c  ...  s
 *                  .       .
 *      P =         .   1   .
 *                  .       .
 *                  -s ...  c
 *                            ...
 *                                1)
 * 
 * N: the dimention of the given matrix A (NxN)
 * i,j: the row and colum which A_ij is the off-diagonal element with the
 *      largenst absolute value
 * A: a poiner to the given real and symmenric matrix
 * c,s: the rotation elements of the rotation matrix P
 * A_tag: a pointer to the output A' matrix (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int computeA_tag(int N, int i, int j, double*** A, double c, double s, double*** A_tag){
    int r;
    double a_tag_ri, a_tag_rj;
    double s_squared, c_squared;
    double a_ii, a_jj, a_ij;

    matDup(N, A, A_tag);

    for(r = 0; r < N; r++){
        if(r != i && r != j){
            a_tag_ri = c*((*A)[r][i]) - s*((*A)[r][j]);
            a_tag_rj = c*((*A)[r][j]) + s*((*A)[r][i]);

            (*A_tag)[r][i] = a_tag_ri;
            (*A_tag)[i][r] = a_tag_ri;
            
            (*A_tag)[r][j] = a_tag_rj;
            (*A_tag)[j][r] = a_tag_rj;
        }
    }

    s_squared = pow(s, 2.0);
    c_squared = pow(c, 2.0);
    a_ii = (*A)[i][i];
    a_jj = (*A)[j][j];
    a_ij = (*A)[i][j];

    (*A_tag)[i][i] = c_squared*a_ii + s_squared*a_jj - 2*s*c*a_ij;
    (*A_tag)[j][j] = s_squared*a_ii + c_squared-a_jj + 2*s*c*a_ij;
    (*A_tag)[i][j] = 0;
    (*A_tag)[j][i] = 0;

    return 0;
}

/* 
 * Function: convergenceTest
 * ---------------------
 * determine if there is a convergence during the iterations of jacobi 
 *      algorithm
 * 
 * N: the dimention of both given matrices (NxN)
 * epsilon: the criterion of convergence
 * mat1: a pointer to the first matrix
 * mat2: a pointer to the second matrix
 * 
 * returns: off(mat1)^2 - off(mat2)^2 <= epsilon
 */
bool convergenceTest(int N, double epsilon, double*** mat1, double*** mat2){
    return (offCalc(N,mat1) - offCalc(N,mat2)) <= epsilon;
}

/* 
 * Function: offCalc
 * ---------------------
 * computes the sum squares of all off-diagonal elements of a given matrix
 *  
 * N: the dimention of the given matrix (NxN)
 * mat: a pointer to the matrix
 * 
 * returns: off(mat)^2  as described above
 */
double offCalc(int N, double*** mat){
    double off = 0;
    int i,j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            if(j != i){
                off += pow((*mat)[i][j], 2.0);
            }
        }
    }

    return off;
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