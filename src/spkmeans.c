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
 *      dimensional array (initiailzed with zeroes)
 * eigenValues: a pointer to the output eigenvalues that arranged as one 
 *      dimensional array (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int jacobi_func(int N, double*** symMat, double*** eigenVectors, double** eigenValues){
    int MAX_ROTATIONS = 100;
    int i = 0, j = 0, countRot = 0, k;
    double c = 0, s = 0;
    double*** P = NULL, ***A = NULL, ***A_tag = NULL, ***tempMat = NULL;
    double EPSILON = 0.00001;

    unityMat(N, eigenVectors);  /* eigenVectors = I(N) */
    
    A = initMat(N);
    if (A == NULL){
        return 1;
    }
    
    A_tag = initMat(N);
    if (A_tag == NULL){
        return 1;
    }
    
    tempMat = initMat(N);
    if (tempMat == NULL){
        return 1;
    }

    matDup(N, symMat, A_tag);   /* A' = symMat (by values) */

    P = initMat(N);
    if (P == NULL){
        return 1;
    }

    do{
        matDup(N, A_tag, A);    /* A = A' (by values)    */
         
        if (find_ij_pivot(N, A, &i, &j)){       
            return 1;
        }

        unityMat(N, P);                         /* P = I(N)                      */
        buildRotMat(N, A, i, j, &c, &s, P);     /* P is rotate matrix wrt A and c and s are updeted     */ /* Liad - maybe to split into 2 diff functions: computes_c&s + bulidP? */
        matMult(N, eigenVectors, P, tempMat);   /* tempMat = eigenVectors * P    */
        matDup(N, tempMat, eigenVectors);       /* eigenVectors = tempMat        */
        computeA_tag(N, i, j, A, c, A, A_tag);  /* A' = P^T*A*P                  */

        countRot++;
    } while ((!convergenceTest(N, EPSILON, A, A_tag)) && (countRot < MAX_ROTATIONS));

    /* Extract the eigenValues from the diagonal of A' */
    for(k = 0; k < N; k++){
        (*eigenValues)[k] = (*A_tag)[k][k];
    }
    
    /* Free all matrixes which allocated during this function */
    if (freeMat(N, P) || freeMat(N, A) || freeMat(N, A_tag) || freeMat(N, tempMat)){
        return 1;
    }
    
    return 0;
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

void minusSqrtD(int N, double*** D){ // Guy - changed return type to void??? Liad - I think we shuld remain the return type as int for scalability
    int i;
    for (i = 0; i < N; i++){
        (*D)[i][i] = 1 / (sqrt((*D)[i][i]));
    }
}

int diagMatMult(int N, int diagPosition, double*** mat1, double*** mat2, double*** outputMat){
    
}

/* 
 * Function: unityMat
 * ------------------
 * sets a given matrix to be unity matrix
 * 
 * N: the dimention of the given matrix mat (NxN)
 * mat: a pointer to the output unity matrix
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int unityMat(int N, double*** mat){
    int i, j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            (*mat)[i][j] = i == j ? 1 : 0;
        }
    }

    return 0;
}

/* 
 * Function: buildRotMat
 * ---------------------
 * builds the rotation matrix P using a given symmetric matrix A &
 *      computes s and c values
 * 
 * N: the dimention of both given matrices A,P (NxN)
 * A: a pointer to the given symmetric matrix
 * i,j: the indexes of the pivot element A_ij of the matrix A
 * c_p, s_p: pointers to the output values of c and s
 * P: a pointer to the output rotation matrix (initializes as unity matrix)
 * 
 * returns: 0 if there is no exception and 1 elsewhere
 */
int buildRotMat(int N, double*** A, int i, int j, int* c_p, int* s_p, double*** P){
    int sign_theta;
    double theta, abs_theta, t, c, s;

    theta = (((*A)[j][j]) - ((*A)[i][i]))) / (2*((*A)[i][j]));
    abs_theta = fabs(theta);
    sign_theta = theta >= 0 ? 1 : -1;

    t = (sign_theta*1.0) / (abs_theta + pow(pow(theta, 2.0) + 1, 0.5));

    c = 1.0 / pow(pow(t, 2.0) + 1, 0.5);
    s = t*c;

    *c_p = c;
    *s_p = s;
    (*P)[i][i] = c;
    (*P)[i][j] = s;
    (*P)[j][j] = c;
    (*P)[j][i] = -1.0*s;
    /* TODO: which is i and j in the matrix P??? */

    return 0;
}

/* 
 * Function: find_ij_pivot
 * -----------------------
 * computes the indexes i,j of the pivot element A_ij
 * 
 * The pivot A_ij is the off-digonal element with the largest absolute value of A
 * 
 * N: the dimention of the given matrix A (NxN)
 * A: a pointer to the given real symmetric matrix
 * i_p: a pointer to the output index i of the pivot A_ij
 * j_p: a pointer to the output index j of the pivot A_ij
 *  
 * returns: 0 if there is no exception and 1 elsewhere
 */
int find_ij_pivot(int N, double*** A, int* i_p, int* j_p){
    int max_i = 0, max_j = 0, i, j;
    double max_offDiag = 0;
    
    for(i = 1; i < N; i++){
        for(j = i + 1; j < N; j++){
            if(fabs((*A)[i][j]) >= max_offDiag){
                max_i = i;
                max_j = j;
            }
        }
    }

    if(max_i == 0 || max_j == 0){
        return 1;
    }

    *i_p = max_i;
    *j_p = max_j;

    return 0;
}

/* 
 * Function: matMult
 * -----------------
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
<<<<<<< HEAD
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
=======
    int i,j,k;
    double*** mat2_T = NULL;
    double m_ij;

    mat2_T = initMat(N);
    if (mat2_T == NULL){
        return 1;
    }

    matTranspose(N, mat2, mat2_T);

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            m_ij = 0;
            for(k = 0; k < N; k++){
                m_ij += ((*mat1)[i][k]) * ((*mat2_T)[j][k]);
                /* multipling using mat2 transposed, mat2_T for decreasing
                 the cash misses */  
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
 * ----------------
 * creates a duplication of a given matrix
 * 
 * N: the dimention of the given matrix (NxN)
 * origMat: a pointer to the given matrix
 * newMat: a pointer to the output dulicated matrix
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
>>>>>>> 42849178e70cc2287390a40ce8afd0f2834a4bb2
}

/* 
 * Function: matTranspose
 * ----------------------
 * computes a transpose matrix of a given nxn matrix
 * 
 * N: the dimention of the given matrix (NxN)
 * mat: a pointer to the given matrix
 * matT: a pointer to the output transpose matrix
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
 * ----------------------
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
 * A_tag: a pointer to the output A' matrix
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
 * -------------------------
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
int convergenceTest(int N, double epsilon, double*** mat1, double*** mat2){
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
    double** newMat = malloc(N * sizeof(double*)); // Guy - should we check if succeeded??? Liad - Yes, and if didn't we should return NULL.
    
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
