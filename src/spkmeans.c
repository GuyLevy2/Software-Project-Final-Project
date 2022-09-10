#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"

/* === Function Declerations === */

int eigenComp(const void*, const void*);
int minusSqrtD(int, double***);
int identityMat(int, double***);
int buildRotMat(double***, int, int, double*, double*, double***);
int find_ij_pivot(int, double***, int*, int*);
int matMult(int, double***, double***, double***);
int matDup(int, double***, double***);
int computeA_tag(int, int, int, double***, double, double, double***);
int convergenceTest(int, double, double***, double***);
double offCalc(int, double***);
void printMat(int, int, double***);
void SwitchColumnsOfMat(int, int, int, double***);
void printError();
int validateAndProcessInput(int, char**, int*, int*, char**, double***, char**);
double calcDistSquared(double*, double*, int);
int assignVectorToClosestCluster(double*, double****, double**, int**, int, int);
void resetClusters(int**, int);
double updateCentroids(int, int, double***, double***, int*);
int freeMemory(double****, int**, int);

/* === End Of Function Declerations === */

int main(int argc, char *argv[]) {
    int dimension, line_count, i, inputError, success;
    char *inputFile, *goal;
    double *eigenValues;
    double **vectorsList, **wamMat, **ddgMat, **LNormMat, **eigenVectors;

    inputError = validateAndProcessInput(argc, argv, &dimension, &line_count,
         &inputFile, &vectorsList, &goal);
    
    if (!inputError){
        printf("Invalid Input!");
        return 1;
    }

    /* in case goal == jacobi:
     computes eigenvectors and eigenvalues     
     free memory and abort.
     else - continue with the next goals */
    if (strcmp(goal, "jacobi") == 0){
        eigenVectors = initMat(line_count);
        if (eigenVectors == NULL){
            printError();
            return 1;
        }

        eigenValues = (double*)malloc(line_count * sizeof(double));
        if (eigenValues == NULL){
            printError();
            return 1;
        }

        /* Jacobi procedure call */
        success = jacobi_func(line_count, &vectorsList, &eigenVectors, &eigenValues);
        if (success == 1){
            printError();
            return 1;
        }

        /* printing output */
        for (i = 0; i < line_count - 1; i++){
            printf("%.4f", eigenValues[i]);
            printf(",");
        }
        printf("%.4f", eigenValues[line_count - 1]);
        printf("\n");
        printMat(line_count, line_count, &eigenVectors);
        
        /* free */
        free(eigenValues);
        freeMat(line_count, &eigenVectors);
        freeMat(line_count, &vectorsList);

        return 0;
    }

    /* beginning of algorithm */
    wamMat = initMat(line_count);
    if (wamMat == NULL){
        printError();
        return 1;
    }

    /* wam procedure call */
    success = wam_func(&vectorsList, line_count, dimension, &wamMat);
    if (success == 1){
        printError();
        return 1;
    }
    

    /* in case goal == wam:
     print the wam matrix:
     free memory and abort.
     else - continue with the next goal */
    if (strcmp(goal, "wam") == 0){
        printMat(line_count, line_count, &wamMat);
        
        freeMat(line_count, &wamMat);
        freeMat(line_count, &vectorsList);

        return 0;
    }

    ddgMat = initMat(line_count);
    if (ddgMat == NULL){
        printError();
        return 1;
    }

    /* ddg procedure */
    success = ddg_func(line_count, &wamMat, &ddgMat);
    if (success == 1){
        printError();
        return 1;
    }

    /* in case goal == ddg:
     print the ddg matrix:
     free memory and abort.
     else - continue with the next goal */
    if (strcmp(goal, "ddg") == 0){
        printMat(line_count, line_count, &ddgMat);

        freeMat(line_count, &wamMat);
        freeMat(line_count, &ddgMat);
        freeMat(line_count, &vectorsList);

        return 0;
    }

    LNormMat = initMat(line_count);
    if (LNormMat == NULL){
        printError();
        return 1;
    }

    /* lNorm procedure */
    success = lNorm_func(line_count, &wamMat, &ddgMat, &LNormMat);
    if (success == 1){
        printError();
        return 1;
    }
    
    /* in case goal == lnorm:
     print the lnorm matrix:
     free memory and abort.
     else - continue with the next goal */
    if (strcmp(goal, "lnorm") == 0){
        printMat(line_count, line_count, &LNormMat);
        
        freeMat(line_count, &wamMat);
        freeMat(line_count, &ddgMat);
        freeMat(line_count, &LNormMat);
        freeMat(line_count, &vectorsList);

        return 0;
    }

    return 0;
}


/* === Function Definitions === */

/* Main Utility Functions */

/* 
 * Function: wam_func
 * ------------------
 * calculates the Weighted Adjacency matrix
 * 
 * vectors_list: the list of all N vectors
 * N: the number of vectors
 * dim: the dimension of the vectors
 * outputWamMat: a pointer to the output matrix WAM
 * 
 * returns: 0
 */
int wam_func(double*** vectors_list, int N, int dim, double*** outputWamMat){
    int i, j;
    double *vec1, *vec2;
    double dist;
    
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){ 
            if (i == j){
                (*outputWamMat)[i][j] = 0;
            }
            else{
                vec1 = (*vectors_list)[i];
                vec2 = (*vectors_list)[j];
                dist = sqrt(calcDistSquared(vec1, vec2, dim));
                (*outputWamMat)[i][j] = exp(-(dist/2.0));
            }
        }
    }

    return 0;
}

/* 
 * Function: ddg_func
 * ------------------
 * calculates the Diagonal Degree matrix
 * 
 * N: the number of vectors
 * wamMat: a pointer to the Weighted Adjacency matrix
 * outputDdgMat: a pointer to the ddg output matrix (initialized with zeros)
 * 
 * returns: 0
 */
int ddg_func(int N, double*** wamMat, double*** outputDdgMat){
    int i, z;
    double sum;

    for (i = 0; i < N; i++){
        sum = 0;

        for (z = 0; z < N; z++){
            sum += (*wamMat)[i][z];
        }
        
        (*outputDdgMat)[i][i] = sum;
    }

    return 0;
}

/* 
 * Function: lNorm_func
 * --------------------
 * calculates the Normalized Graph Laplacian matrix
 * 
 * N: the number of vectors
 * wamMat: a pointer to the Weighted Adjacency matrix
 * ddgMat: a pointer to the Diagonal Degree matrix
 * lNormMat: a pointer to the output lNorm matrix
 * 
 * returns: 0
 */
int lNorm_func(int N, double*** wamMat, double*** ddgMat,
     double*** outputLNormMat){
    
    int i, j;
    double **minusSqrtMat = NULL, **tmpMat = NULL;

    /* creating the D^-0.5 matrix */
    minusSqrtMat = initMat(N);
    if (minusSqrtMat == NULL){
        return 1;
    }
    matDup(N, ddgMat, &minusSqrtMat);
    minusSqrtD(N, &minusSqrtMat);
    
    /* calculating (D^-0.5*W*D^-0.5) */
    tmpMat = initMat(N);
    if (tmpMat == NULL){
        return 1;
    }

    matMult(N, &minusSqrtMat, wamMat, &tmpMat);
    matMult(N, &tmpMat, &minusSqrtMat, outputLNormMat);

    /* freeing utility matrices */
    freeMat(N, &tmpMat);
    freeMat(N, &minusSqrtMat);

    /* subtraction from identity matrix */
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                (*outputLNormMat)[i][j] = 1 - (*outputLNormMat)[i][j];
            }
            else{
                (*outputLNormMat)[i][j] = 0 - (*outputLNormMat)[i][j];
            }
        }
    }

    return 0;
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
 * eigenVectors: a pointer to the output eigenvectors that are arranged as a two
 *      dimensional array (initiailzed with zeroes)
 * eigenValues: a pointer to the output eigenvalues that are arranged as a one 
 *      dimensional array (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 otherwise
 */
int jacobi_func(int N, double*** symMat, double*** eigenVectors, 
                                                double** eigenValues){
    int MAX_ROTATIONS = 100;
    int i = 0, j = 0, countRot = 0, k;
    double c = 0, s = 0;
    double EPSILON = 0.00001;
    double **P = NULL, **A = NULL, **A_tag = NULL, **tempMat = NULL;

    /* initializing utility matrices */
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

    P = initMat(N);
    if (P == NULL){
        return 1;
    }

    identityMat(N, eigenVectors);   /* eigenVectors = I(N)      */
    matDup(N, symMat, &A_tag);       /* A' = symMat (by values)  */

    do{
        matDup(N, &A_tag, &A);        /* A = A' (by values)       */
         
        if (find_ij_pivot(N, &A, &i, &j)){ /* if A is a diagonal matrix */    
            break; 
        }

        identityMat(N, &P);                      /* P = I(N)                                             */
        buildRotMat(&A, i, j, &c, &s, &P);     /* P is rotate matrix wrt A. Also c and s are updeted   */
        matMult(N, eigenVectors, &P, &tempMat);   /* tempMat = eigenVectors * P                           */
        matDup(N, &tempMat, eigenVectors);       /* eigenVectors = tempMat                               */
        computeA_tag(N, i, j, &A, c, s, &A_tag);  /* A' = P^T*A*P                                         */

        countRot++;
    } while ((!convergenceTest(N, EPSILON, &A, &A_tag)) 
                && (countRot < MAX_ROTATIONS));

    /* Extract the eigenValues from the diagonal of A' */
    for(k = 0; k < N; k++){
        (*eigenValues)[k] = (A_tag)[k][k];
    }
    
    /* Free all matrixes which allocated during this function */
    freeMat(N, &A);
    freeMat(N, &A_tag);
    freeMat(N, &tempMat);
    freeMat(N, &P);
    
    return 0;
}

/* 
 * Function: kmeans_c
 * --------------------
 * preforms the k-means algorithm
 * 
 * k: the number of centroids
 * dimension: the dimension of the vectors
 * line_count: the number of vectors
 * maxIter: the maximum number of iterations
 * EPSILON: the value for convergence
 * vectorsList: the list of vectors
 * centroidsList: a pointer to a list to be used for centroids
 * 
 * returns: 0 in case of correct run and 1 otherwise
 */
int kmeans_c(int k, int dimension, int line_count, int maxIter, double EPSILON,
     double** vectorsList, double*** centroidsList){
    
    double ***clusterList;
    double **cluster;
    double *vec, *centroid_i;
    double deltaMiu;
    int *index_to_insert;
    int i, j, iteration_number;
    
    index_to_insert = (int*)malloc(k * sizeof(int));
    if (index_to_insert == NULL){
        return 1;
    }

    clusterList = (double***)malloc(k * sizeof(double**));
    if (clusterList == NULL){
        return 1;
    }
    
    /* initializing cluster list, each cluster, and centroids list 
        with first k vectors - the vectors of kmeans++ algorithem */
    for (i = 0; i < k; i++){ 
        cluster = (double**)malloc(line_count * sizeof(double*));
        if (cluster == NULL){
            return 1;
        }

        clusterList[i] = cluster;
        
        vec = vectorsList[i];
        
        centroid_i = (double*)malloc(dimension * sizeof(double));
        if (centroid_i == NULL){
            return 1;
        }

        for (j = 0; j < dimension; j++){
            centroid_i[j] = vec[j];
        }
        
        (*centroidsList)[i] = centroid_i;
        index_to_insert[i] = 0;
    }
    
    deltaMiu = EPSILON*2;
    iteration_number = 0;

    /* iterations of kmeans */
    while(deltaMiu >= EPSILON && iteration_number < maxIter){
        for(i = 0; i < line_count; ++i){
            vec = vectorsList[i];
            assignVectorToClosestCluster(vec, &clusterList, (*centroidsList), &index_to_insert, k, dimension);
        }

        deltaMiu = updateCentroids(dimension, k, centroidsList, clusterList, index_to_insert);
        resetClusters(&index_to_insert, k);
        iteration_number++;
    }

    /* Free memory that allocated */
    freeMemory(&clusterList, &index_to_insert, k);

    return 0;
}


/* Granular Utility Functions */

/* 
 * Function: eigenGap
 * ------------------
 * calculates the Eigengap Heuristic measure
 * 
 * vectors_list: the list of all N vectors
 * N: the number of eigenvalues
 * eigenValues: a pointer to an array containing the eigenvalues.
 * 
 * returns: argmax_i(delta_i), i in [1,...,floor(n/2)] and -1 in case of error
 */
int eigenGap(int N, double** eigenValues){
    int k = -1;
    int i;
    double delta;
    double maxDelta = -1;
    double *eigenValuesDup = NULL;

    /* Creating a copy of eigenValues - voinding currapt the given eigenValues */
    eigenValuesDup = (double*)malloc(N * sizeof(double));
    if (eigenValuesDup == NULL){
        return -1;
    }
    for (i = 0; i < N; i++){
        eigenValuesDup[i] = (*eigenValues)[i];
    }

    qsort(eigenValuesDup, N, sizeof(double), eigenComp);

    for (i = 0; i < floor(N/2); i++){
        delta = fabs(eigenValuesDup[i+1] - eigenValuesDup[i]);
        
        if (delta > maxDelta){
            maxDelta = delta;
            k = i + 1; /* we added 1 to comply with project requirements */
        }
    }
    
    free(eigenValuesDup);
    
    return k;
}

/* 
 * Function: eigenComp
 * ------------------
 * comparator for double values for decreasing ordered qsort
 *
 * a: the first double
 * b: the second double
 * 
 * returns: >0 if (a<b), 0 if (a==b), <0 otherwise
 */
int eigenComp(const void* a, const void* b){
    double diff = ((*(double*)b) - (*(double*)a));

    if (diff > 0){
        return 1;
    }
    else if (diff < 0){
        return -1;
    }

    return 0;
}

/* 
 * Function: minusSqrtD
 * --------------------
 * processes some D matrix to its minus sqrt form
 * 
 * N: the dimention of the given matrices (NxN)
 * D: a pointer to the matrix
 * 
 * returns: 0
 */
int minusSqrtD(int N, double*** D){
    int i;
    
    for (i = 0; i < N; i++){
        (*D)[i][i] = 1 / (sqrt((*D)[i][i]));
    }

    return 0;
}

/* 
 * Function: identityMat
 * ------------------
 * sets a given matrix to be an identity matrix
 * 
 * N: the dimention of the given matrix mat (NxN)
 * mat: a pointer to the given and output matrix with some values
 * 
 * returns: 0
 */
int identityMat(int N, double*** mat){
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
 * A: a pointer to the given symmetric matrix
 * i,j: the indexes of the pivot element A_ij of the matrix A
 * c_p, s_p: pointers to the output values of c and s
 * P: a pointer to the output rotation matrix (initializes as identity matrix)
 * 
 * returns: 0 if there is no exception and 1 otherwise
 */
int buildRotMat(double*** A, int i, int j, double* c_p,
     double* s_p, double*** P){
    
    int sign_theta;
    double theta, abs_theta, t, c, s;

    theta = (((*A)[j][j]) - ((*A)[i][i])) / (2*((*A)[i][j]));
    abs_theta = fabs(theta);
    sign_theta = theta >= 0 ? 1 : -1;

    t = (sign_theta*1.0) / (abs_theta + sqrt((theta * theta) + 1));

    c = 1.0 / sqrt((t * t) + 1);
    s = t*c;

    *c_p = c;
    *s_p = s;
    (*P)[i][i] = c;
    (*P)[i][j] = s;
    (*P)[j][j] = c;
    (*P)[j][i] = -1.0*s;

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
 * returns: 0 if A is not diagonal and 1 otherwise
 */
int find_ij_pivot(int N, double*** A, int* i_p, int* j_p){
    int max_i = 0, max_j = 0, i, j;
    double max_offDiag = 0;
    
    for(i = 0; i < N; i++){
        for(j = i + 1; j < N; j++){
            if(fabs((*A)[i][j]) > max_offDiag){
                max_i = i;
                max_j = j;
                max_offDiag = fabs((*A)[i][j]);
            }
        }
    }

    if(max_i == 0 && max_j == 0){ /* if all of the off-diagonal elements are zeros */
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
 * returns: 0 if there is no exception and 1 otherwise
 */
int matMult(int N, double*** mat1, double*** mat2, double*** outputMat){
    int i,j,k;
    double m_ij;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            m_ij = 0;
            for(k = 0; k < N; k++){
                m_ij += ((*mat1)[i][k]) * ((*mat2)[k][j]);
            }
            (*outputMat)[i][j] = m_ij;
        }
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
 * returns: 0
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
 * returns: 0
 */
int computeA_tag(int N, int i, int j, double*** A, double c,
     double s, double*** A_tag){
    
    int r;
    double a_tag_ri, a_tag_rj;
    double s_squared, c_squared;
    double a_ii, a_jj, a_ij;

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

    s_squared = (s * s);
    c_squared = (c * c);
    a_ii = (*A)[i][i];
    a_jj = (*A)[j][j];
    a_ij = (*A)[i][j];

    (*A_tag)[i][i] = c_squared*a_ii + s_squared*a_jj - 2*s*c*a_ij;
    (*A_tag)[j][j] = s_squared*a_ii + c_squared*a_jj + 2*s*c*a_ij;
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
 * -----------------
 * computes the sum squares of all off-diagonal elements of a given matrix
 *  
 * N: the dimention of the given matrix (NxN)
 * mat: a pointer to the matrix
 * 
 * returns: off(mat)^2 as described above
 */
double offCalc(int N, double*** mat){
    double off = 0;
    int i,j;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            if(j != i){
                off += ((*mat)[i][j] * (*mat)[i][j]);
            }
        }
    }

    return off;
}

/* 
 * Function: initMat
 * ---------------------
 * initizlises a zero - matrix in dimension N
 *  
 * N: the dimention of the given matrix (NxN)
 * 
 * returns: a new initializes matrix or NULL if an error has occurred
 */
double** initMat(int N){
    int i, j;
    double *vec = NULL;
    double **newMat = NULL;
    
    newMat = (double**)malloc(N * sizeof(double*)); 
    if(newMat == NULL){
        return NULL;
    }

    for (i = 0; i < N; i++){
        vec = (double*)malloc(N * sizeof(double));
        if (vec == NULL){
            return NULL;
        }

        for (j = 0; j < N; j++){
            vec[j] = 0;
        }

        newMat[i] = vec;
    }

    return newMat;
}

/* 
 * Function: initMatMN
 * ---------------------
 * initizlises a zero - matrix in dimension MxN
 *  
 * N: the dimention of the given matrix (MxN)
 * M: the dimention of the given matrix (MxN)
 * 
 * returns: a new initializes matrix or NULL if an error has occurred
 */
double** initMatMN(int rows, int cols){
    int i, j;
    double *vec = NULL;
    double **newMat = NULL;
    
    newMat = (double**)malloc(rows * sizeof(double*));
    if(newMat == NULL){
        return NULL;
    }

    for (i = 0; i < rows; i++){
        vec = (double*)malloc(cols * sizeof(double));
        if (vec == NULL){
            return NULL;
        }

        for (j = 0; j < cols; j++){
            vec[j] = 0;
        }

        newMat[i] = vec;
    }

    return newMat;
}

/* 
 * Function: freeMat
 * ---------------------
 * frees memory of a given matrix
 *  
 * N: the number of rows of the given matrix (Nx*)
 * mat: the matrix to be freed
 * 
 * returns: - 
 * 
 */
void freeMat(int N, double*** mat){
    int i;
    for (i = 0; i < N; i++){
        free((*mat)[i]);
    }
    free(*mat);
}

/* 
 * Function: printMat
 * ---------------------
 * prints the input matrix
 *  
 * m: the number of rows of the given matrix
 * n: the number of columns of the given matrix
 * mat: the matrix to be printed
 * 
 * returns: NULL 
 * 
 */
void printMat(int m, int n, double*** mat){
    int i,j;

    for(i = 0; i < m; ++i){
        for (j = 0; j < n - 1; ++j){
            printf("%.4f", (*mat)[i][j]);
            printf(",");
        }
        printf("%.4f", (*mat)[i][n-1]);
        printf("\n");
    }
}

/* 
 * Function: printError
 * ---------------------
 * prints the error message
 * 
 * returns: NULL 
 * 
 */
void printError(){
    printf("An Error Has Occurred");
}

/* Functions for spk proccedure */

/* 
 * Function: sortEigenValuesAndEigenVectors
 * ---------------------------------
 * sorting the eigenvalues and corresponding eigenvectors in decreasing order
 * 
 * N: the number of eigenvalues and eigenvectors and the dimension of eigenvectors
 * eigenValues: a pointer to the list of eigenvalues
 * eigenVectos: a pointer to the list of eigenvectors
 * 
 * returns: 0
 */
int sortEigenValuesAndEigenVectors(int N, double **eigenValues,
     double ***eigenVectors){
    
    int i;
    double temp;
    int flag = 1;

    while (flag)
    {
        flag = 0;
        for (i = 0; i < N-1; i++){
            if ((*eigenValues)[i] < (*eigenValues)[i+1]){
                /* Switch operation */
                flag = 1;
                
                temp = (*eigenValues)[i];
                (*eigenValues)[i] = (*eigenValues)[i+1];
                (*eigenValues)[i+1] = temp;

                SwitchColumnsOfMat(N, i, i+1, eigenVectors);
            }
        }
    }
    return 0;
}

/* 
 * Function: SwitchColumnsOfMat
 * ---------------------------------
 * switching the columns of a given matrix
 * 
 * numOfRows: the number of rows in the matrix
 * i,j: the index of columns to be switched
 * mat: a pointer to the matrix
 * 
 * returns: -
 */
void SwitchColumnsOfMat(int numOfRows, int i, int j, double ***mat){
    double temp;
    int m;

    for (m = 0; m < numOfRows; m++){
        temp = (*mat)[m][i];
        (*mat)[m][i] = (*mat)[m][j];
        (*mat)[m][j] = temp;
    }
}

/* 
 * Function: Fill_K_LargestEigenVectors
 * ---------------------------------
 * initialize matrix U with the values in the first k columns of given matrix
 * 
 * N: dimensions of the matrices (NxN)
 * K: the number of columns to copy
 * eigenVectors: a pointer to the matrix from which the columns will be copied
 * U: a pointer to the matrix to copy into
 * 
 * returns: 0
 */
int Fill_K_LargestEigenVectors(int N, int K, double ***eigenVectors, double ***U){
    int i, j;

    for (i = 0; i < N; i++){
        for (j = 0; j < K; j++){
            (*U)[i][j] = (*eigenVectors)[i][j];
        }
    }

    return 0;
}

/* 
 * Function: ReNormalizedRows
 * ---------------------------------
 * re-normalizes the rows of a given matrix
 * 
 * N: number of rows in the matrices
 * K: number of columns in the matrices
 * U: input matrix
 * T: output matrix
 * 
 * returns: 0
 */
int ReNormalizedRows(int N, int K, double ***U, double ***T){
    int i, j;
    double rowSum = 0;

    for (i = 0; i < N; i++){
        rowSum = 0;
        for (j = 0; j < K; j++){
            rowSum += ((*U)[i][j] * (*U)[i][j]);
        } 
        rowSum = sqrt(rowSum);

        if (rowSum != 0){
            for (j = 0; j < K; j++){
                (*T)[i][j] = (*U)[i][j] / rowSum;
            }
        }
    }

    return 0;
}

/* 
 * Function: validateAndProcessInput
 * ---------------------------------
 * processes the CMD input and file
 *  
 * argc: argc
 * argv: argv
 * dimension: a pointer to the dimension
 * line_count: a pointer to N
 * inputFile: a pointer to the input file path
 * vectors_list: a pointer to the list of all N vectors
 * goal: a pointer to the goal string
 * 
 * returns: 0 if the input is invalid and 1 if the input is valid
 */
int validateAndProcessInput(int argc, char* argv[], int* dimension,
     int* line_count, char** inputFile, double*** vectorsList, char** goal){

    int after_firstline;
    double *vec;
    int vectors_index;
    int vectorList_index;
    FILE *fp1, *fp2;
    char *strNum;
    char line[2048];


    if (argc != 3){
        return 0;
    }
    
    *goal = argv[1];
    if (strcmp(*goal, "wam") != 0 &&
        strcmp(*goal, "ddg") != 0 &&
        strcmp(*goal, "lnorm") != 0 &&
        strcmp(*goal, "jacobi") != 0){
    
        return 0;
    }


    *inputFile = argv[2];
    
    /* First open - for dimension and number of vectors */
    fp1 = fopen(*inputFile, "r");
    if (fp1 == NULL){ /* if can't open the file. */
        return 0;
    }
    
    *dimension = 0;
    *line_count = 0;
    after_firstline = 0;
 
    while (fgets(line, 2048, fp1)){
        strNum = strtok(line, ",");

        while (strNum){
            if (after_firstline == 0){
                *dimension += 1;
            }

            strNum = strtok(NULL, ",");
        }
        
        after_firstline = 1;
        *line_count += 1;
    }
    
    fclose(fp1);

    /* Second open - reading the values and allocating vectors in memory */
    fp2 = fopen(*inputFile, "r"); 
    if (fp2 == NULL){ /* if can't open the file. */
        return 0;
    }

    *vectorsList = malloc(*line_count * sizeof(double*));
    vectors_index = 0;
    vectorList_index = 0;
    
    while (fgets(line, 2048, fp2)){
        vec = malloc(*dimension * sizeof(double));
        strNum = strtok(line, ",");

        while (strNum){
            vec[vectors_index] = atof(strNum);
            vectors_index++;
            strNum = strtok(NULL, ",");
        }

        (*vectorsList)[vectorList_index] = vec;
        vectorList_index++;
        vectors_index = 0;
    }

    fclose(fp2);

    return 1;
}

/* 
 * Function: calcDistSquared
 * ------------------------
 * computes the Euclidian distance between to vectors, squared
 *  
 * v1: the first vector
 * v2: the second vector
 * d: the dimension of the vectors
 * 
 * returns: the Euclidian distance between to vectors, squared
 */
double calcDistSquared(double* v1, double* v2, int d){
    double dist = 0;
    int i = 0;
    
    for(i = 0; i < d; ++i){
        dist += ((v1[i]-v2[i]) * (v1[i]-v2[i]));
    }
    return dist;
}

/* 
 * Function: assignVectorToClosestCluster
 * ------------------------
 * finds and assignes a given vector to the closest cluster by its centroid
 *  
 * 
 * vector: the vector to assign
 * clusterList: a pointer to the list of clusters
 * centroidsList: the list of centroids
 * index_to_insert: a pointer to the list of last insert indexes for each cluster
 * k: the number of centroids
 * d: the dimension of the vectors
 * 
 * returns: 1
 */
int assignVectorToClosestCluster(double* vector, double**** clusterList,
     double** centroidsList, int** index_to_insert, int k, int d){
    
    double minDist = -1.0, dist = 0.0;
    int minDistIndex = -1, i = 0;
    double* cent_vec;
    int last_index;
    
    for (i = 0; i < k; ++i){
        cent_vec = centroidsList[i];
        dist = calcDistSquared(vector, cent_vec, d);
        
        if ((dist < minDist) || (minDist == -1.0)){
            minDist = dist;
            minDistIndex = i;
        }
    }
    
    /* insert vector to the closest cluster list at the last index to insert 
        of the list and update last index */
    last_index = (*index_to_insert)[minDistIndex];
    (*clusterList)[minDistIndex][last_index] = vector;
    (*index_to_insert)[minDistIndex] = last_index + 1;

    return 1;   
}

/* 
 * Function: resetClusters
 * ------------------------
 * resets the last indexes for insert in the clusters 
 *  
 * index_to_insert: a pointer to the list of last insert indexes for each cluster
 * k: the number of centroids
 * 
 * returns: -
 */
void resetClusters(int** index_to_insert, int k){
   int i;
   
   for(i = 0; i < k; i++){
       (*index_to_insert)[i] = 0;
   }
}

/* 
 * Function: updateCentroids
 * ------------------------
 * updates the centroids for each cluster
 * 
 * dimension: the dimension of the vectors 
 * k: the number of centroids
 * centroidsList: a pointer to the list of centroids 
 * clusterList: the list of clusters
 * index_to_insert: a pointer to the list of last insert indexes for each cluster
 * 
 * returns: the maximum distance between an old and a new centroid for some centroid
 */
double updateCentroids(int dimension, int k, double*** centroidsList,
     double*** clusterList, int* index_to_insert){
    
    double maxDelta = -1.0, centroidDelta = 0.0, indexDelta = 0.0;
    int i, j, m;

    double indexAverage = 0.0, indexSum = 0.0;
    int lastIndex = 0;

    for (i = 0; i < k; i++){ /* for every cluster */
        centroidDelta = 0.0;
        
        for (j = 0; j < dimension; j++){ /* for every index in dimension */
            indexAverage = 0.0;
            lastIndex = index_to_insert[i];
            indexSum = 0.0;

            for (m = 0; m < lastIndex; m++){ /* for every vector */
                indexSum += clusterList[i][m][j];
            }

            indexAverage = indexSum / lastIndex;
            indexDelta = indexAverage - (*centroidsList)[i][j];
            (*centroidsList)[i][j] = indexAverage;
            centroidDelta += indexDelta * indexDelta;
        }
        
        centroidDelta = sqrt(centroidDelta);

        if (centroidDelta > maxDelta || maxDelta == -1.0){
            maxDelta = centroidDelta;
        }
    }

    return maxDelta;
}

/* 
 * Function: freeMemory
 * ------------------------
 * frees memory for clusters, and index to insert list
 * 
 * clusterList: a pointer to the list of cluster 
 * index_to_insert: a pointer to the list of last insert indexes for each cluster
 * k: the number of centroids
 * 
 * returns: 1
 */
int freeMemory(double**** clusterList, int** index_to_insert, int k){
    int i;
    /*free clusters */
    for(i = 0; i < k; i++){
        free((*clusterList)[i]);
    }

    /*free clusters list */
    free(*clusterList);

    /*free index_to_insert */
    free(*index_to_insert);

    return 1;
}

/* === End Of Function Definitions === */