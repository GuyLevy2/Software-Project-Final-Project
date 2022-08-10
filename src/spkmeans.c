#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
    int k = 0, dimension, line_count;
    int maxIter, iteration_number;
    char* inputFile, *outputFile;
    double** vectorsList, **cluster, **centroidsList;;
    double EPSILON = 0.001;
    int inputError, writeSuccess, freeSuccess;
    double*** clusterList;
    int* index_to_insert;
    int i, j;
    double* vec, *centroid_i;
    double deltaMiu;

    inputError = validateAndProcessInput(argc, argv, &k, &dimension, &line_count, &maxIter, &inputFile, &outputFile, &vectorsList);
    
    if (!inputError){
        printf("Invalid Input!");
        return 1;
    }

    // Guy - remember to free!!???
}


/* Main Utility Functions */

/* 
 * Function: wam_func
 * ------------------
 * calculates the Weighted Adjacency matrix
 * 
 * vectors_list: the list of all N vectors
 * N: the number of vectors
 * dim: the dimension of the vectors
 * outputWamMat: a pointer to the output matrix
 * 
 * returns: 0
 */
int wam_func(double*** vectors_list, int N, int dim, double*** outputWamMat){
    int i, j;
    double **vec1, **vec2;
    double dist;
    
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            
            if (i == j){
                (*outputWamMat)[i][j] = 0;
            }
            
            else{
                vec1 = (*vectors_list)[i];
                vec2 = (*vectors_list)[j];
                dist = vectorDist(vec1, vec2, dim);
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
 * wamMat: a pointer to the WAM
 * outputDdgMat: a pointer to the ddg output matrix (initialized with zeros)
 * 
 * returns: 0
 */
int ddg_func(int N, double*** wamMat, double*** outputDdgMat){
    int i, j;
    double sum;

    for (i = 0; i < N; i++){
        sum = 0;

        for (j = 0; j < N; j++){
            sum += (*wamMat)[i][j];
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
 * wamMat: a pointer to the WAM
 * ddgMat: a pointer to the ddg matrix
 * lNormMat: a pointer to the output lNorm matrix
 * 
 * returns: NULL
 */
int lNorm_func(int N, double*** wamMat, double*** ddgMat, double*** outputLNormMat){
    int i, j;
    double ***minusSqrtMat, ***tmpMat;

    /* creating the D^-0.5 matrix */
    minusSqrtMat = initMat(N);
    if (minusSqrtMat == NULL){
        return 1;
    }
    matDup(N, ddgMat, minusSqrtMat);
    minusSqrtD(N, minusSqrtMat);
    
    /* calculating (D^-0.5*W*D^-0.5) */
    tmpMat = initMat(N);
    if (tmpMat == NULL){
        return 1;
    }
    matMult(N, minusSqrtMat, wamMat, tmpMat);
    matMult(N, tmpMat, minusSqrtMat, outputLNormMat);
    
    /* freeing utility matrices */
    freeMat(N, tmpMat);
    freeMat(N, minusSqrtMat);

    /* subtraction from identity matrix */
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                (*outputLNormMat)[i][j] = 1 - (*outputLNormMat)[i][j];
            }
            else{
                (*outputLNormMat)[i][j] = - (*outputLNormMat)[i][j];
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
 * eigenVectors: a pointer to the output eigenvectors that are arranged as a two
 *      dimensional array (initiailzed with zeroes)
 * eigenValues: a pointer to the output eigenvalues that are arranged as a one 
 *      dimensional array (initialized with zeroes)
 * 
 * returns: 0 if there is no exception and 1 otherwise
 */
int jacobi_func(int N, double*** symMat, double*** eigenVectors, double** eigenValues){
    int MAX_ROTATIONS = 100;
    int i = 0, j = 0, countRot = 0, k;
    double c = 0, s = 0;
    double ***P = NULL, ***A = NULL, ***A_tag = NULL, ***tempMat = NULL;
    double EPSILON = 0.00001;

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

    identityMat(N, eigenVectors);  /* eigenVectors = I(N) */

    matDup(N, symMat, A_tag);   /* A' = symMat (by values) */

    do{
        matDup(N, A_tag, A);    /* A = A' (by values)    */
         
        if (find_ij_pivot(N, A, &i, &j)){ /* if A is a diagonal matrix */    
            break; 
        }

        identityMat(N, P);                      /* P = I(N)                                             */
        buildRotMat(N, A, i, j, &c, &s, P);     /* P is rotate matrix wrt A and c and s are updeted     */
        matMult(N, eigenVectors, P, tempMat);   /* tempMat = eigenVectors * P                           */
        matDup(N, tempMat, eigenVectors);       /* eigenVectors = tempMat                               */
        computeA_tag(N, i, j, A, c, s, A_tag);  /* A' = P^T*A*P                                         */

        countRot++;
    } while ((!convergenceTest(N, EPSILON, A, A_tag)) 
                && (countRot < MAX_ROTATIONS));

    /* Extract the eigenValues from the diagonal of A' */
    for(k = 0; k < N; k++){
        (*eigenValues)[k] = (*A_tag)[k][k];
    }
    
    /* Free all matrixes which allocated during this function */
    freeMat(N, P);
    freeMat(N, A);
    freeMat(N, A_tag);
    freeMat(N, tempMat);
    
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
 * returns: argmax_i(delta_i), i in [1,...,floor(n/2)]
 */
int eigenGap(int N, double** eigenValues){
    int k = -1;
    int i;
    double delta;
    double maxDelta = -1;

    qsort(*eigenValues, N, sizeof(double), eigenComp);

    for (i = 0; i < floor(N/2); i++){
        delta = fabs((*eigenValues)[i+1] - (*eigenValues)[i]);
        
        if (delta > maxDelta){
            maxDelta = delta;
            k = i;
        }
    }
    
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
    return *((double*)b) - *((double*)a);
}

/* 
 * Function: vectorDist
 * --------------------
 * calculates the euclid norm of the subtraction of 2 vectors
 * 
 * v1: the first vector
 * v2: the second vector
 * dim: the dimension of the vectors
 * 
 * returns: the norm
 */
double vectorDist(double* v1, double* v2, int dim){
    double dist = 0;
    int i = 0;
    
    for(i = 0; i < dim; ++i){
        dist += pow(v1[i]-v2[i], 2.0);
    }

    return pow(dist, 0.5);
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
 * mat: a pointer to the output identity matrix
 * 
 * returns: 0 if there is no exception and 1 otherwise
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
 * N: the dimention of both given matrices A,P (NxN)
 * A: a pointer to the given symmetric matrix
 * i,j: the indexes of the pivot element A_ij of the matrix A
 * c_p, s_p: pointers to the output values of c and s
 * P: a pointer to the output rotation matrix (initializes as identity matrix)
 * 
 * returns: 0 if there is no exception and 1 otherwise
 */
int buildRotMat(int N, double*** A, int i, int j, int* c_p, int* s_p, double*** P){
    int sign_theta;
    double theta, abs_theta, t, c, s;

    theta = (((*A)[j][j]) - ((*A)[i][i])) / (2*((*A)[i][j]));
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
 * returns: 0 if A is not diagonal and 1 otherwise
 */
int find_ij_pivot(int N, double*** A, int* i_p, int* j_p){
    int max_i = 0, max_j = 0, i, j;
    double max_offDiag = 0;
    
    for(i = 1; i < N; i++){
        for(j = i + 1; j < N; j++){
            if(fabs((*A)[i][j]) > max_offDiag){
                max_i = i;
                max_j = j;
                max_offDiag = fabs((*A)[i][j]);
            }
        }
    }

    if(max_i == 0 || max_j == 0){ /* if all of the off-diagonal elements are zeros */
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

    if(freeMat(N, mat2_T)){
        return 1;
    }
    return 0; // Guy - I think there is no need for transpose and it complicates the function
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
 * returns: 0 if there is no exception and 1 otherwise
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
 * ----------------------
 * computes a transpose matrix of a given nxn matrix
 * 
 * N: the dimention of the given matrix (NxN)
 * mat: a pointer to the given matrix
 * matT: a pointer to the output transpose matrix
 * 
 * returns: 0 if there is no exception and 1 otherwise
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
 * returns: 0 if there is no exception and 1 otherwise
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
                off += pow((*mat)[i][j], 2.0);
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
 * returns: a pointer to the new matrix
 */
double*** initMat(int N){
    int i, j;
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

/* 
 * Function: freeMat
 * ---------------------
 * frees memory of a given matrix
 *  
 * N: the dimention of the given matrix (NxN)
 * mat: the matrix to be freed
 * 
 * returns: NULL 
 * 
 * Guy - there is no point in returning anything - we cannot check for errors in the functions as well
 */
int freeMat(int N, double*** mat){
    int i;

    for (i = 0; i < N; i++){
        free((*mat)[i]);
    }
    
    free(*mat);
}


/* Functions Imported From Previous Excercises */

/* 
 * Function: validateAndProcessInput
 * ---------------------------------
 * processes the CMD input and file
 *  
 * argc: argc
 * argv: argv
 * k: a pointer to k
 * dimension: a pointer to the dimension
 * line_count: a pointer to N
 * maxIter: a pointer to the max iterations number
 * inputFile: a pointer to the input file path
 * outputfile: a pointer to the output file path
 * vectors_list: a pointer to the list of all N vectors
 * 
 * returns: 0 if the input is invalid and 1 if the input is valid
 */
int validateAndProcessInput(int argc, char* argv[], int* k, int* dimension, int* line_count, int* maxIter, char** inputFile, char** outputFile, double*** vectorsList){
    char *k_str;
    char *maxIter_str;
    int outputStrLen;
    char c;
    int after_firstline;
    int curr_dimension;
    double *vec;
    int vectors_index;
    int vectorList_index;
    int i, j;
    FILE *fp1, *fp2;
    int ch; /* Guy added because of fgetc func giving error EOF??? */

    if (argc > 5 || argc < 4){
        return 0;
    }

    if (argc == 4){
        k_str = argv[1];
        maxIter_str = "200";
        *inputFile = argv[2];
        *outputFile = argv[3];
    }

    if (argc == 5){
        k_str = argv[1];
        maxIter_str = argv[2];
        *inputFile = argv[3];
        *outputFile = argv[4];
    }

    if(!isValidInteger(k_str)){ /* checking that k_str is a valid number - may cause problems and might not be needed? */
        return 0;
    }
    
    *k = atoi(k_str);
    
    if(!isValidInteger(maxIter_str)){ /* checking that str is a valid number - is it needed? */
        return 0;
    }
    
    *maxIter = atoi(maxIter_str);

    if (*k <= 1 || *maxIter < 0){
        return 0;
    }

    /* checking that the output file is in the right format */
    outputStrLen = 0;
    while ((*outputFile)[outputStrLen] != '\0'){
        outputStrLen++;
    }

    if (outputStrLen < 5){
        return 0;
    }

    if (((*outputFile)[outputStrLen - 1] != 't' 
        || (*outputFile)[outputStrLen - 2] != 'x' 
        || (*outputFile)[outputStrLen - 3] != 't' 
        || (*outputFile)[outputStrLen - 4] != '.')

        && ((*outputFile)[outputStrLen - 1] != 'c' 
        || (*outputFile)[outputStrLen - 2] != 's' 
        || (*outputFile)[outputStrLen - 3] != 'v' 
        || (*outputFile)[outputStrLen - 4] != '.'
        )){
        
        return 0;
    }

    
    fp1 = fopen(*inputFile, "r"); /* First open */
    if (fp1 == NULL){ /* if can't open the file. */
        return 0;
    }
    
    *dimension = 1;
    *line_count = 0;
    
    after_firstline = 0;
    curr_dimension = 1;
 
    for (ch = getc(fp1); ch != EOF; ch = getc(fp1)){ /* First run - dimension and number of vectors. */
        c = (char) ch; /* Guy added because of fgetc func giving error EOF - mem leak??? */
        if (c == ','){
            if (after_firstline == 0){
                *dimension += 1;
            }
            else{
                curr_dimension += 1;
            }
        }

        if (c == '\n'){
            *line_count = *line_count + 1;
            after_firstline = 1;
            
            if(*dimension != curr_dimension && curr_dimension != 1){
                fclose(fp1); /* Guy Memory leak */
                return 0;
            }
            
            curr_dimension = 1;
        }

        /* if there is a char that is not a number or , or \n. */
        if ((c != '.') && (c != '-') && (c != ',') && (c != '\n') && (isdigit(c) == 0)){
            fclose(fp1); /* `Guy Memory leak problem */
            return 0;
        }
    }
    
    fclose(fp1);
    
    if (*k >= *line_count){ /* if K >= N */
            return 0;
    } 

    fp2 = fopen(*inputFile, "r"); /* Second open */
    if (fp2 == NULL){ /* if can't open the file. */
        return 0;
    }

    *vectorsList = malloc(*line_count * sizeof(double*));
    vec = malloc(*dimension * sizeof(double));
    vectors_index = 0;
    vectorList_index = 0;
    
    for(i = 0; i < *line_count; i++){
        for (j = 0; j < *dimension - 1; j++){
            fscanf(fp2, "%lf,", &(vec[vectors_index]));
            vectors_index++;
        }

        fscanf(fp2, "%lf\n", &(vec[vectors_index]));
        vectors_index++;

        (*vectorsList)[vectorList_index] = vec;
        vectorList_index++;
        
        vectors_index = 0;
        vec = malloc(*dimension * sizeof(double));
    }

    free(vec); /* Guy Mem Testing */

    fclose(fp2);

    return 1;
}

/* 
 * Function: isValidInteger
 * ------------------------
 * checks if a string represents a valid int
 *  
 * str: the string to be checked
 * 
 * returns: 1 if the string is a valid int, and 0 otherwise
 */
int isValidInteger(char* str){
    int i;
    char c;

    if (str[0] == '\0'){
        return 0;
    }
    
    i = 0;
    c = str[0];
    while (c != '\0'){ /* checking that str is a valid number */
      if (c < 48 || c > 57){
         return 0; /* str is not an integer */
      }
      else{
         i++;
      }
    c = str[i];
   }

   return 1;
}