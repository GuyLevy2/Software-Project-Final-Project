#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

/* API functions */

/* 
 * Function: eigengapHeuristic_fit
 * -------------------------------
 * API function for Eigengap Heuristic Algorithem
 * 
 * Input: PyObject arguments:
 * (1) The number of eigen-values (int)
 * (2) A list of eigen-values (as float)
 * 
 * Output: the number k of the eigengap heuristic
 * returns: PyObject of the number k (int)
 */
static PyObject* eigengapHeuristic_fit(PyObject *self, PyObject *args){
    int N, k, i;
    double *eigenValues = NULL;
    PyObject *eigenValues_obj;
    
    /* Input */
    if(!PyArg_ParseTuple(args, "iO", &N, &eigenValues_obj)) {
        return Py_BuildValue("");
    }

    eigenValues = (double*)malloc(N * sizeof(double));
    if(eigenValues == NULL){
        return Py_BuildValue("");
    }
    
    for (i = 0; i < N; i++){
        eigenValues[i] = PyFloat_AsDouble(PyList_GetItem(eigenValues_obj, i));
    }

    /* Body */
    k = eigenGap(N, &eigenValues);

    /* Free eigenValues */
    free(eigenValues);

    /* Output */
    if(k == -1){
        return Py_BuildValue("");
    }

    return Py_BuildValue("i", k);
}

/* 
 * Function: kmeans_fit
 * --------------------
 * API function for kmeans Algorithem
 * 
 * Input: PyObject arguments:
 * (1) the number of centroids (int)
 * (2) the dimension of each vector (int)
 * (3) the number of vectors (int)
 * (4) the maximum iterations number during the algorithem
 * (5) a value of epsilon for convergence during the algorithem
 * (6) A list of vectors (as float)
 * 
 * Output: the finel k cetroids
 * returns: PyObject as list of floats
 */
static PyObject* kmeans_fit(PyObject *self, PyObject *args){
    int k, dimension, line_count, maxIter, kmeans_success; 
    int i,j;
    double EPSILON;
    double **vectorsList, **centroids_list;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */

    PyObject *vec_list_obj;
    PyObject *centroids_list_obj; /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iiiidO", &k, &dimension, &line_count, &maxIter, &EPSILON, &vec_list_obj)) {
        return Py_BuildValue("");
    }

    retInitMat = Create_C_Mat_From_PyObj(N, dimension, vec_list_obj);
    if(retInitMat == NULL){
        return Py_BuildValue("");
    }
    vectorsList = *retInitMat;
    
    /* Body */
    centroids_list = malloc(k * sizeof(double*));   /* The inner vectors are initialized in kmeans_c function !!! */
    if (centroids_list == NULL){
        return Py_BuildValue("");
    }

    kmeans_success = kmeans_c(k, dimension, line_count, maxIter, EPSILON, vectorsList, &centroids_list);
    if (kmeans_success == 0){
        return Py_BuildValue("");
    }
    
    /* Free vectorsList */
    freeMat(&vectorsList);

    /* Output */
    centroids_list_obj = Create_PyObj_Mat_From_C(k, dimension, &centroids_list);

    /* Free centroids_list */
    freeMat(k, &centroids_list);

    return Py_BuildValue("O", centroids_list_obj);
}

/* 
 * Function: spk_fit
 * -----------------
 * API function for executing the spk akgorithem steps 1-5
 * 
 * Input: PyObject arguments:
 * (1) the dimension of each vector (int)
 * (2) the number of vectors (int)
 * (3) the number of clusters k (or k=0 for eigengup heuristic method)
 * (3) A list of vectors (as float)
 * 
 * Output: the finel T matrix (step 5)
 * returns: PyObject as 2-dimensional list of floats (matrix)
 */
static PyObject* spk_fit(PyObject *self, PyObject *args){
    int N, dimension, K, success;
    double **vectorsList = NULL, **wamMat = NULL, **ddgMat = NULL;
    double **lnormMat = NULL, **eigenVectors = NULL, *eigenValues = NULL;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */
    double **U, **T;
    PyObject *vec_list_obj;     /* input */
    PyObject *T_Matrix_obj;    /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iiIO", &dimension, &N, &K &vec_list_obj)) {
        return Py_BuildValue("");
    }

    retInitMat = Create_C_Mat_From_PyObj(N, dimension, vec_list_obj);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    vectorsList = *retInitMat;

    /* Body */
    /* WAM */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    wamMat = *retInitMat;
    success = wam_func(&vectorsList, N, dimension, &wamMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* DDG */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    ddgMat = *retInitMat;
    success = ddg_func(N, &wamMat, &ddgMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* lnorm */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    lnormMat = *retInitMat;
    success = lNorm_func(N, &wamMat, &ddgMat, &lnormMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* eigenValues & eigenVectors (Jacobi algo) */
    eigenValues = (double*)malloc(N * sizeof(double));
    if (eigenValues == NULL){
        return Py_BuildValue("");
    }
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    eigenVectors = *retInitMat;
    success = jacobi_func(N, &lnormMat, &eigenVectors, &eigenValues);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* eigenGap Heuristic method (in case of k=0) */
    if (K == 0){
        K = eigenGap(N, &eigenValues);
    }

    /* U = largest_k_eigenvectors */
    success = sortEigenValuesAndEigenVectors(N, &eigenValues, &eigenVectors);
    if(success == 1){
        return Py_BuildValue("");
    }
    retInitMat = initMatMN(N, K);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    U = *retInitMat;
    success = Fill_K_LargestEigenVectors(N, K, &eigenVectors, &U);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* T = renormlized(U) */
    retInitMat = initMatMN(N, K);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    T = *retInitMat;
    success = ReNormalizedRows(N, K, &U, &T);
    if(success == 1){
        return Py_BuildValue("");
    }

    /* Free */
    freeMat(N, &vectorsList);
    freeMat(N, &wamMat);
    freeMat(N, &ddgMat);
    freeMat(N, &lnormMat);
    free(eigenValues);
    freeMat(N, &eigenVectors);
    freeMat(N, &U);

    /* Build output */
    T_Matrix_obj = Create_PyObj_Mat_From_C(N, K, &T);
    
    /* free T */
    freeMat(N, &T);

    return Py_BuildValue("O", T_Matrix_obj);
}

/* 
 * Function: wam_fit
 * -----------------
 * API function for computing the Weighted Adjacency matrix
 * 
 * Input: PyObject arguments:
 * (1) the dimension of each vector (int)
 * (2) the number of vectors (int)
 * (3) A list of vectors (as float)
 * 
 * Output: the finel WAM of the given vectors list
 * returns: PyObject as 2-dimensional list of floats (matrix)
 */
static PyObject* wam_fit(PyObject *self, PyObject *args){
    int N, dimension, success;
    double **vectorsList = NULL, **wamMat = NULL;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */

    PyObject *vec_list_obj;     /* input */
    PyObject *wamMatrix_obj;    /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iiO", &N, &dimension, &vec_list_obj)) {
        return Py_BuildValue("");
    }

    retInitMat = Create_C_Mat_From_PyObj(N, dimension, vec_list_obj);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    vectorsList = *retInitMat;

    /* Body */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    wamMat = *retInitMat;

    success = wam_func(&vectorsList, N, dimension, &wamMat);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free vectorsList */
    freeMat(N, &vectorsList);

    /* Build output */
    wamMatrix_obj = Create_PyObj_Mat_From_C(N, N, &wamMat);
    
    /* free wamMat */
    freeMat(N, &wamMat);

    return Py_BuildValue("O", wamMatrix_obj);
}

/* 
 * Function: ddg_fit
 * -----------------
 * API function for computing the Diagonal Degree matrix
 * 
 * Input: PyObject arguments:
 * (1) the dimension N of each given matrix (NxN) (int)
 * (2) A 2-dimensional list of Weighted Adjacency matrix (as float)
 * 
 * Output: the finel DDG of the given WAM
 * returns: PyObject as 2-dimensional list of floats (matrix)
 */
static PyObject* ddg_fit_B(PyObject *self, PyObject *args){ /* Option B */
    int N, success;
    double **ddgMat = NULL, **wamMat = NULL;
    PyObject *wamMat_obj;     /* input */
    PyObject *ddgMat_obj;    /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iO", &N, &wamMat_obj)) {
        return Py_BuildValue("");
    }

    wamMat = *(Create_C_Mat_From_PyObj(N, N, wamMat_obj));
    if (wamMat == NULL){
        return Py_BuildValue("");
    }

    /* Body */
    ddgMat = *(initMat(N));
    if (ddgMat == NULL){
        return Py_BuildValue("");
    }

    success = ddg_func(N, &wamMat, &ddgMat);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free wamMat */
    freeMat(N, &wamMat);

    /* Output */
    ddgMat_obj = Create_PyObj_Mat_From_C(N, N, &ddgMat);
    
    /* Free ddgMat */
    freeMat(N, &ddgMat);

    return Py_BuildValue("O", ddgMat_obj);
}
static PyObject* ddg_fit(PyObject *self, PyObject *args){ /* Option A */
    int N, dimension, success;
    double **vectorsList = NULL, **ddgMat = NULL, **wamMat = NULL;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */

    PyObject *vec_list_obj;     /* input */
    PyObject *ddgMat_obj;    /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iiO", &N, &dimension, &vec_list_obj)) {
        return Py_BuildValue("");
    }
    
    retInitMat = Create_C_Mat_From_PyObj(N, dimension, vec_list_obj);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    vectorsList = *retInitMat;

    /* Body */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    wamMat = *retInitMat;
    success = wam_func(&vectorsList, N, dimension, &wamMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    ddgMat = *retInitMat;
    success = ddg_func(N, &wamMat, &ddgMat);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free */
    freeMat(N, &vectorsList);
    freeMat(N, &wamMat);

    /* Output */
    ddgMat_obj = Create_PyObj_Mat_From_C(N, N, &ddgMat);
    
    /* Free ddgMat */
    freeMat(N, &ddgMat);

    return Py_BuildValue("O", ddgMat_obj);
}

/* 
 * Function: lnorm_fit
 * -------------------
 * API function for computing the Normalized Graph Laplacian matrix
 * 
 * Input: PyObject arguments:
 * (1) the dimension N of each given matrix (NxN) (int)
 * (2) A 2-dimensional list of Weighted Adjacency matrix (as float)
 * (3) A 2-dimensional list of Diagonal Degree matrix (as float)
 * 
 * Output: the L-norm matrix of the given WAM and DDG
 * returns: PyObject as 2-dimensional list of floats (matrix)
 */
static PyObject* lnorm_fit_B(PyObject *self, PyObject *args){ /* Option B */
    int N, success;
    double **ddgMat = NULL, **wamMat = NULL, **lnormMat = NULL;
    PyObject *wamMat_obj, *ddgMat_obj;     /* input */
    PyObject *lnormMat_obj;    /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iOO", &N, &wamMat_obj, &ddgMat_obj)) {
        return Py_BuildValue("");
    }

    wamMat = *(Create_C_Mat_From_PyObj(N, N, wamMat_obj));
    if (wamMat == NULL){
        return Py_BuildValue("");
    }

    ddgMat = *(Create_C_Mat_From_PyObj(N, N, ddgMat_obj));
    if (ddgMat == NULL){
        return Py_BuildValue("");
    }

    /* Body */
    lnormMat = *(initMat(N));
    if (lnormMat == NULL){
        return Py_BuildValue("");
    }

    success = lNorm_func(N, &wamMat, &ddgMat, &lnormMat);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free wamMat and ddgMat */
    freeMat(N, &wamMat);
    freeMat(N, &ddgMat);

    /* Output */
    lnormMat_obj = Create_PyObj_Mat_From_C(N, N, &lnormMat);
    
    /* Free lnormMat */
    freeMat(N, &lnormMat);

    return Py_BuildValue("O", lnormMat_obj);
}
static PyObject* lnorm_fit(PyObject *self, PyObject *args){ /* Option A */
    int N, dimension, success;
    double **vectorsList = NULL, **ddgMat = NULL, **wamMat = NULL;
    double **lnormMat = NULL;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */

    PyObject *vec_list_obj;     /* input */
    PyObject *lnormMat_obj;    /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iiO", &N, &dimension, &vec_list_obj)) {
        return Py_BuildValue("");
    }
    
    retInitMat = Create_C_Mat_From_PyObj(N, dimension, vec_list_obj);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    vectorsList = *retInitMat;

    /* Body */
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    wamMat = *retInitMat;
    success = wam_func(&vectorsList, N, dimension, &wamMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    ddgMat = *retInitMat;
    success = ddg_func(N, &wamMat, &ddgMat);
    if(success == 1){
        return Py_BuildValue("");
    }

    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    lnormMat = *retInitMat;
    success = lNorm_func(N, &wamMat, &ddgMat, &lnormMat);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free */
    freeMat(N, &vectorsList);
    freeMat(N, &wamMat);
    freeMat(N, &ddgMat);

    /* Output */
    lnormMat_obj = Create_PyObj_Mat_From_C(N, N, &lnormMat);
    
    /* Free lnormMat */
    freeMat(N, &lnormMat);

    return Py_BuildValue("O", lnormMat_obj);
}

/* 
 * Function: jacobi_fit
 * --------------------
 * API function for the Jacobian algorithem in order to find 
 *  eigenvalues and eigenvectors of a given symmetric real-valued matrix
 * 
 * Input: PyObject arguments:
 * (1) the dimension N of the given matrix (NxN) (int)
 * (2) A 2-dimensional list: symmetric matrix (as float)
 * 
 * Output: the eigenvalues and eigenvectors of the given matrix
 * returns: PyObject as tuple of 2 arguments:
 *          (1) A list of floats - eigenvalues
 *          (2) A 2-dimensional list of floats - eigenvectors
 */
static PyObject* jacobi_fit(PyObject *self, PyObject *args){
    int N, success;
    double **symMat = NULL;
    double **eigenVectors = NULL;
    double *eigenValues = NULL;
    double ***retInitMat = NULL; /* The return pointer to matrix which back after call to initMat */

    PyObject *symMat_obj;     /* input */
    PyObject *eigenVecors_obj, *eigenValues_obj;    /* output */

    /* Input */
    if(!PyArg_ParseTuple(args, "iO", &N, &symMat_obj)) {
        return Py_BuildValue("");
    }

    retInitMat = Create_C_Mat_From_PyObj(N, N, symMat_obj);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    symMat = *retInitMat;

    /* Body */
    eigenValues = (double*)malloc(N * sizeof(double));
    if (eigenValues == NULL){
        return Py_BuildValue("");
    }
    retInitMat = initMat(N);
    if (retInitMat == NULL){
        return Py_BuildValue("");
    }
    eigenVectors = *retInitMat;

    success = jacobi_func(N, &symMat, &eigenVectors, &eigenValues);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free symMat */
    freeMat(N, &symMat);

    /* Output */
    eigenVecors_obj = Create_PyObj_Mat_From_C(N, N, &eigenVectors);
    eigenValues_obj = Create_PyObj_Arr_From_C(N, &eigenValues);
    
    /* Free eigenVectors and eigenValues */
    freeMat(N, &eigenVectors);
    free(eigenValues);

    return Py_BuildValue("(OO)", eigenValues_obj, eigenVecors_obj); /* As Tuple of 2 objects (OO) */
}

/* Utility functions*/

/* 
 * Function: Create_C_Mat_From_PyObj
 * ---------------------------------
 * creates, allocates and duplicates a new C matrix (2-dimensional array) from
 *  corresponding PyObject (as double)
 * 
 * numOfRows: the number of rows in the given matrix
 * numOfCols: the number of columns in the given matrix
 * py_mat: PyObject of matrix
 * 
 * returns: a pointer to the new aloocated matrix (NULL in case of error)
 */
double*** Create_C_Mat_From_PyObj(int numOfRows, int numOfCols, PyObject* py_mat){
    double **mat = NULL;
    double *row = NULL;
    PyObject *row_obj, *value;
    
    mat = (double**)malloc(numOfRows * sizeof(double*));
    if (mat == NULL){
        return NULL;
    }

    for (i = 0; i < numOfRows; i++){
        row = (double*)malloc(numOfCols * sizeof(double));
        if (row == NULL){
            return NULL;
        }
        mat[i] = row;   
    }

    for (i = 0; i < numOfRows; i++) {
        row_obj = PyList_GetItem(py_mat, i);
        for (j = 0; j < numOfCols; j++){
            value = PyList_GetItem(row_obj, j);
            mat[i][j] = PyFloat_AsDouble(value);
        }
    }

    return &mat;
}

/* 
 * Function: Create_PyObj_Mat_From_C
 * ---------------------------------
 * creates and duplicates a PyObject matrix (2-dimensional array) from
 *  corresponding C matrix (as float)
 * 
 * numOfRows: the number of rows in the given matrix
 * numOfCols: the number of columns in the given matrix
 * c_mat: a pointer to 2-dimensional array of double in C
 * 
 * returns: a pointer to the new PyObject matrix
 */
PyObject* Create_PyObj_Mat_From_C(int numOfRows, int numOfCols, double*** c_mat){
    PyObject *outMatrix_obj, *matRow_obj; 
    int i, j;

    outMatrix_obj = PyList_New(numOfRows);
    for (i = 0; i < numOfRows; i++) { 
        matRow_obj = PyList_New(numOfCols);
        for (j = 0; j < numOfCols; j++){
            PyList_SetItem(matRow_obj, j, Py_BuildValue("d", (*c_mat)[i][j]));
        }
        PyList_SetItem(outMatrix_obj, i, matRow_obj);
    }

    return outMatrix_obj;
}

/* 
 * Function: Create_C_Arr_From_PyObj
 * ---------------------------------
 * creates, allocates and duplicates a new C array (1-dimensional array) from
 *  corresponding PyObject (as double)
 * 
 * numOfElements: the number of elements in the given array
 * py_arr: PyObject of array
 * 
 * returns: a pointer to the new aloocated array (NULL in case of error)
 */
double*** Create_C_Arr_From_PyObj(int numOfElements, PyObject* py_arr){
    double *arr = NULL;
    
    arr = (double*)malloc(numOfElements * sizeof(double));
    if (arr == NULL){
        return NULL;
    }

    for (i = 0; i < numOfElements; i++){
        arr[i] = PyFloat_AsDouble(List_GetItem(py_arr, i));   
    }

    return &arr;
}

/* 
 * Function: Create_PyObj_Arr_From_C
 * ---------------------------------
 * creates and duplicates a PyObject array (1-dimensional array) from
 *  corresponding 1-dimensional array in C (as float)
 * 
 * numOfElements: the number of elements in the given array
 * c_arr: a pointer to 1-dimensional array of double in C
 * 
 * returns: a pointer to the new PyObject array
 */
PyObject* Create_PyObj_Arr_From_C(int numOfElements, double** c_arr){
    PyObject *outArr_obj; 
    int i;

    outArr_obj = PyList_New(numOfElements);
    for (i = 0; i < numOfElements; i++) { 
        PyList_SetItem(outArr_obj, i, Py_BuildValue("d", (*c_arr)[i]));
    }

    return outArr_obj;
}

int sortEigenValuesAndEigenVectors(int N, double **eigenValues, double ***eigenVectors){
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

void SwitchColumnsOfMat(int numOfRows, int i, int j, double ***mat){
    double temp;
    int m;

    for (m = 0; m < numOfRows; m++){
        temp = (*mat)[m][i];
        (*mat)[m][i] = (*mat)[m][j];
        (*mat)[m][j] = temp;
    }
}

int Fill_K_LargestEigenVectors(int N, int K, double ***eigenVectors, double ***U){
    int i, j;

    for (i = 0; i < N; i++){
        for (j = 0; j < K; j++){
            (*U)[i][j] = (*eigenVectors)[i][j];
        }
    }

    return 0;
}

int ReNormalizedRows(int N, int K, double ***U, double ***T){
    int i, j;
    double rowSum = 0;

    for (i = 0; i < N; i++){
        rowSum = 0;
        for (j = 0; j < K; j++){
            rowSum += pow((*U)[i][j], 2.0);
        } 
        rowSum = pow(rowSum, 0.5);

        for (j = 0; j < K; j++){
            (*T)[i][j] = (*U)[i][j] / rowSum;
        }
    }

    return 0;
}

/* Setup Area */

static PyMethodDef capiMethods[] = {
    {"eigengapHeuristic_fit",
    (PyCFunction) eigengapHeuristic_fit,
    METH_VARARGS,
    PyDoc_STR("Eigengap Heuristic Algorithem")},
    {"kmeans_fit",
    (PyCFunction) kmeans_fit,
    METH_VARARGS,
    PyDoc_STR("kmeans Algorithem")},
    {"wam_fit",
    (PyCFunction) wam_fit,
    METH_VARARGS,
    PyDoc_STR("Weighted Adjacency Matrix Computaion")},
    {"ddg_fit",
    (PyCFunction) ddg_fit,
    METH_VARARGS,
    PyDoc_STR("Diagonal Degree Matrix Computaion")},
    {"lnorm_fit",
    (PyCFunction) lnorm_fit,
    METH_VARARGS,
    PyDoc_STR("Normalized Graph Laplacian Matrix Computaion")},
    {"jacobi_fit",
    (PyCFunction) jacobi_fit,
    METH_VARARGS,
    PyDoc_STR("jacobi Alogorithem")},
    {NULL, NULL, 0, NULL}

};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}

