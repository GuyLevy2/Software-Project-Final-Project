#define PY_SSIZE_T_CLEAN
#include "spkmeans.c"
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

static PyObject* eigengapHeuristic_fit(PyObject *self, PyObject *args){
    int N, k, i;
    double *eigenValues = NULL;
    PyObject *eigenValues_obj;
    
    /* Get input */
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

    /* body */
    k = eigenGap(N, &eigenValues);

    /* Free eigenValues */
    free(eigenValues);

    /* output */
    if(k == -1){
        return Py_BuildValue("");
    }

    return Py_BuildValue("i", k);
}

static PyObject* kmeans_fit(PyObject *self, PyObject *args){
    int k, dimension, line_count, maxIter, kmeans_success; 
    int i,j;
    double EPSILON;
    double **vectorsList;
    double **centroids_list;
    PyObject *vec_list_obj;
    PyObject *centroids_list_obj; /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iiiidO", &k, &dimension, &line_count, &maxIter, &EPSILON, &vec_list_obj)) {
        return Py_BuildValue("");
    }

    vectorsList = *(Create_C_Mat_From_PyObj(N, dimension, vec_list_obj));
    if(vectorsList == NULL){
        return Py_BuildValue("");
    }
    
    /* body */
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

    /* output */
    centroids_list_obj = Create_PyObj_Mat_From_C(k, dimension, &centroids_list);

    /* free centroids_list */
    freeMat(k, &centroids_list);

    return Py_BuildValue("O", centroids_list_obj);
}

static PyObject* wam_fit(PyObject *self, PyObject *args){
    int N, dimension, success;
    double **vectorsList = NULL, **wamMat = NULL;
    PyObject *vec_list_obj;     /* input */
    PyObject *wamMatrix_obj;    /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iiO", &dimension, &N, &vec_list_obj)) {
        return Py_BuildValue("");
    }

    vectorsList = *(Create_C_Mat_From_PyObj(N, dimension, vec_list_obj));
    if (vectorsList == NULL){
        return Py_BuildValue("");
    }

    /* Body */
    wamMat = *(initMat(N));
    if (wamMat == NULL){
        return Py_BuildValue("");
    }

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

static PyObject* ddg_fit(PyObject *self, PyObject *args){
    int N, success;
    double **ddgMat = NULL, **wamMat = NULL;
    PyObject *wamMat_obj;     /* input */
    PyObject *ddgMat_obj;    /* output */

    /* Get input */
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

    /* Build output */
    ddgMat_obj = Create_PyObj_Mat_From_C(N, N, &ddgMat);
    
    /* free ddgMat */
    freeMat(N, &ddgMat);

    return Py_BuildValue("O", ddgMat_obj);
}

static PyObject* lnorm_fit(PyObject *self, PyObject *args){
    int N, success;
    double **ddgMat = NULL, **wamMat = NULL, **lnormMat = NULL;
    PyObject *wamMat_obj, *ddgMat_obj;     /* input */
    PyObject *lnormMat_obj;    /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iOO", &N, &wamMat_obj, &ddgMat_obj)) {
        return Py_BuildValue("");
    }

    wamMat = *(Create_C_Mat_From_PyObj(N, N, wamMat_obj));
    if (wamMat == NULL){
        return Py_BuildValue("");
    }

    ddgMat = Create_C_Mat_From_PyObj(N, N, ddgMat_obj)
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

    /* Build output */
    lnormMat_obj = Create_PyObj_Mat_From_C(N, N, &lnormMat);
    
    /* free lnormMat */
    freeMat(N, &lnormMat);

    return Py_BuildValue("O", lnormMat_obj);
}

static PyObject* jacobi_fit(PyObject *self, PyObject *args){
    int N, success;
    double **symMat = NULL;
    double **eigenVectors = NULL;
    double *eigenValues = NULL;
    PyObject *symMat_obj;     /* input */
    PyObject *eigenVecors_obj, *eigenValues_obj;    /* output */

    /* Get input */
    if(!PyArg_ParseTuple(args, "iO", &N, &symMat_obj)) {
        return Py_BuildValue("");
    }

    symMat = *(Create_C_Mat_From_PyObj(N, N, symMat_obj));
    if (symMat == NULL){
        return Py_BuildValue("");
    }

    /* Body */
    eigenValues = (double*)malloc(N * sizeof(double));
    if (eigenValues == NULL){
        return Py_BuildValue("");
    }
    eigenvectors = *(initMat(N));
    if (eigenvectors == NULL){
        return Py_BuildValue("");
    }

    success = jacobi_func(N, &symMat, &eigenVectors, &eigenValues);
    if(success == 1){
        return Py_BuildValue("");
    }
    
    /* Free symMat */
    freeMat(N, &symMat);

    /* Build output */
    eigenVecors_obj = Create_PyObj_Mat_From_C(N, N, &eigenVectors);
    eigenValues_obj = Create_PyObj_Arr_From_C(N, &eigenValues);
    /* free eigenVectors and eigenValues */
    freeMat(N, &eigenVectors);
    free(eigenValues);

    return Py_BuildValue("(OO)", eigenValues_obj, eigenVecors_obj); /* Returns Tuple object */
}

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

PyObject* Create_PyObj_Arr_From_C(int numOfElements, double** c_arr){
    PyObject *outArr_obj; 
    int i;

    outArr_obj = PyList_New(numOfElements);
    for (i = 0; i < numOfElements; i++) { 
        PyList_SetItem(outArr_obj, i, Py_BuildValue("d", (*c_arr)[i]));
    }

    return outArr_obj;
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

