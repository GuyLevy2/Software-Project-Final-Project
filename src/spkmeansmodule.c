#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

static PyObject* eigengapHeuristic_fit(PyObject *self, PyObject *args){

}

static PyObject* kmeans_fit(PyObject *self, PyObject *args){
    
}

static PyObject* wam_fit(PyObject *self, PyObject *args){
    
}

static PyObject* ddg_fit(PyObject *self, PyObject *args){
    
}

static PyObject* lnorm_fit(PyObject *self, PyObject *args){
    
}

static PyObject* jacobi_fit(PyObject *self, PyObject *args){
    
}

static PyObject* fit(PyObject *self, PyObject *args)
{
    int k;
    int dimension;
    int line_count;
    int maxIter;
    double EPSILON;
    int kmeans_success;
    
    double** vectorsList;
    double* vector;
    double** centroids_list; /* output */
    
    PyObject* vec_list_obj;
    PyObject* vector_obj;
    PyObject* coord_obj;
    
    PyObject* centroids_list_obj; /* output */
    PyObject* centroid_obj; /* output */

    int i,j;

    if(!PyArg_ParseTuple(args, "iiiidO", &k, &dimension, &line_count, &maxIter, &EPSILON, &vec_list_obj)) {
        return Py_BuildValue("");
    }

    vectorsList = (double**)malloc(sizeof(double**)*line_count);
    
    if(vectorsList == NULL)
        return Py_BuildValue("");
    
    for(i = 0; i < line_count; i++){
        vector = (double*)malloc(sizeof(double*)*dimension);
        
        if (vector == NULL)
            return Py_BuildValue("");
        
        vectorsList[i] = vector;
    }
    
    for (i = 0; i < line_count; i++) {
        vector_obj = PyList_GetItem(vec_list_obj, i);
        for (j = 0; j < dimension; j++){
            coord_obj = PyList_GetItem(vector_obj, j);
            vectorsList[i][j] = PyFloat_AsDouble(coord_obj);
        }
    }
    
    centroids_list = malloc(k * sizeof(double*));
    if (centroids_list == NULL){
        return Py_BuildValue("");
    }


    kmeans_success = kmeans_c(k, dimension, line_count, maxIter, EPSILON, vectorsList, &centroids_list);
    if (kmeans_success == 0){
        return Py_BuildValue("");
    }
    
    centroids_list_obj = PyList_New(k);

    for (i = 0; i < k; i++) { 
        centroid_obj = PyList_New(dimension);
        for (j = 0; j < dimension; j++){
            PyList_SetItem(centroid_obj, j, Py_BuildValue("d", centroids_list[i][j]));
        }
        PyList_SetItem(centroids_list_obj, i, centroid_obj);
    }
    
    /* free centroids_list */
    for(i = 0; i < k; i++){
        free(centroids_list[i]);
    }
    free(centroids_list);    

    return Py_BuildValue("O", centroids_list_obj);
}

static PyMethodDef capiMethods[] = {
    {"eigengapHeuristic_fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("Eigengap Heuristic Algorithem")},
    {"kmeans_fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("kmeans Algorithem")},
    {"wam_fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("Weighted Adjacency Matrix Computaion")},
    {"ddg_fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("Diagonal Degree Matrix Computaion")},
    {"lnorm_fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("Normalized Graph Laplacian Matrix Computaion")},
    {"jacobi_fit",
    (PyCFunction) fit,
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

