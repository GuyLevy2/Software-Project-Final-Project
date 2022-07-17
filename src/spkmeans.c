#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>


double calcDistSqrt(double* v1, double* v2, int d){
    double dist = 0;
    int i = 0;
    
    for(i = 0; i < d; ++i){
        dist += pow(v1[i]-v2[i], 2.0);
    }
    return dist;
}

int assignVectorToClosestCluster(double* vector, double**** clusterList, double** centroidsList, int** index_to_insert, int k, int d){
    double minDist = -1.0;
    int minDistIndex = -1;
    int i = 0;
    double dist = 0.0;
    double* cent_vec;
    int last_index;
    
    for (i = 0; i < k; ++i){
        cent_vec = centroidsList[i];
        dist = calcDistSqrt(vector, cent_vec, d);
        
        if ((dist < minDist) || (minDist == -1.0)){
            minDist = dist;
            minDistIndex = i;
        }
    }
    
    /* insert vector to the closest cluster list at the last index to insert of the list and update last index */
    last_index = (*index_to_insert)[minDistIndex];
    (*clusterList)[minDistIndex][last_index] = vector;
    (*index_to_insert)[minDistIndex] = last_index + 1;

    return 1;   
}

void resetClusters(int** index_to_insert, int k){
   int i;
   
   for(i = 0; i < k; i++){
       (*index_to_insert)[i] = 0;
   }
}

double updateCentroids(int dimension, int k, double*** centroidsList, double*** clusterList, int* index_to_insert){
    double maxDelta = -1.0;
    double centroidDelta = 0.0;
    double indexDelta = 0.0;
    int i;
    int j;
    int m;
    double indexAverage = 0.0;
    double indexSum = 0.0;

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
        
        centroidDelta = pow(centroidDelta, 0.5);

        if (centroidDelta > maxDelta || maxDelta == -1.0){
            maxDelta = centroidDelta;
        }
    }

    return maxDelta;
}

int freeMemory(double*** vectorsList, double**** clusterList, int** index_to_insert, int line_count, int k){
    int i;
    /* free vectors */
    for (i = 0; i < line_count; i++){
        free((*vectorsList)[i]);
    }
    
    /* free vectors list */
    free(*vectorsList);

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

/* CAPI */

static int kmeans_c(int k, int dimension, int line_count, int maxIter, double EPSILON, double** vectorsList, double*** centroidsList)
{
    double*** clusterList;
    int* index_to_insert;
    double** cluster;
    int i;
    double* vec;
    double* centroid_i;
    int j;
    double deltaMiu;
    int iteration_number;
    int freeSuccess;
    
    index_to_insert = malloc(k * sizeof(int));
    if (index_to_insert == NULL){
        return 0;
    }

    clusterList = malloc(k * sizeof(double**));
    if (clusterList == NULL){
        return 0;
    }
    

    for (i = 0; i < k; i++){ /* initializing cluster list, each cluster, and centroids list with first k vectors - the vectors of kmeans++ algorithem. */
        cluster = malloc(line_count * sizeof(double*));
        if (cluster == NULL){
            return 0;
        }

        clusterList[i] = cluster;
        
        vec = vectorsList[i];
        
        centroid_i = malloc(dimension * sizeof(double));
        if (centroid_i == NULL){
            return 0;
        }

        for (j = 0; j < dimension; j++){
            centroid_i[j] = vec[j];
        }
        
        (*centroidsList)[i] = centroid_i;
        index_to_insert[i] = 0;
    }
    
    deltaMiu = EPSILON*2;
    iteration_number = 0;

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
    freeSuccess = freeMemory(&vectorsList, &clusterList, &index_to_insert, line_count, k);
    if (freeSuccess == 0){
        return 0;
    }

    return 1;
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
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("kmeans algo")},
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


