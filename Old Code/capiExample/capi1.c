#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>


static double geo_c(double z, int n)
{
    double geo_sum = 0;
    geo_sum = n*z;
    return geo_sum;
}

static PyObject* geo_capi(PyObject *self, PyObject *args)
{
    double z;
    int n;

    if(!PyArg_ParseTuple(args, "di", &z, &n)) {
        return NULL;
    }

    return Py_BuildValue("d", geo_c(z, n));
}

static PyMethodDef capiMethods[] = {
    {"geo",
    (PyCFunction) geo_capi,
    METH_VARARGS,
    PyDoc_STR("A geometric ....")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_demo1",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_capi_demo1(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}

