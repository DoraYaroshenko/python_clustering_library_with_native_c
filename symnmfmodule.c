# define PY_SSIZE_T_CLEAN
# include <Python.h>

static PyMethodDef symnmfMethods[] = {
    {"fit",
     (PyCFunction) fit,
     METH_VARARGS,
     PyDoc_STR("run kmeans on centroids")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_mysymnmfpp(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}