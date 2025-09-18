#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

all_vecs *parseDataPoints(PyObject *points)
{
    all_vecs *all_vectors = (all_vecs*)malloc(sizeof(all_vecs));
    int N = PyList_Size(points);
    if (N == 0)
        return NULL;
    PyObject *first_point = PyList_GetItem(points, 0);
    if (!PyList_Check(first_point))
        return NULL;
    int dim = PyList_Size(first_point);

    all_vectors->num_vectors = N;
    all_vectors->all_vectors = (vector *)malloc(sizeof(vector) * N);
    if (all_vectors->all_vectors == NULL)
        return NULL;
    for (int i = 0; i < N; i++)
    {
        PyObject *point = PyList_GetItem(points, i);
        if (!PyList_Check(point) || PyList_Size(point) != dim)
            return NULL;
        all_vectors->all_vectors[i].dimension = dim;
        all_vectors->all_vectors[i].coordinates = (double *)malloc(sizeof(double) * dim);
        if (all_vectors->all_vectors[i].coordinates == NULL)
            return NULL;
        for (int j = 0; j < dim; j++)
        {
            PyObject *coord = PyList_GetItem(point, j);
            if (!PyFloat_Check(coord))
                return NULL;
            all_vectors->all_vectors[i].coordinates[j] = PyFloat_AsDouble(coord);
        }
    }
    return all_vectors;
}

double **parseMatrix(PyObject *matrix)
{
    int n = PyList_Size(matrix);
    double **parsedMatrix = (double **)malloc(sizeof(double *) * n);
    if (parsedMatrix == NULL)
        return NULL;
    for (int i = 0; i < n; i++)
    {
        PyObject *row = PyList_GetItem(matrix, i);
        if (!PyList_Check(row))
            return NULL;
        int k = PyList_Size(row);
        parsedMatrix[i] = (double *)malloc(sizeof(double) * k);
        for (int j = 0; j < k; j++)
        {
            PyObject *entry = PyList_GetItem(row, j);
            if (!PyFloat_Check(entry))
                return NULL;
            parsedMatrix[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return parsedMatrix;
}

PyObject *parseResultMatrix(double **matrix, int n, int k)
{
    PyObject *resultMat = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *row = PyList_New(k);
        for (int j = 0; j < k; j++)
        {
            PyList_SetItem(row, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(resultMat, i, row);
    }
    return resultMat;
}

static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;
    if (!PyArg_ParseTuple(args, "OO", &H, &W))
        return NULL;
    if (!PyList_Check(H) || !PyList_Check(W))
        return NULL;
    double **parsedH = parseMatrix(H);
    double **parsedW = parseMatrix(W);
    int rowsH = PyList_Size(H);
    int colsH = PyList_Size(PyList_GetItem(H, 0));
    /*int rowsW = PyList_Size(W);
    int colsW = PyList_Size(PyList_GetItem(W, 0));
    */
    double** result = iterateAlgorithm(parsedH, parsedW, rowsH, colsH);
    return parseResultMatrix(result, rowsH, colsH);
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PyList_Check(points))
    //     return NULL;
    all_vecs *dataPoints = parseDataPoints(points);
    int n = dataPoints->num_vectors;
    double **similarityMat = similarityMatrix(*dataPoints);
    PyObject *result = parseResultMatrix(similarityMat, n, n);
    freeMatrix(similarityMat, n);
    freeMemory(dataPoints, n);
    return result;
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PyList_Check(points))
    //     return NULL;
    all_vecs *dataPoints = parseDataPoints(points);
    int n = dataPoints->num_vectors;
    double **diagonalDegreeMat = diagonalDegreeMatrix(*dataPoints);
    PyObject *result = parseResultMatrix(diagonalDegreeMat, n, n);
    freeMatrix(diagonalDegreeMat, n);
    freeMemory(dataPoints,n);
    return result;
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PyList_Check(points))
    //     return NULL;
    all_vecs *dataPoints = parseDataPoints(points);
    int n = dataPoints->num_vectors;
    double **normalizedMat = normalizedSimilarityMatrix(*dataPoints);
    PyObject *result = parseResultMatrix(normalizedMat, n, n);
    freeMatrix(normalizedMat, n);
    freeMemory(dataPoints, n);
    return result;
}

static PyMethodDef symnmfMethods[] = {
    {"symnmf",
     (PyCFunction)symnmf,
     METH_VARARGS,
     PyDoc_STR("returns the final H")},
    {"sym",
     (PyCFunction)sym,
     METH_VARARGS,
     PyDoc_STR("returns the similarity matrix")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR("returns the diagonal degree matrix")},
    {"norm",
     (PyCFunction)norm,
     METH_VARARGS,
     PyDoc_STR("returns the normalized similarity matrix")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",
    NULL,
    -1,
    symnmfMethods};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}