#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

dataPoints parseDataPoints(PyObject *points)
{
    int N = PySequence_Size(points);
    PyObject *first_point = PySequence_GetItem(points, 0);
    // if (!PySequence_Check(first_point))
    //     return NULL
    int dim = PySequence_Size(first_point);
    dataPoints all_vectors = createDataPoints(N,dim);
    for (int i = 0; i < all_vectors.num_vectors; i++)
    {
        PyObject *point = PySequence_GetItem(points, i);
        // if (!PySequence_Check(point) || PySequence_Size(point) != dim)
        //     return NULL;
        for (int j = 0; j < dim; j++)
        {
            PyObject *coord = PySequence_GetItem(point, j);
            // if (!PyFloat_Check(coord))
            //     return NULL;
            all_vectors.all_vectors[i].coordinates[j] = PyFloat_AsDouble(coord);
        }
    }
    return all_vectors;
}

matrix parseMatrix(PyObject *mat)
{
    int n = PySequence_Size(mat);
    PyObject *row = PySequence_GetItem(mat, 0);
    // if (!PySequence_Check(row))
    //     return NULL;
    int dim = PySequence_Size(row);
    matrix parsedMatrix = createMatrix(n,dim);
    for (int i = 0; i < n; i++)
    {
        row = PySequence_GetItem(mat, i);
        for (int j = 0; j < dim; j++)
        {
            PyObject *entry = PySequence_GetItem(row, j);
            // if (!PyFloat_Check(entry))
            //     return NULL;
            parsedMatrix.matrixEntries[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return parsedMatrix;
}

PyObject *parseResultMatrix(matrix mat)
{
    PyObject *resultMat = PyList_New(mat.numOfRows);
    for (int i = 0; i < mat.numOfRows; i++)
    {
        PyObject *row = PyList_New(mat.numOfCols);
        for (int j = 0; j < mat.numOfCols; j++)
        {
            PyList_SetItem(row, j, PyFloat_FromDouble(mat.matrixEntries[i][j]));
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
    // if (!PySequence_Check(H) || !PySequence_Check(W))
    //     return NULL;
    matrix parsedH = parseMatrix(H);
    matrix parsedW = parseMatrix(W);
    matrix result = iterateAlgorithm(parsedH, parsedW);
    return parseResultMatrix(result);
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PySequence_Check(points))
    //     return NULL;
    dataPoints vectors = parseDataPoints(points);
    printDataPoints(vectors);
    matrix similarityMat = similarityMatrix(vectors);
    PyObject *result = parseResultMatrix(similarityMat);
    freeMatrix(similarityMat);
    freeDataPoints(vectors);
    return result;
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PySequence_Check(points))
    //     return NULL;
    dataPoints vectors = parseDataPoints(points);
    matrix diagonalDegreeMat = diagonalDegreeMatrix(vectors);
    PyObject *result = parseResultMatrix(diagonalDegreeMat);
    freeMatrix(diagonalDegreeMat);
    freeDataPoints(vectors);
    return result;
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    // if (!PySequence_Check(points))
    //     return NULL;
    dataPoints vectors = parseDataPoints(points);
    matrix normalizedMat = normalizedSimilarityMatrix(vectors);
    PyObject *result = parseResultMatrix(normalizedMat);
    freeMatrix(normalizedMat);
    freeDataPoints(vectors);
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
    "_symnmf",
    NULL,
    -1,
    symnmfMethods};

PyMODINIT_FUNC PyInit__symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}