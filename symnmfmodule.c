#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

/*parseDataPoints saves the given datapoints in a strucr dataPoints, so it can be used by other functions*/
int parseDataPoints(PyObject *points, dataPoints *allVectors)
{
    int i, j, dim, N;
    PyObject *firstPoint = PySequence_GetItem(points, 0);
    PyObject *point;
    PyObject *coordinate;
    if (!PySequence_Check(firstPoint))
        return FAIL;
    dim = PySequence_Size(firstPoint);
    N = PySequence_Size(points);
    if (initializeDataPoints(N, dim, allVectors) == FAIL)
        return FAIL;
    for (i = 0; i < allVectors->numOfPoints; i++)
    {
        point = PySequence_GetItem(points, i);
        if (!PySequence_Check(point) || PySequence_Size(point) != dim)
        {
            goto cleanup;
        }
        for (j = 0; j < dim; j++)
        {
            coordinate = PySequence_GetItem(point, j);
            if (!PyFloat_Check(coordinate))
                goto cleanup;
            allVectors->points[i].coordinates[j] = PyFloat_AsDouble(coordinate);
        }
    }
    return SUCCESS;
cleanup:
    freeDataPoints(*allVectors);
    return FAIL;
}

/*parseMatrix saves given matrix in a struct matrix*/
int parseMatrix(PyObject *mat, matrix *parsedMatrix)
{
    int n, dim, i, j;
    PyObject *row = PySequence_GetItem(mat, 0);
    PyObject *entry;
    if (!PySequence_Check(row))
        return FAIL;
    dim = PySequence_Size(row);
    n = PySequence_Size(mat);
    if (initializeMatrix(n, dim, parsedMatrix) == FAIL)
        return FAIL;
    for (i = 0; i < n; i++)
    {
        row = PySequence_GetItem(mat, i);
        if (!PySequence_Check(row))
        {
            goto cleanup;
        }
        for (j = 0; j < dim; j++)
        {
            entry = PySequence_GetItem(row, j);
            if (!PyFloat_Check(entry))
                goto cleanup;
            parsedMatrix->matrixEntries[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return SUCCESS;
cleanup:
    freeMatrix(*parsedMatrix);
    return FAIL;
}

/*parseResultMatrix converts the data from the struct matrix to a Python list of lists*/
PyObject *parseResultMatrix(matrix mat)
{
    int i, j;
    PyObject *row;
    PyObject *resultMat = PyList_New(mat.numOfRows);
    for (i = 0; i < mat.numOfRows; i++)
    {
        row = PyList_New(mat.numOfCols);
        for (j = 0; j < mat.numOfCols; j++)
        {
            PyList_SetItem(row, j, PyFloat_FromDouble(mat.matrixEntries[i][j]));
        }
        PyList_SetItem(resultMat, i, row);
    }
    return resultMat;
}

/*symnmf performs symNMF using iterateAlgorithm function from symnmf.c*/
static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject *H, *W, *result = NULL;
    matrix parsedH, parsedW;
    ;
    if (!PyArg_ParseTuple(args, "OO", &H, &W))
        return NULL;
    if (!PySequence_Check(H) || !PySequence_Check(W))
        return NULL;
    if (parseMatrix(H, &parsedH) == SUCCESS && parseMatrix(W, &parsedW) == SUCCESS && iterateAlgorithm(&parsedH, parsedW) == SUCCESS)
        result = parseResultMatrix(parsedH);
    freeMatrix(parsedH);
    freeMatrix(parsedW);
    return result;
}

/*sym calculates the similarity matrix*/
static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *points, *result = NULL;
    dataPoints vectors;
    matrix similarityMat;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    if (!PySequence_Check(points))
        return NULL;
    if (parseDataPoints(points, &vectors) == SUCCESS && similarityMatrix(vectors, &similarityMat) == SUCCESS)
        result = parseResultMatrix(similarityMat);
    freeMatrix(similarityMat);
    freeDataPoints(vectors);
    return result;
}


/*ddg calculates diagonal degree matrix*/
static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *points, *result = NULL;
    dataPoints vectors;
    matrix diagonalDegreeMat;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    if (!PySequence_Check(points))
        return NULL;
    if (parseDataPoints(points, &vectors) == SUCCESS && diagonalDegreeMatrix(vectors, &diagonalDegreeMat) == SUCCESS)
        result = parseResultMatrix(diagonalDegreeMat);
    freeMatrix(diagonalDegreeMat);
    freeDataPoints(vectors);
    return result;
}

/*calculates normalized similarity matrix*/
static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *points, *result = NULL;
    dataPoints vectors;
    matrix normalizedMat;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    if (!PySequence_Check(points))
        return NULL;
    if (parseDataPoints(points, &vectors) == SUCCESS && normalizedSimilarityMatrix(vectors, &normalizedMat) == SUCCESS)
        result = parseResultMatrix(normalizedMat);
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
    "_symnmf", /* TODO: maybe change python module name to symnmf_c */
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