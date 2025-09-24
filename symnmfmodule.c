#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

/* TODO: keep conventions!!! (blank lines are not a bad idea also)*/

int parseDataPoints(PyObject *points, dataPoints *allVectors)
{
    int N = PySequence_Size(points);
    PyObject *first_point = PySequence_GetItem(points, 0);
    if (!PySequence_Check(first_point))
        return 0;
    int dim = PySequence_Size(first_point);
    if(!initializeDataPoints(N,dim, allVectors)) return 0;
    for (int i = 0; i < allVectors->num_vectors; i++)
    {
        PyObject *point = PySequence_GetItem(points, i);
        if (!PySequence_Check(point) || PySequence_Size(point) != dim)
            return 0;
        for (int j = 0; j < dim; j++)
        {
            PyObject *coord = PySequence_GetItem(point, j);
            allVectors->all_vectors[i].coordinates[j] = PyFloat_AsDouble(coord);
        }
    }
    return 1;
}

int parseMatrix(PyObject *mat, matrix *parsedMatrix)
{
    int n = PySequence_Size(mat);
    PyObject *row = PySequence_GetItem(mat, 0);
    if (!PySequence_Check(row))
        return 0;
    int dim = PySequence_Size(row);
    if(!initializeMatrix(n,dim,parsedMatrix)) return 0;;
    for (int i = 0; i < n; i++)
    {
        row = PySequence_GetItem(mat, i);
        for (int j = 0; j < dim; j++)
        {
            PyObject *entry = PySequence_GetItem(row, j);
            if (!PyFloat_Check(entry))
                return 0;
            parsedMatrix->matrixEntries[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return 1;
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
    if (!PySequence_Check(H) || !PySequence_Check(W))
        return NULL;
    matrix parsedH; 
    if(!parseMatrix(H,&parsedH)) return NULL;
    matrix parsedW;
    if(!parseMatrix(W, &parsedW)) return NULL;
    if(!iterateAlgorithm(&parsedH, parsedW)) return NULL;
    return parseResultMatrix(parsedH);
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;
    if (!PySequence_Check(points))
        return NULL;
    dataPoints vectors;
    if(!parseDataPoints(points,&vectors)) return NULL;
    matrix similarityMat;
    if(!similarityMatrix(vectors,&similarityMat)) return NULL;
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
    if (!PySequence_Check(points))
        return NULL;
    dataPoints vectors;
    if(!parseDataPoints(points,&vectors)) return NULL;
    matrix diagonalDegreeMat;
    if(!diagonalDegreeMatrix(vectors,&diagonalDegreeMat)) return NULL;
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
    if (!PySequence_Check(points))
        return NULL;
    dataPoints vectors;
    if(!parseDataPoints(points,&vectors)) return NULL;
    matrix normalizedMat;
    if(!normalizedSimilarityMatrix(vectors,&normalizedMat)) return NULL;
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