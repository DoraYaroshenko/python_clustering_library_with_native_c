#define PY_SSIZE_T_CLEAN
#define EPS 0.0001
#include <Python.h>
#include "symnmf.h"
#include "symnmf.c"

all_vecs parseDataPoints(PyObject *points)
{
    all_vecs all_vectors;
    int N = PyList_Size(points);
    if (N == 0)
        return NULL;

    PyObject *first_point = PyList_GetItem(points, 0);
    if (!PyList_Check(first_point))
        return NULL;
    int dim = PyList_Size(first_point);

    all_vectors.num_vectors = N;
    all_vectors.all_vectors = (vector *)malloc(sizeof(vector) * N);
    if (all_vectors.all_vectors == NULL)
        return NULL;
    for (int i = 0; i < N; i++)
    {
        PyObject *point = PyList_GetItem(points, i);
        if (!PyList_Check(point) || PyList_Size(point) != dim)
            return NULL;
        all_vectors.all_vectors[i].dimension = dim;
        all_vectors.all_vectors[i].coordinates = (double *)malloc(sizeof(double) * dim);
        if (all_vectors.all_vectors[i].coordinates == NULL)
            return NULL;
        for (int j = 0; j < dim; j++)
        {
            PyObject *coord = PyList_GetItem(point, j);
            if (!PyFloat_Check(coord))
                return NULL;
            all_vectors.all_vectors[i].coordinates[j] = PyFloat_AsDouble(coord);
        }
    }
    return all_vectors;
}

float **parseMatrix(PyObject *matrix)
{
    int n = PyList_Size(matrix);
    float **parsedMatrix = (float **)malloc(sizeof(float *) * n);
    if (parsedMatrix == NULL)
        return NULL;
    for (int i = 0; i < n; i++)
    {
        PyObject *row = PyList_GetItem(matrix, i);
        if (!PyList_Check(row) || PyList_Size(row) != n)
            return NULL;
        int k = PyList_Size(row);
        parsedMatrix[i] = (float *)malloc(sizeof(float) * k);
        for (int j = 0; j < k; j++)
        {
            float entry = PyList_GetItem(row, j);
            if (!PyFloat_Check(entry))
                return NULL;
            parsedMatrix[i][j] = PyFloat_AsDouble(entry);
        }
    }
    return parsedMatrix;
}

PyObject *parseResultMatrix(float **matrix, int n, int k)
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

void freeDataPoints(all_vecs *points)
{
    int n = points->num_vectors;
    for (int i = 0; i < n; i++)
    {
        free(points->all_vectors[i]);
    }
    free(points);
}

float **transpose(float **matrix, int rows, int cols)
{
    float **transposedMatrix = (float **)malloc(sizeof(float *) * cols);
    for (int i = 0; i < cols; i++)
    {
        float *row = (float *)malloc(sizeof(float) * rows);
        for (int j = 0; j < rows; j++)
        {
            transposedMatrix[i][j] = matrix[j][i];
        }
    }
    return transposedMatrix;
}

float trace(float **matrix, int n){
    float sum = 0;
    for(int i=0;i<n;i++){
        sum+=matrix[i][i];
    }
    return sum;
}

float **substractMatrices(float **A, float **B, int rows, int cols){
    float** result = (float**)malloc(sizeof(float*)*rows);
    for(int i=0;i<rows;i++){
        result[i] = (float*)malloc(sizeof(float)*cols);
        for(int j=0;j<cols;j++){
            result[i][j] = A[i][j]-B[i][j];
        }
    }
    return result;
}

float **updateH(float **H, float **W, int n, int k)
{
    float beta = 0.5;
    float **updatedH = (float **)malloc(sizeof(float *) * n);
    float **WH = matrixMultiplication(W, H, n, n, k);
    float **transposedH = transpose(H);
    float **HMulTransposedH = matrixMultiplication(H, transposedH, n, k, n);
    float **HMulTransposedHMulH = matrixMultiplication(HMulTransposedH, H, n, n, k);
    for (int i = 0; i < n; i++)
    {
        updatedH[i] = (float *)malloc(sizeof(float) * k);
        for (int j = 0; j < k; j++)
        {
            updatedH[i][j] = H[i][j]*(1-beta+beta*(WH[i][j]/HMulTransposedHMulH[i][j]));
        }
    }
    freeMatrix(WH);
    freeMatrix(transposedH);
    freeMatrix(HMulTransposedH);
    freeMatrix(HMulTransposedHMulH);
    return updatedH;
}

float **iterateAlgorithm(float** H, float** W, int n, int k){
    int max_iter = 300;
    for(int i=0;i<max_iter;i++){
        float** updatedH = updatedH(H,W,n,k);
        float** updatedHMinusH = substractMatrices(updatedH, H, n, k);
        float frobeniusNorm = 0;
        for(int j=0;j<n;j++){
            for(int l=0;l<k;l++){
                frobeniusNorm+=pow(abs(updatedHMinusH[j][l]),2);
            }
        }
        if(frobeniusNorm<EPS) break;
        freeMatrix(updatedHMinusH);
        freeMatrix(H);
        H = updatedH;
    }
    return H;
}

static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;
    if (!PyArg_ParseTuple(args, "OOid", &H, &W))
        return NULL;
    if (!PyList_Check(H) || !PyList_Check(W))
        return NULL;
    float **parsedH = parseMatrix(H);
    float **parsedW = parseMatrix(W);
    int rowsH = PyList_Size(H);
    int colsH = PyList_Size(PyList_GetItem(H,0));
    int rowsH = PyList_Size(W);
    int colsH = PyList_Size(PyList_GetItem(W,0));
    return iterateAlgorithm(parsedH,parsedW,rowsH,colsH);
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "OOid", &points))
        return NULL;
    if (!PyList_Check(points))
        return NULL;
    all_vecs dataPoints = parseDataPoints(points);
    int n = dataPoints.num_vectors;
    float **similarityMat = similarityMatrix(dataPoints);
    PyObject *result = parseResultMatrix(similarityMat, n, n);
    freeMatrix(similarityMat);
    freeDataPoints(dataPoints);
    return result;
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "OOid", &points))
        return NULL;
    if (!PyList_Check(points))
        return NULL;
    all_vecs dataPoints = parseDataPoints(points);
    int n = dataPoints.num_vectors;
    float **diagonalDegreeMat = diagonalDegreeMatrix(dataPoints);
    PyObject *result = parseResultMatrix(diagonalDegreeMat, n, n);
    freeMatrix(diagonalDegreeMat);
    freeDataPoints(dataPoints);
    return result;
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *points;
    if (!PyArg_ParseTuple(args, "OOid", &points))
        return NULL;
    if (!PyList_Check(points))
        return NULL;
    all_vecs dataPoints = parseDataPoints(points);
    int n = dataPoints.num_vectors;
    float **normalizedMat = normalizedSimilarityMatrix(dataPoints);
    PyObject *result = parseResultMatrix(normalizedMat, n, n);
    freeMatrix(normalizedMat);
    freeDataPoints(dataPoints);
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

PyMODINIT_FUNC PyInit_mysymnmfpp(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}