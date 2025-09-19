#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/*checked and workes -> printMatrix works, matrixMultiplication works*/
void testMatrixMultiplication()
{
    printf("Testing matrix multiplication and printMatrix");
    matrix A = createMatrix(3,4);
    matrix B = createMatrix(4,3);
    matrix res;
    int i;
    int j;
    for (i = 0; i < A.numOfRows; i++)
    {
        for (j = 0; j < A.numOfCols; j++)
        {
            A.matrixEntries[i][j] = i * j + j;
        }
    }
    printMatrix(A);
    for (i = 0; i < B.numOfRows; i++)
    {
        for (j = 0; j < B.numOfCols; j++)
        {
            B.matrixEntries[i][j] = 1;
        }
    }
    printMatrix(B);
    res = matrixMultiplication(A, B);
    printMatrix(res);
    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(res);
}

/*checked and works -> distance works and print vector works*/
void testDistance()
{
    vector a = createVector(3);
    vector b = createVector(3);
    int i;
    for (i = 0; i < a.dimension; i++)
    {
        a.coordinates[i] = i;
        b.coordinates[i] = 2 * i;
    }
    printVector(a);
    printVector(b);
    printf("%f", distance(a, b));
    freeVector(a);
    freeVector(b);
}

/*tested -> similarityMatrix works*/
void testSimilarityMatrix()
{
    int i;
    int j;
    matrix similarityMat;
    dataPoints points = createDataPoints(3,4);
    for (i = 0; i < points.num_vectors; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    similarityMat = similarityMatrix(points);
    printMatrix(similarityMat);
    freeMemory(points);
    freeMatrix(similarityMat);
}

/*tested -> diagonalDegreeMatrix works*/
void testDiagonalDegreeMatrix()
{
    int i;
    int j;
    matrix diagonalDegreeMat;
    dataPoints points = createDataPoints(3,4);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    diagonalDegreeMat = diagonalDegreeMatrix(points);
    printMatrix(diagonalDegreeMat);
    freeMatrix(diagonalDegreeMat);
}

void testNormalizedSimilarityMatrix()
{
    int i;
    int j;
    double **normalizedMat;
    all_vecs points;
    points.all_vectors = (vector *)malloc(sizeof(vector) * 3);
    points.num_vectors = 3;
    for (i = 0; i < 3; i++)
    {
        points.all_vectors[i].dimension = 4;
        points.all_vectors[i].coordinates = (double *)malloc(sizeof(double) * 4);
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(&points.all_vectors[i]);
    }
    normalizedMat = normalizedSimilarityMatrix(points);
    printMatrix(normalizedMat, 3, 3);
    freeMemory(&points, 3);
    freeMatrix(normalizedMat, 3);
}
/*tested -> transpose works*/
void testTranspose()
{
    double **A = (double **)malloc(sizeof(double *) * 3);
    double **res;
    int i;
    int j;
    for (i = 0; i < 3; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * 4);
        for (j = 0; j < 4; j++)
        {
            A[i][j] = i * j + j;
        }
    }
    printf("Original matrix:\n");
    printMatrix(A, 3, 4);
    res = transpose(A, 3, 4);
    printf("Transposed matrix:\n");
    printMatrix(res, 4, 3);
    freeMatrix(res, 4);
    freeMatrix(A, 3);
}

void testTrace()
{
    double **A = (double **)malloc(sizeof(double *) * 3);
    double t;
    int i;
    int j;
    for (i = 0; i < 3; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * 3);
        for (j = 0; j < 3; j++)
        {
            A[i][j] = i * j + j;
        }
    }
    printf("Original matrix:\n");
    printMatrix(A, 3, 3);
    t = trace(A, 3);
    printf("%f\n", t);
    freeMatrix(A, 3);
}

/*tested - substractMatrices works*/
void testSubstractMatrices()
{
    double **A = (double **)malloc(sizeof(double *) * 3);
    double **B = (double **)malloc(sizeof(double *) * 3);
    double **res;
    int i;
    int j;
    for (i = 0; i < 3; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * 4);
        for (j = 0; j < 4; j++)
        {
            A[i][j] = i * j + j;
        }
    }
    printMatrix(A, 3, 4);
    for (i = 0; i < 3; i++)
    {
        B[i] = (double *)malloc(sizeof(double) * 4);
        for (j = 0; j < 4; j++)
        {
            B[i][j] = 1;
        }
    }
    printMatrix(B, 3, 4);
    res = substractMatrices(A, B, 3, 4);
    printMatrix(res, 3, 4);
    freeMatrix(res, 3);
    freeMatrix(A, 3);
    freeMatrix(B, 3);
}

/*tester -> updateH works, iterateAlgorithm has segmentation fault*/
void testUpdateH()
{
    int i;
    int j;
    double **H;
    double **updatedH;
    double **W;
    double **result;
    all_vecs points;
    points.all_vectors = (vector *)malloc(sizeof(vector) * 3);
    points.num_vectors = 3;
    printf("Datapoints:\n");
    for (i = 0; i < 3; i++)
    {
        points.all_vectors[i].dimension = 4;
        points.all_vectors[i].coordinates = (double *)malloc(sizeof(double) * 4);
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(&points.all_vectors[i]);
    }
    W = normalizedSimilarityMatrix(points);
    printf("W:\n");
    printMatrix(W,3,3);
    H = (double**)malloc(sizeof(double*)*3);
    for(i=0;i<3;i++){
        H[i] = (double*)malloc(sizeof(double)*2);
        for(j=0;j<2;j++){
            H[i][j] = 0.7929;
        }
    }
    printf("H:\n");
    printMatrix(H,3,2);
    updatedH = updateH(H,W,3,2);
    printf("Updated H:\n");
    printMatrix(updatedH,3,2);
    result = iterateAlgorithm(H,W,3,2);
    freeMatrix(updatedH,3);
    printf("Final H:\n");
    printMatrix(result,3,2);
    freeMatrix(updatedH,3);
    freeMatrix(H,3);
    freeMatrix(W,3);
}