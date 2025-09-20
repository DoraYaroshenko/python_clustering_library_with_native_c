#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/*checked and workes -> printMatrix works, matrixMultiplication works*/
void testMatrixMultiplication()
{
    matrix A = createMatrix(3,4);
    matrix B = createMatrix(4,3);
    matrix res;
    int i;
    int j;
    printf("Testing matrix multiplication and printMatrix\n");
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
    printf("Testing distance and printVector\n");
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
    printf("Testing similarityMatrix\n");
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
    freeDataPoints(points);
    freeMatrix(similarityMat);
}

/*tested -> diagonalDegreeMatrix works*/
void testDiagonalDegreeMatrix()
{
    int i;
    int j;
    matrix diagonalDegreeMat;
    dataPoints points = createDataPoints(3,4);
    printf("Testing diagonalDegreeMatrix\n");
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
    matrix normalizedMat;
    dataPoints points = createDataPoints(3,4);
    printf("Testing normalizedSimilarityMatrix\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    normalizedMat = normalizedSimilarityMatrix(points);
    printMatrix(normalizedMat);
    freeDataPoints(points);
    freeMatrix(normalizedMat);
}
/*tested -> transpose works*/
void testTranspose()
{
    matrix A = createMatrix(3,4);
    matrix res;
    int i;
    int j;
    printf("Testing transpose\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            A.matrixEntries[i][j] = i * j + j;
        }
    }
    printf("Original matrix:\n");
    printMatrix(A);
    res = transpose(A);
    printf("Transposed matrix:\n");
    printMatrix(res);
    freeMatrix(res);
    freeMatrix(A);
}

void testTrace()
{
    matrix A = createMatrix(3,3);
    double t;
    int i;
    int j;
    printf("Testing trace\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            A.matrixEntries[i][j] = i * j + j;
        }
    }
    printf("Original matrix:\n");
    printMatrix(A);
    t = trace(A);
    printf("%f\n", t);
    freeMatrix(A);
}

/*tested - substractMatrices works*/
void testSubstractMatrices()
{
    matrix A = createMatrix(3,4);
    matrix B = createMatrix(3,4);
    matrix res;
    int i;
    int j;
    printf("Testing substractMatrices\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            A.matrixEntries[i][j] = i * j + j;
        }
    }
    printMatrix(A);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            B.matrixEntries[i][j] = 1;
        }
    }
    printMatrix(B);
    res = substractMatrices(A, B);
    printMatrix(res);
    freeMatrix(res);
    freeMatrix(A);
    freeMatrix(B);
}

/*tester -> updateH works, iterateAlgorithm has segmentation fault*/
void testUpdateH()
{
    int i;
    int j;
    matrix H;
    matrix updatedH;
    matrix W;
    matrix result;
    dataPoints points = createDataPoints(3,4);
    printf("Testing updateH and iterateAlgorithm\n");
    printf("Datapoints:\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    W = normalizedSimilarityMatrix(points);
    printf("W:\n");
    printMatrix(W);
    H = createMatrix(3,2);
    for(i=0;i<3;i++){
        for(j=0;j<2;j++){
            H.matrixEntries[i][j] = 0.7929;
        }
    }
    printf("H:\n");
    printMatrix(H);
    updatedH = updateH(H,W);
    printf("Updated H:\n");
    printMatrix(updatedH);
    result = iterateAlgorithm(H,W);
    freeMatrix(updatedH);
    printf("Final H:\n");
    printMatrix(result);
    freeMatrix(H);
    freeMatrix(W);
    freeDataPoints(points);
}