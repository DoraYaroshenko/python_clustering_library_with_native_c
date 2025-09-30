#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/*checked and workes -> printMatrix works, matrixMultiplication works*/
void testMatrixMultiplication()
{
    matrix A;
    matrix B;
    matrix res;
    int i;
    int j;
    if(!initializeMatrix(3,4, &A)) errorHandling();
    if(!initializeMatrix(4,3,&B)) errorHandling();
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
    if(!matrixMultiplication(A, B,&res)) errorHandling();
    printMatrix(res);
    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(res);
}

/*checked and works -> distance works and print vector works*/
void testDistance()
{
    vector a;
    vector b;
    int i;
    if(!initializeVector(3,&a)) errorHandling();
    if(!initializeVector(3,&b)) errorHandling();
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
    dataPoints points;
    if(!initializeDataPoints(3,4,&points)) errorHandling();
    printf("Testing similarityMatrix\n");
    for (i = 0; i < points.num_vectors; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    if(!similarityMatrix(points,&similarityMat)) errorHandling();
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
    dataPoints points;
    if(!initializeDataPoints(3,4,&points)) errorHandling();
    printf("Testing diagonalDegreeMatrix\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    if(!diagonalDegreeMatrix(points,&diagonalDegreeMat)) errorHandling();
    printMatrix(diagonalDegreeMat);
    freeMatrix(diagonalDegreeMat);
}

void testNormalizedSimilarityMatrix()
{
    int i;
    int j;
    matrix normalizedMat;
    dataPoints points;
    if(!initializeDataPoints(3,4,&points)) errorHandling();
    printf("Testing normalizedSimilarityMatrix\n");
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
        {
            points.all_vectors[i].coordinates[j] = i + j;
        }
        printVector(points.all_vectors[i]);
    }
    if(!normalizedSimilarityMatrix(points,&normalizedMat)) errorHandling();
    printMatrix(normalizedMat);
    freeDataPoints(points);
    freeMatrix(normalizedMat);
}
/*tested -> transpose works*/
void testTranspose()
{
    matrix A;
    matrix res;
    int i;
    int j;
    if(!initializeMatrix(3,4,&A)) errorHandling();
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
    if(!transpose(A,&res)) errorHandling();
    printf("Transposed matrix:\n");
    printMatrix(res);
    freeMatrix(res);
    freeMatrix(A);
}

void testTrace()
{
    matrix A;
    double t;
    int i;
    int j;
    if(!initializeMatrix(3,3,&A)) errorHandling();
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
    matrix A;
    matrix B;
    matrix res;
    int i;
    int j;
    if(!initializeMatrix(3,4,&A)) errorHandling();
    if(!initializeMatrix(3,4,&B)) errorHandling();
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
    if(!substractMatrices(A, B,&res)) errorHandling();
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
    dataPoints points;
    if(!initializeDataPoints(3,4,&points)) errorHandling();
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
    if(!normalizedSimilarityMatrix(points,&W)) errorHandling();
    printf("W:\n");
    printMatrix(W);
    if(!initializeMatrix(3,2,&H)) errorHandling();
    for(i=0;i<3;i++){
        for(j=0;j<2;j++){
            H.matrixEntries[i][j] = 0.7929;
        }
    }
    printf("H:\n");
    printMatrix(H);
    if(!updateH(H,W,&updatedH)) errorHandling();
    printf("Updated H:\n");
    printMatrix(updatedH);
    if(!iterateAlgorithm(&H,W)) errorHandling();
    freeMatrix(updatedH);
    printf("Final H:\n");
    printMatrix(H);
    freeMatrix(H);
    freeMatrix(W);
    freeDataPoints(points);
}