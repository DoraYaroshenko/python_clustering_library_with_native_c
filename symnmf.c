#define EPS 0.0001
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

vector createVector(int dim)
{
    vector newVec;
    newVec.dimension = dim;
    newVec.coordinates = (double *)calloc(newVec.dimension, sizeof(double));
    return newVec;
}

void freeVector(vector vec)
{
    free(vec.coordinates);
}

matrix createMatrix(int rows, int cols)
{
    int i;
    matrix newMatrix;
    newMatrix.numOfRows = rows;
    newMatrix.numOfCols = cols;
    newMatrix.matrixEntries = (double **)calloc(newMatrix.numOfRows, sizeof(double *));
    for (i = 0; i < newMatrix.numOfRows; i++)
    {
        newMatrix.matrixEntries[i] = (double *)calloc(newMatrix.numOfCols, sizeof(double));
    }
    return newMatrix;
}

void freeMatrix(matrix mat)
{
    int i;
    for (i = 0; i < mat.numOfRows; i++)
    {
        free(mat.matrixEntries[i]);
    }
    free(mat.matrixEntries);
}

dataPoints createDataPoints(int numOfVectors, int dim)
{
    int i;
    dataPoints allVectors;
    allVectors.num_vectors = numOfVectors;
    allVectors.all_vectors = (vector *)malloc(sizeof(vector) * numOfVectors);
    for (i = 0; i < allVectors.num_vectors; i++)
    {
        allVectors.all_vectors[i] = createVector(dim);
    }
    return allVectors;
}

void freeDataPoints(dataPoints all_vectors)
{
    int i;
    for (i = 0; i < all_vectors.num_vectors; i++)
    {
        freeVector(all_vectors.all_vectors[i]);
    }
    free(all_vectors.all_vectors);
}

matrix matrixMultiplication(matrix a, matrix b)
{
    matrix resMatrix = createMatrix(a.numOfRows, b.numOfCols);
    int i;
    int j;
    double sum;
    int k;
    if (a.numOfCols != b.numOfRows)
        errorHandling();
    for (i = 0; i < a.numOfRows; i++)
    {
        for (j = 0; j < b.numOfCols; j++)
        {
            sum = 0;
            for (k = 0; k < b.numOfRows; k++)
            {
                sum += a.matrixEntries[i][k] * b.matrixEntries[k][j];
            }
            resMatrix.matrixEntries[i][j] = sum;
        }
    }
    return resMatrix;
}

double distance(vector v1, vector v2)
{
    double sum = 0;
    int i;
    for (i = 0; i < v1.dimension; i++)
    {
        sum += pow(v1.coordinates[i] - v2.coordinates[i], 2);
    }
    return sqrt(sum);
}

void printMatrix(matrix mat)
{
    int i;
    int j;
    for (i = 0; i < mat.numOfRows; i++)
    {
        for (j = 0; j < mat.numOfCols; j++)
        {
            printf("%.4f", mat.matrixEntries[i][j]);
            if (j != mat.numOfCols - 1)
                printf(",");
            else
                printf("\n");
        }
    }
}

void printVector(vector vec)
{
    int i;
    for (i = 0; i < vec.dimension; i++)
    {
        if (i == vec.dimension - 1)
        {
            printf("%.4f", (vec.coordinates)[i]);
        }
        else
            printf("%.4f,", (vec.coordinates)[i]);
    }
    printf("\n");
}

void printDataPoints(dataPoints points){
    int i;
    for(i=0;i<points.num_vectors;i++){
        printVector(points.all_vectors[i]);
    }
}

matrix transpose(matrix mat)
{
    matrix transposedMatrix = createMatrix(mat.numOfCols, mat.numOfRows);
    int i;
    int j;
    for (i = 0; i < mat.numOfCols; i++)
    {
        for (j = 0; j < mat.numOfRows; j++)
        {
            transposedMatrix.matrixEntries[i][j] = mat.matrixEntries[j][i];
        }
    }
    return transposedMatrix;
}

double trace(matrix mat)
{
    double sum = 0;
    int i;
    for (i = 0; i < mat.numOfRows; i++)
    {
        sum += mat.matrixEntries[i][i];
    }
    return sum;
}

matrix substractMatrices(matrix A, matrix B)
{
    matrix result = createMatrix(A.numOfRows, A.numOfCols);
    int i;
    int j;
    if (A.numOfRows != B.numOfRows || A.numOfCols != B.numOfCols)
        errorHandling();
    for (i = 0; i < A.numOfRows; i++)
    {
        for (j = 0; j < A.numOfCols; j++)
        {
            result.matrixEntries[i][j] = A.matrixEntries[i][j] - B.matrixEntries[i][j];
        }
    }
    return result;
}

matrix updateH(matrix H, matrix W)
{
    double beta = 0.5;
    matrix updatedH = createMatrix(H.numOfRows, H.numOfCols);
    matrix WH = matrixMultiplication(W, H);
    matrix transposedH = transpose(H);
    matrix HMulTransposedH = matrixMultiplication(H, transposedH);
    matrix HMulTransposedHMulH = matrixMultiplication(HMulTransposedH, H);
    int i;
    int j;
    if (H.numOfRows != W.numOfRows)
        errorHandling();
    for (i = 0; i < H.numOfRows; i++)
    {
        for (j = 0; j < H.numOfCols; j++)
        {
            updatedH.matrixEntries[i][j] = H.matrixEntries[i][j] * (1 - beta + beta * (WH.matrixEntries[i][j] / HMulTransposedHMulH.matrixEntries[i][j]));
        }
    }
    freeMatrix(WH);
    freeMatrix(transposedH);
    freeMatrix(HMulTransposedH);
    freeMatrix(HMulTransposedHMulH);
    return updatedH;
}

matrix iterateAlgorithm(matrix H, matrix W)
{
    int max_iter = 300;
    int i, j, l;
    matrix updatedH, updatedHMinusH, prevH;
    double frobeniusNorm;
    for (i = 0; i < max_iter; i++)
    {
        updatedH = updateH(H, W);
        updatedHMinusH = substractMatrices(updatedH, H);
        frobeniusNorm = 0;
        for (j = 0; j < H.numOfRows; j++)
        {
            for (l = 0; l < H.numOfCols; l++)
            {
                frobeniusNorm += pow(fabs(updatedHMinusH.matrixEntries[j][l]), 2);
            }
        }
        if (frobeniusNorm < EPS)
            break;
        freeMatrix(updatedHMinusH);
        prevH = H;
        H = updatedH;
        freeMatrix(prevH);
    }
    return H;
}

void errorHandling()
{
    printf("An Error Has Occurred\n");
}

int calculateNumOfPoints(char *filename)
{
    int numOfPoints = 0;
    double n;
    char c;
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        errorHandling();
        exit(1);
    }
    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
            numOfPoints++;
    }
    fclose(fp);
    return numOfPoints;
}

int calculateDimension(char *filename)
{
    int dimension=0;
    FILE *fp = fopen(filename, "r");
    double n;
    char c;
    if (!fp)
    {
        errorHandling();
        exit(1);
    }
    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        dimension++;
        if (c == '\n')
            break;
    }
    fclose(fp);
    return dimension;
}

/*tested -> getInput works*/
dataPoints getInput(char *filename)
{
    double coordinate;
    char c;
    int numOfPoints = calculateNumOfPoints(filename);
    int dimension = calculateDimension(filename);
    int i = 0, j = 0;
    dataPoints points = createDataPoints(numOfPoints, dimension);
    FILE *fp = fopen(filename, "r");
    /*printf("num of points: %d\n",numOfPoints);
    printf("dimension: %d\n",dimension);*/
    if (!fp)
    {
        errorHandling();
        exit(1);
    }
    while (fscanf(fp, "%lf%c", &coordinate, &c) == 2)
    {
        points.all_vectors[i].coordinates[j] = coordinate;
        j++;
        if (c == '\n')
        {
            i++;
            j = 0;
        }
    }
    fclose(fp);
    return points;
}

matrix similarityMatrix(dataPoints points)
{
    int n = points.num_vectors;
    matrix outputMatrix = createMatrix(n, n);
    int i;
    int j;
    for (i = 0; i < outputMatrix.numOfRows; i++)
    {
        for (j = 0; j < outputMatrix.numOfCols; j++)
        {
            if (j != i)
            {
                double exponent = -pow(distance(points.all_vectors[i], points.all_vectors[j]), 2) / 2;
                /*printf("%.4f\n",exponent);*/ 
                outputMatrix.matrixEntries[i][j] = exp(exponent);
            }
            else
                outputMatrix.matrixEntries[i][j] = 0;
        }
    }
    return outputMatrix;
}

matrix diagonalDegreeMatrix(dataPoints points)
{
    int n = points.num_vectors;
    matrix outputMatrix = createMatrix(n, n);
    matrix similarityMat = similarityMatrix(points);
    int i;
    int j;
    double rowSum;
    for (i = 0; i < similarityMat.numOfRows; i++)
    {
        rowSum = 0;
        for (j = 0; j < similarityMat.numOfCols; j++)
        {
            rowSum += similarityMat.matrixEntries[i][j];
        }
        outputMatrix.matrixEntries[i][i] = rowSum;
    }
    freeMatrix(similarityMat);
    return outputMatrix;
}

matrix normalizedSimilarityMatrix(dataPoints points)
{
    matrix similarityMat = similarityMatrix(points);
    matrix degreeMat = diagonalDegreeMatrix(points);
    int i;
    matrix mulMat;
    matrix resMatrix;
    for (i = 0; i < degreeMat.numOfRows; i++)
    {
        degreeMat.matrixEntries[i][i] = 1 / sqrt(degreeMat.matrixEntries[i][i]);
    }
    mulMat = matrixMultiplication(degreeMat, similarityMat);
    resMatrix = matrixMultiplication(mulMat, degreeMat);
    freeMatrix(mulMat);
    freeMatrix(similarityMat);
    freeMatrix(degreeMat);
    return resMatrix;
}

int main(int argc, char **argv)
{
    char *goal;
    char *path;
    dataPoints points;
    /*int i;
    int n;*/
    matrix outputMatrix;
    if (argc != 3)
        return (1);
    goal = argv[1];
    /* printf("Goal: %s\n",goal);*/
    path = argv[2];
    /* printf("Path: %s\n",path);*/
    points = getInput(path);
    /*n = points.num_vectors;
    for(i=0;i<n;i++){
        printVector(points.all_vectors[i]);
        printf("\n");
    }*/
    if (strcmp(goal, "sym")==0)
    {
        /* printf("The goal is sym\n");*/
        outputMatrix = similarityMatrix(points);
    }
    if (strcmp(goal, "ddg")==0)
    {
        /* printf("The goal is ddg\n");*/
        outputMatrix = diagonalDegreeMatrix(points);
    }
    if (strcmp(goal, "norm")==0)
    {
        /* printf("The goal is norm\n");*/
        outputMatrix = normalizedSimilarityMatrix(points);
    }
    printMatrix(outputMatrix);
    freeMatrix(outputMatrix);
    freeDataPoints(points);
    return (0);
}

/*int main()
{
    testMatrixMultiplication();
    testDistance();
    testSimilarityMatrix();
    testDiagonalDegreeMatrix();
    testNormalizedSimilarityMatrix();
    testTranspose();
    testTrace();
    testSubstractMatrices();
    testUpdateH();
    return (0);
}*/