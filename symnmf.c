#define EPS 0.0001
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/*
    TODO: in all function make a cleanup goto if you fail to malloc stuff
    make sure that your 'free<*>' checks that the pointers are not NULL and frees if they are not
*/

/*
    TODO: make 0 (or -1) and 1 globals for FAIL and SUCCESS and return those,
    then maybe change !function() to function() != SUCCESS
*/

int initializeVector(int dim, vector *vec)
{
    vec->dimension = dim;
    vec->coordinates = (double *)calloc(dim, sizeof(double));
    if (vec->coordinates == NULL)
        return 0;
    return 1;
}

void freeVector(vector vec)
{
    free(vec.coordinates);
}

int initializeMatrix(int rows, int cols, matrix *newMatrix)
{
    int i;
    newMatrix->numOfRows = rows;
    newMatrix->numOfCols = cols;
    newMatrix->matrixEntries = (double **)calloc(newMatrix->numOfRows, sizeof(double *));
    if (newMatrix->matrixEntries == NULL)
        return 0;
    for (i = 0; i < newMatrix->numOfRows; i++)
    {
        newMatrix->matrixEntries[i] = (double *)calloc(newMatrix->numOfCols, sizeof(double));
        if (newMatrix->matrixEntries[i] == NULL)
            return 0;
    }
    return 1;
}

void freeMatrix(matrix mat)
{
    int i;
    for (i = 0; i < mat.numOfRows; i++)
    {
        if (mat.matrixEntries[i] != NULL)
        {
            free(mat.matrixEntries[i]);
        }
    }
    free(mat.matrixEntries);
}

int initializeDataPoints(int numOfVectors, int dim, dataPoints *points)
{
    int i;
    points->num_vectors = numOfVectors;
    points->all_vectors = (vector *)malloc(sizeof(vector) * numOfVectors);
    if (points->all_vectors == NULL)
        return 0;
    for (i = 0; i < points->num_vectors; i++)
    {
        vector new_vec; /* TODO: declare this at the start of the function */
        if (!initializeVector(dim, &new_vec))
            return 0;
        points->all_vectors[i] = new_vec;
    }
    return 1;
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

int matrixMultiplication(matrix a, matrix b, matrix *resMatrix)
{
    int i; /* TODO: declare i,j,k together? */
    int j;
    double sum;
    int k;
    if (!initializeMatrix(a.numOfRows, b.numOfCols, resMatrix))
        return 0;
    if (a.numOfCols != b.numOfRows)
        return 0;
    for (i = 0; i < a.numOfRows; i++)
    {
        for (j = 0; j < b.numOfCols; j++)
        {
            sum = 0;
            for (k = 0; k < b.numOfRows; k++)
            {
                sum += a.matrixEntries[i][k] * b.matrixEntries[k][j];
            }
            resMatrix->matrixEntries[i][j] = sum;
        }
    }
    return 1;
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
            printf("%.4f\n", (vec.coordinates)[i]);
        }
        else
        {
            printf("%.4f,", (vec.coordinates)[i]);
        }
    }
}

void printDataPoints(dataPoints points)
{
    int i;
    for (i = 0; i < points.num_vectors; i++)
    {
        printVector(points.all_vectors[i]);
    }
}

int transpose(matrix mat, matrix *transposedMatrix)
{
    int i;
    int j;
    if (!initializeMatrix(mat.numOfCols, mat.numOfRows, transposedMatrix))
        return 0;
    for (i = 0; i < mat.numOfCols; i++)
    {
        for (j = 0; j < mat.numOfRows; j++)
        {
            transposedMatrix->matrixEntries[i][j] = mat.matrixEntries[j][i];
        }
    }
    return 1;
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

int substractMatrices(matrix A, matrix B, matrix *result)
{
    int i;
    int j;
    if (!initializeMatrix(A.numOfRows, A.numOfCols, result))
        return 0;
    if (A.numOfRows != B.numOfRows || A.numOfCols != B.numOfCols)
        return 0;
    for (i = 0; i < A.numOfRows; i++)
    {
        for (j = 0; j < A.numOfCols; j++)
        {
            result->matrixEntries[i][j] = A.matrixEntries[i][j] - B.matrixEntries[i][j];
        }
    }
    return 1;
}

int updateH(matrix H, matrix W, matrix *updatedH)
{
    double beta = 0.5; /* TODO: make global constant in header */
    matrix WH;         /* TODO: better names */
    matrix transposedH;
    matrix HMulTransposedH;
    matrix HMulTransposedHMulH;
    int i;
    int j;
    if (!initializeMatrix(H.numOfRows, H.numOfCols, updatedH))
        goto cleanup;
    if (!matrixMultiplication(W, H, &WH))
        return 0;
    if (!transpose(H, &transposedH))
        return 0;
    if (!matrixMultiplication(H, transposedH, &HMulTransposedH))
        return 0;
    if (!matrixMultiplication(HMulTransposedH, H, &HMulTransposedHMulH))
        return 0;
    for (i = 0; i < H.numOfRows; i++)
    {
        for (j = 0; j < H.numOfCols; j++)
        {
            double denom = HMulTransposedHMulH.matrixEntries[i][j];
            if (denom == 0)
                denom += +0.000001;
            updatedH->matrixEntries[i][j] = H.matrixEntries[i][j] * (1 - beta + beta * (WH.matrixEntries[i][j] / denom));
        }
    }
    return 1;

cleanup:
    freeMatrix(WH);
    freeMatrix(transposedH);
    freeMatrix(HMulTransposedH);
    freeMatrix(HMulTransposedHMulH);
    return 0;
}

int iterateAlgorithm(matrix *H, matrix W)
{
    int max_iter = 300; /* make global constant */
    int i, j, l;
    matrix updatedH, updatedHMinusH, prevH;
    double frobeniusNorm;
    for (i = 0; i < max_iter; i++)
    {
        if (!updateH(*H, W, &updatedH))
            return 0;
        if (!substractMatrices(updatedH, *H, &updatedHMinusH))
            return 0;
        frobeniusNorm = 0;
        for (j = 0; j < H->numOfRows; j++)
        {
            for (l = 0; l < H->numOfCols; l++)
            {
                frobeniusNorm += pow(fabs(updatedHMinusH.matrixEntries[j][l]), 2);
            }
        }
        freeMatrix(updatedHMinusH);
        prevH = *H;
        *H = updatedH;
        freeMatrix(prevH);
        if (frobeniusNorm < EPS)
            break;
    }
    return 1;
}

void errorHandling()
{
    printf("An Error Has Occurred\n");
    exit(1);
}

int calculateNumOfPoints(char *filename)
{
    int numOfPoints = 0;
    double n;
    char c;
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        return -1;
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
    int dimension = 0;
    FILE *fp = fopen(filename, "r");
    double n;
    char c;
    if (!fp)
    {
        return -1;
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

int getInput(char *filename, dataPoints *points)
{
    double coordinate;
    char c;
    int numOfPoints = calculateNumOfPoints(filename);
    int dimension = calculateDimension(filename);
    int i = 0, j = 0;
    FILE *fp = fopen(filename, "r");
    if (numOfPoints == -1)
        return 0;
    if (dimension == -1)
        return 0;
    if (!initializeDataPoints(numOfPoints, dimension, points))
        return 0;
    if (!fp)
    {
        return 0;
    }
    while (fscanf(fp, "%lf%c", &coordinate, &c) == 2)
    {
        points->all_vectors[i].coordinates[j] = coordinate;
        j++;
        if (c == '\n')
        {
            i++;
            j = 0;
        }
    }
    fclose(fp);
    return 1;
}

int similarityMatrix(dataPoints points, matrix *outputMatrix)
{
    int n = points.num_vectors;
    int i;
    int j;
    if (!initializeMatrix(n, n, outputMatrix))
        return 0;
    for (i = 0; i < outputMatrix->numOfRows; i++)
    {
        for (j = 0; j < outputMatrix->numOfCols; j++)
        {
            if (j != i)
            {
                double exponent = -pow(distance(points.all_vectors[i], points.all_vectors[j]), 2) / 2;
                outputMatrix->matrixEntries[i][j] = exp(exponent);
            }
            else
                outputMatrix->matrixEntries[i][j] = 0;
        }
    }
    return 1;
}

int diagonalDegreeMatrix(dataPoints points, matrix *outputMatrix)
{
    int n = points.num_vectors;
    matrix similarityMat;
    int i;
    int j;
    double rowSum;
    if (!initializeMatrix(n, n, outputMatrix))
        return 0;
    if (!similarityMatrix(points, &similarityMat))
        return 0;
    for (i = 0; i < similarityMat.numOfRows; i++)
    {
        rowSum = 0;
        for (j = 0; j < similarityMat.numOfCols; j++)
        {
            rowSum += similarityMat.matrixEntries[i][j];
        }
        outputMatrix->matrixEntries[i][i] = rowSum;
    }
    freeMatrix(similarityMat);
    return 1;
}

int normalizedSimilarityMatrix(dataPoints points, matrix *resMatrix)
{
    matrix similarityMat;
    matrix degreeMat;
    int i;
    matrix mulMat;
    if (!similarityMatrix(points, &similarityMat))
        return 0;
    if (!diagonalDegreeMatrix(points, &degreeMat))
        return 0;
    for (i = 0; i < degreeMat.numOfRows; i++)
    {
        degreeMat.matrixEntries[i][i] = 1 / (sqrt(degreeMat.matrixEntries[i][i]));
    }
    if (!matrixMultiplication(degreeMat, similarityMat, &mulMat))
        return 0;
    if (!matrixMultiplication(mulMat, degreeMat, resMatrix))
        return 0;
    freeMatrix(mulMat);
    freeMatrix(similarityMat);
    freeMatrix(degreeMat);
    return 1;
}

int main(int argc, char **argv)
{
    char *goal;
    char *path;
    dataPoints points;
    matrix outputMatrix;
    if (argc != 3)
        return (1);
    goal = argv[1];
    path = argv[2];
    /* TODO: rename to 'readPointsFromFile' since this is what it does */
    if (!getInput(path, &points))
        errorHandling();

    if (strcmp(goal, "sym") == 0)
    {
        if (!similarityMatrix(points, &outputMatrix))
            errorHandling();
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        if (!diagonalDegreeMatrix(points, &outputMatrix))
            errorHandling();
    }
    else if (strcmp(goal, "norm") == 0)
    {
        if (!normalizedSimilarityMatrix(points, &outputMatrix))
            errorHandling();
    }
    else
        errorHandling();
    printMatrix(outputMatrix);
    freeMatrix(outputMatrix);
    freeDataPoints(points);
    return (0);
}