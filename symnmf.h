typedef struct
{
    double *coordinates;
    int dimension;
} vector;

typedef struct
{
    vector *all_vectors;
    int num_vectors;
} dataPoints;

typedef struct
{
    double **matrixEntries;
    int numOfRows;
    int numOfCols;
} matrix;

vector createVector(int dim);
void freeVector(vector vec);
dataPoints createDataPoints(int n, int d);
void freeDataPoints(dataPoints all_vectors);
matrix createMatrix(int rows, int cols);
void freeMatrix(matrix mat);
matrix matrixMultiplication(matrix a, matrix b);
double distance(vector v1, vector v2);
void printMatrix(matrix mat);
void printVector(vector vec);
void errorHandling();
dataPoints getInput();
matrix similarityMatrix(dataPoints points);
matrix diagonalDegreeMatrix(dataPoints points);
matrix normalizedSimilarityMatrix(dataPoints points);
matrix transpose(matrix mat);
double trace(matrix mat);
matrix substractMatrices(matrix a, matrix b);
matrix updateH(matrix H, matrix W);
matrix iterateAlgorithm(matrix H, matrix W);
int calculateNumOfPoints(char *filename);
int calculateDimension(char *filename);

void testMatrixMultiplication();
void testDistance();
void testSimilarityMatrix();
void testDiagonalDegreeMatrix();
void testNormalizedSimilarityMatrix();
void testTranspose();
void testTrace();
void testSubstractMatrices();
void testUpdateH();