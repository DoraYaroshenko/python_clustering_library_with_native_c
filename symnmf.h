#define FAIL -1
#define SUCCESS 1
#define EPS 0.0001
#define BETA 0.5
#define MAXITER 300

typedef struct
{
    double *coordinates;
    int dimension;
} vector;

typedef struct
{
    vector *points;
    int numOfPoints;
} dataPoints;

typedef struct
{
    double **matrixEntries;
    int numOfRows;
    int numOfCols;
} matrix;

int initializeVector(int dim, vector * vec);
void freeVector(vector vec);
int initializeDataPoints(int numOfVectors, int dim, dataPoints* points);
void freeDataPoints(dataPoints allVectors);
int initializeMatrix(int rows, int cols, matrix *newMatrix);
void freeMatrix(matrix mat);
int matrixMultiplication(matrix a, matrix b, matrix *resMatrix);
double distance(vector v1, vector v2);
void printMatrix(matrix mat);
void printVector(vector vec);
void errorHandling();
int readPointsFromFile(char *filename, dataPoints *points);
int similarityMatrix(dataPoints points, matrix *outputMatrix);
int diagonalDegreeMatrix(dataPoints points, matrix *outputMatrix);
int normalizedSimilarityMatrix(dataPoints points, matrix *resMatrix);
int transpose(matrix mat, matrix *transposedMatrix);
double trace(matrix mat);
int substractMatrices(matrix A, matrix B, matrix *result);
int updateH(matrix H, matrix W, matrix *updatedH);
int iterateAlgorithm(matrix *H, matrix W);
int calculateNumOfPoints(char *filename);
int calculateDimension(char *filename);
void printDataPoints(dataPoints points);