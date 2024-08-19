#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>

typedef struct
{
    double eigenValue;
    double *eigenVector;
} eigenPair;

double** waj_calc(double *orgMat, int vicCnt, int vicDim);

double** ddg_calc(double** waj, int vicCnt, int colsCnt);

double** lnorm_calc(double ** ddg,double** waj ,  int vicCnt);

double*** jacobi_AlgVer2(double** lnorm, int vicCnt);

double * getEigenValues(double ** eigenValuesMat, int vicCnt);

eigenPair* buildEigenPairs(double * eigenValuesArr, double** eigenVectors, int arrLen, int clusterVecLen);

int getNumOfClusters(eigenPair* sortedEigenPairs, int arrLen);

double** getTMat(eigenPair* eigenPairs, int vicCnt, int numOfClusters);

void print2DMAt(double ** mat, int rows, int cols);

int md_comparator(const void *v1, const void *v2);

void free_2dMat(double** mat, int rows);

double ** create_2d_mat(int rows, int cols);