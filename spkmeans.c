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

enum goal
{
    wam,
    ddg,
    lnorm,
    jacobi
};

void print2DMAt(double ** mat, int rows, int cols)
{
    int i, j;
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            if (j !=  cols - 1){
                printf("%.4f,", mat[i][j]);
            }
            else{
                printf("%.4f", mat[i][j]);

            }
        }
        printf("\n");
    }
}

void printorgMat(double * orgMat, int vic_diminsion, int vicCnt){
    int i;
    for(i = 0; i < vic_diminsion*vicCnt; i++){
        if ((i % vic_diminsion) == 0){
            printf("\n");
        }
        printf("%f ,", orgMat[i]);
    }
}

double euclidean_distance_calculator(double *v1, double *v2, int diminsion)
{
    double distance;
    int c;
    double coordinate;
    distance = 0;
    for(c = 0; c < diminsion ;c++)
    {
        coordinate = pow((v1[c] - v2[c]), 2);
        distance += coordinate;
    }
    return distance;
}

void fillMat(FILE *mainFile, double* orgMat){
    double vic = 0.0;
    int i;
    i = 0;
    while (fscanf(mainFile, "%lf,", &vic) != EOF)
    {
        orgMat[i] = vic;
        i++;
    }

}

int getVicDim(FILE* mainFile){
    char character;
    int vic_size;
    vic_size = 0;

    while((character = getc(mainFile)) !=  '\n' )
    {
        if(character == ','){
            vic_size++;
        }
    }
    rewind(mainFile);
    vic_size ++; /* inOrder to have the right value */

    return vic_size;
}

int getLineCnt(FILE* mainFile){
    int lineCnt;
    char character;

    lineCnt = 0;
    while( (character = getc(mainFile)) !=  EOF )
    {
        if(character == '\n'){
            lineCnt++;
        }
    }
    rewind(mainFile);
    return lineCnt;
}

double *fileToArr(FILE *mainFile, int vic_diminsion, int vicCnt){
    double *orgMat;
    orgMat = (double*) calloc(vic_diminsion * (vicCnt), sizeof(double));
    fillMat(mainFile, orgMat);
    return orgMat;
}

double ** create_2d_mat(int rows, int cols){
    double ** mat;
    int i;
    mat = (double**) calloc(rows, sizeof(double*));
    if (mat == NULL){
        return NULL;
    }
    for(i = 0; i < rows; i++){
        mat[i] = (double*) calloc(cols, sizeof(double));
        if (mat[i] == NULL)
        {
            return NULL;
        }
    }
    return mat;
}

void free_2dMat(double** mat, int rows){
    int i;
    for(i = 0; i < rows; i++){
        free(mat[i]);
    }
    free(mat);
}

double** squareIdentityMat(int vicCnt){
    double** identity_mat;
    int i;
    identity_mat = create_2d_mat(vicCnt, vicCnt);
    for (i = 0; i < vicCnt; i++){
        identity_mat[i][i] = 1;
    }
    return identity_mat;
}

double** transposeMat(double** mat, int rows, int cols){
    double** res;
    int i;
    int j;
    res = create_2d_mat(cols, rows);
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            res[j][i] = mat[i][j];
        }
    }
    return res;
}

double sign(double x){
    if (x>= 0){
        return 1;
    }
    return -1;
}

int* getIndexesForPMatrix(double** mat, int vicCnt){
    int* resIndexes;
    int index_i;
    int index_j;
    double tmp;
    int i;
    int j;
    resIndexes = (int*) calloc(2, sizeof(int));
    tmp = 0;
    index_i = 0;
    index_j = 1;
    for(i = 0; i < vicCnt; i++){
        for(j = i; j < vicCnt; j++){
            if ((i != j) && (fabs(mat[i][j]) > fabs(tmp))){
                tmp = mat[i][j];
                index_i = i;
                index_j = j;
            }
        }
    }
    resIndexes[0] = index_i;
    resIndexes[1] = index_j;
    return resIndexes;
}

/*code neeed to be improved*/
double** rotationMat(double** mat, int vicCnt)
{
    double** res;
    int * valIndexes;
    int index_i;
    int index_j;
    double realVal;
    double theta;
    double t;
    double c;
    double s;

    valIndexes = getIndexesForPMatrix(mat, vicCnt);
    index_i = valIndexes[0];
    index_j = valIndexes[1];
    realVal = mat[index_i][index_j];

    theta = (mat[index_j][index_j] - mat[index_i][index_i])/ (2*realVal);
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2)+1)); /*there was abug in here*/
    c = 1/(sqrt(pow(t, 2) + 1));
    s = t* c;

    res = squareIdentityMat(vicCnt);

    res[index_i][index_j] = s;
    res[index_j][index_i] = -1 * s;
    res[index_j][index_j] = c;
    res[index_i][index_i] = c;
    free(valIndexes);
    return res;
}
/*
weighted_Adjacency_Matrix Calculator
*/
double** waj_calc(double *orgMat, int vicCnt, int vicDim){
    double ** waj;
    int i;
    int j;
    double val;
    double distance;

    waj = create_2d_mat(vicCnt, vicCnt);

    for (i = 0; i < vicCnt; ++i){
        for(j = 0; j < vicCnt; ++j){
            if (i == j){
                waj[i][j] = 0;
            }
            else{
                distance = euclidean_distance_calculator(orgMat + (i * vicDim), orgMat + (j * vicDim), vicDim);
                val = -1* (sqrt(distance) / 2);
                waj[i][j] = exp(val);
            }
        }
    }
    return waj;
}

double** matrixMulti(double** mat1, double** mat2,  int vicCnt){
    double** res;
    int i;
    int j;
    int k;
    double sum;

    res = create_2d_mat(vicCnt, vicCnt);

    for (i = 0; i < vicCnt; i++){
        for (j = 0; j < vicCnt; j++){
            sum = 0;
            for (k = 0; k < vicCnt; k++){
                sum += mat1[i][k] * mat2[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

double** ddg_calc(double** waj, int vicCnt, int colsCnt){
    double ** ddg;
    int i;
    int j;
    double d_i;
    ddg = create_2d_mat(vicCnt, colsCnt);
    for (i = 0; i < vicCnt; i++){
        d_i = 0;
        for (j = 0; j < colsCnt; j++){
            d_i += waj[i][j];
        }
        ddg[i][i] = d_i;
    }
    return ddg;
}

void squareFilter(double** ddg, int vecCnt)
{
    int i;
    double d_i;
    for (i = 0; i < vecCnt; i++)
    {
        d_i = ddg[i][i];
        ddg[i][i] = pow(d_i, -0.5);
    }
}

double** matrix_substract(double** mat1, double** mat2, int vicCnt){
    double** res;
    int i;
    int j;

    res = create_2d_mat(vicCnt, vicCnt);

    for(i = 0; i < vicCnt; i++ ){
        for(j = 0; j < vicCnt; j++){
            res[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return res;
}

double** lnorm_calc(double ** ddg, double** waj,  int vicCnt){
    double ** identity_mat;
    double ** res1;
    double ** res2;
    double ** lnorm;

    squareFilter(ddg, vicCnt);
    /* build Identity matrix */
    identity_mat = squareIdentityMat(vicCnt);

    /* matrix multiplication */
    res1 = matrixMulti(ddg, waj, vicCnt);
    res2 = matrixMulti(res1, ddg, vicCnt);

    /* matrix substraction */
    lnorm = matrix_substract(identity_mat, res2, vicCnt);

    /* free the memory*/
    free_2dMat(res1, vicCnt);
    free_2dMat(res2, vicCnt);
    free_2dMat(identity_mat, vicCnt);

    return lnorm;
}

double off(double** mat, int rows, int cols){
    double sum;
    int i;
    int j;
    sum = 0;
    for(i = 0; i < rows; i++){
        for(j = i + 1; j < cols; j++){
            if (i != j){
                sum += pow(mat[i][j], 2);
            }
        }
    }
    return 2 * sum;
}

int converge(double ** mat1, double ** mat2, double epsilon, int vicCnt){
    if ((off(mat1, vicCnt, vicCnt) - off(mat2, vicCnt, vicCnt)) <= epsilon){
        return 1;
    }
    return 0;
}

void rotateP(double ** mat, int index_i, int index_j){
    double tmp;
    tmp = mat[index_i][index_j];
    mat[index_i][index_j] = mat[index_j][index_i];
    mat[index_i][index_j] = tmp;
}

double*** jacobi_AlgVer2(double** lnorm, int vicCnt){
    double *** res;
    double ** rotationMa_P;
    double ** rotationMa_P_T;
    double ** firstMultiRes;
    double ** secondMultiRes;
    double ** v;
    double ** tmp_matrix;
    double ** tmp_v;
    int i;
    int max_iter;
    int conv;
    double epsilon;
    res = (double ***) calloc(2, sizeof(double **));
    i = 0;
    max_iter = 100;
    epsilon = 1.0 * pow(10,-5);
    conv = 0;
    tmp_matrix = lnorm;
    v = squareIdentityMat(vicCnt);
    /*P matrix intiazation */

    while (i < max_iter && (conv != 1))
    {
        /*intilize */
        rotationMa_P = rotationMat(tmp_matrix, vicCnt);
        rotationMa_P_T = transposeMat(rotationMa_P, vicCnt, vicCnt);

        firstMultiRes = matrixMulti(rotationMa_P_T, tmp_matrix, vicCnt);
        secondMultiRes = matrixMulti(firstMultiRes, rotationMa_P, vicCnt);
        tmp_v = v;
        v = matrixMulti(tmp_v, rotationMa_P, vicCnt);

        if (i > 0){ free_2dMat(tmp_matrix, vicCnt); }
        free_2dMat(tmp_v, vicCnt);
        free_2dMat(firstMultiRes, vicCnt);
        free_2dMat(rotationMa_P, vicCnt);
        free_2dMat(rotationMa_P_T, vicCnt);


        conv = converge(tmp_matrix, secondMultiRes, epsilon, vicCnt);


        tmp_matrix = secondMultiRes;
        i++;
    }
    res[0] = tmp_matrix;
    res[1] = v;
    return res;
}

double * getEigenValues(double ** eigenValuesMat, int vicCnt){
    double* res;
    int i;

    res = (double*) calloc(vicCnt, sizeof(double));
    for (i = 0; i < vicCnt; i++){
        res[i] = eigenValuesMat[i][i];
    }
    return res;
}

int md_comparator(const void *v1, const void *v2)
{
    const eigenPair *p1 = (eigenPair *)v1;
    const eigenPair *p2 = (eigenPair *)v2;
    if (p1-> eigenValue< p2->eigenValue)
        return -1;
    else if (p1->eigenValue > p2->eigenValue)
        return +1;
    else
        return 0;
}

int getNumOfClusters(eigenPair* sortedEigenPairs, int arrLen){
    double res;
    int i;
    int index;
    res = -1;
    index = -1;
    for (i = 0; i < floor((arrLen/2)); i++){
        double val = sortedEigenPairs[i + 1].eigenValue - sortedEigenPairs[i].eigenValue;
        if (val > res){
            res = val;
            index = i;
        }
    }
    return index + 1;
}

eigenPair* buildEigenPairs(double * eigenValuesArr, double** eigenVectors, int arrLen, int clusterVecLen){
    eigenPair* res;
    int j;
    int i;
    res = (eigenPair*) calloc(arrLen, sizeof(eigenPair));
    for (i = 0; i < arrLen; i++){
        res[i].eigenVector = (double*) calloc(clusterVecLen, sizeof(double));
        res[i].eigenValue = eigenValuesArr[i];
        for (j = 0; j < clusterVecLen; j++){
            res[i].eigenVector[j] = eigenVectors[j][i];
        }
    }
    return res;
}

/* return U matrix containing getting the first k columns of V by eigen value order*/
double** getUMatrix(eigenPair* eigenPairs, int vicCnt, int numOfClusters){
    int i,j;
    double** uMat = create_2d_mat(vicCnt, numOfClusters);
    for(i = 0; i < vicCnt; i++){
        for(j = 0; j < numOfClusters; j++){
            uMat[i][j] = eigenPairs[j].eigenVector[i];
        }
    }
    return uMat;
}

/* calc the norm of a vector*/
double getRowNorm(double* rowVector,int vecDim){
    int i;
    double sum = 0;
    for(i = 0; i < vecDim; i++){
        sum += pow(rowVector[i], 2);
    }
    return sqrt(sum);
}

/* normalizing U matrix using getRowNorm function*/
void normalize(double** uMat, int vicCnt, int numOfClusters){
    int i, j;
    double rowNorm;
    for(i = 0; i < vicCnt; i++){
        rowNorm = getRowNorm(uMat[i], numOfClusters); /*changed*/
        if (rowNorm != 0) {
            for(j = 0; j < numOfClusters; j++){
                uMat[i][j] /= rowNorm;
            }
        }
    }
}

double** getTMat(eigenPair* eigenPairs, int vicCnt, int numOfClusters){
    double** tMat;
    tMat = getUMatrix(eigenPairs, vicCnt, numOfClusters);
    normalize(tMat, vicCnt, numOfClusters);
    return tMat;
}

int main(int argc, char** argv) {
    FILE *mainFile;
    double *orgMat;
    double **weighted_Adjacency_Matrix;
    double **diagonal_Degree_Matrix;
    double **normalized_Graph_Laplacian;
    double ***jacobi_alg_output;
    double **eigenValuesMAt;
    double **eigenVectorsMAT;
    double **tMatrix;
    double *eigenValuesArr;
    eigenPair *eigenPairs;
    int numOfClusters;
    char *fileName;
    int vic_diminsion;
    int vicCnt;
    int i;

    if (argc != 4) {
        printf("Invalid Input");
        return 9;
    }

    /* TODO: free up memory in every step in the process and not only in the end */
    /* file reading*/
    fileName = argv[1];
    mainFile = fopen(fileName, "r");
    vic_diminsion = getVicDim(mainFile);
    vicCnt = getLineCnt(mainFile);
    rewind(mainFile);
    orgMat = fileToArr(mainFile, vic_diminsion, vicCnt);
    fclose(mainFile);

    weighted_Adjacency_Matrix = waj_calc(orgMat, vicCnt, vic_diminsion);

    /* need to add another func */
    /* note here that it D^(-0.5)*/
    diagonal_Degree_Matrix = ddg_calc(weighted_Adjacency_Matrix, vicCnt, vicCnt);

    printf("\n");
    normalized_Graph_Laplacian = lnorm_calc(diagonal_Degree_Matrix, weighted_Adjacency_Matrix, vicCnt);

    /* geting the eigenValues and EigenVectors*/

    jacobi_alg_output = jacobi_AlgVer2(normalized_Graph_Laplacian, vicCnt);

    eigenValuesMAt = jacobi_alg_output[0];
    eigenVectorsMAT = jacobi_alg_output[1];

    eigenValuesArr = getEigenValues(eigenValuesMAt, vicCnt);


    eigenPairs = buildEigenPairs(eigenValuesArr, eigenVectorsMAT, vicCnt, vicCnt);
    qsort(eigenPairs, vicCnt, sizeof(eigenPair), md_comparator);

    numOfClusters = getNumOfClusters(eigenPairs, vicCnt);

    tMatrix = getTMat(eigenPairs, vicCnt, numOfClusters);
    print2DMAt(tMatrix, vicCnt, numOfClusters);

    free(orgMat);
    free_2dMat(weighted_Adjacency_Matrix, vicCnt);
    free_2dMat(diagonal_Degree_Matrix, vicCnt);
    free_2dMat(normalized_Graph_Laplacian, vicCnt);
    free(jacobi_alg_output);
    free_2dMat(eigenValuesMAt, vicCnt);
    free_2dMat(eigenVectorsMAT, vicCnt);
    free(eigenValuesArr);
    free_2dMat(tMatrix, vicCnt);

    for (i = 0; i < vicCnt; i++) {
        free(eigenPairs[i].eigenVector);
    }
    free(eigenPairs);
    printf("Iam here");
    return 0;
}