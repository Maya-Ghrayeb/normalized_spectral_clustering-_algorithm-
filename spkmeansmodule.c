#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include "spkmeans.h"

typedef struct
{
    double *centroid; /* the centroid of this specific cluster */
    int numberOfPointsInThisCluster;
    double *sumOfAllPointsInCluster; /* the sum of all the points in this specific cluster */
} cluster;


PyObject* cMatToPyMat(double** mat, int row, int cols){
    PyObject *res, *matRow;
    int i, j;
    res =  PyList_New(row);
    if (res == NULL)
        return NULL;
    
    for (i = 0; i < row; i++)
    {
        matRow = PyList_New(cols);
        if (matRow == NULL)
            return NULL;
        for (j = 0; j < cols; j++)
        {
            PyList_SET_ITEM(matRow, j, Py_BuildValue("d", mat[i][j]));
        }
        PyList_SetItem(res, i, Py_BuildValue("O", matRow));
    }
    return res;
}

double* pyObjToC1DArr(PyObject *pythonVectors, int vecCnt, int vecDim){
    int i, j, casting_failed;
    double *points;
    PyObject *vector, *coordinate;
    casting_failed = 0;
    points = (double*) calloc(vecCnt * vecDim, sizeof(double));
    /* fill the points */
    for (i = 0; i < vecCnt; i++) {
        vector = PyList_GetItem(pythonVectors, i);
        if (!PyList_Check(vector)){
            continue;
        }
        for (j = 0; j < vecDim; j++) {
            coordinate = PyList_GetItem(vector, j);
            points[i*vecDim + j] = PyFloat_AsDouble(coordinate); /** convert a Python float to C double **/
            /** check casting from python float to C double **/
            if (PyErr_Occurred() && points[i*vecDim + j]  == -1.0){
                casting_failed = 1;
                return NULL;
            }
        }
    }
    return points;

}

PyObject* jacobiOutputToPy3DMat(double*** mat, int row, int cols)
{
    PyObject *res, *mat1, *mat2;
    res = PyList_New(2);
    mat1 = cMatToPyMat(mat[0], row, cols);
    mat2 = cMatToPyMat(mat[1], row, cols);
    PyList_SetItem(res, 0, Py_BuildValue("O", mat1));
    PyList_SetItem(res, 1, Py_BuildValue("O", mat2));
    return res;
}

static PyObject* calcTMarix(PyObject* self, PyObject *args)
{
    int vecCnt, vecDim, casting_failed, numOfClusters;
    int numOfClustersInp;
    double *points;
    double * eigenValuesArr;    
    double** waj;
    double** ddg;
    double** lnorm;
    double ** eigenValuesMAt;
    double ** eigenVectorsMAT;
    double ** tMatrix;
    double*** jacobi_alg_output;
    eigenPair* eigenPairs;
   

    PyObject *pythonVectors, *res;
    casting_failed = 0;

    if (!PyArg_ParseTuple(args, "Oiii", &pythonVectors, &vecCnt, &vecDim, &numOfClustersInp))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    points = pyObjToC1DArr(pythonVectors, vecCnt, vecDim);

    waj = waj_calc(points, vecCnt, vecDim);
    ddg = ddg_calc(waj, vecCnt, vecCnt);
    lnorm = lnorm_calc(ddg, waj, vecCnt);
    jacobi_alg_output = jacobi_AlgVer2(lnorm, vecCnt);
    eigenValuesMAt = jacobi_alg_output[0];
    eigenVectorsMAT = jacobi_alg_output[1];
    eigenValuesArr = getEigenValues(eigenValuesMAt, vecCnt);
    
    
    eigenPairs = buildEigenPairs(eigenValuesArr, eigenVectorsMAT, vecCnt, vecCnt);
    qsort(eigenPairs, vecCnt, sizeof(eigenPair), md_comparator);
    
    if (numOfClustersInp == -1)
    {
        numOfClusters = getNumOfClusters(eigenPairs, vecCnt);
    }
    else
    {
        numOfClusters = numOfClustersInp;
    }
    tMatrix = getTMat(eigenPairs, vecCnt, numOfClusters);
    res = cMatToPyMat(tMatrix, vecCnt, numOfClusters);

    free(points);
    free_2dMat(waj, vecCnt);
    free_2dMat(ddg, vecCnt);
    free_2dMat(lnorm, vecCnt);
    free(jacobi_alg_output);
    free_2dMat(eigenValuesMAt, vecCnt);
    free_2dMat(eigenVectorsMAT, vecCnt);  
    free(eigenValuesArr);
    free_2dMat(tMatrix, vecCnt);

    for (int i = 0; i < vecCnt; i++){
        free(eigenPairs[i].eigenVector);
    }
    free(eigenPairs);
    //printf("Iam here");
    return res;

}

static double distance_calculator(double *v1, double *v2, int d)
{
    double distance;
    int c;
    double coordinate;
    distance = (double) 0;
    for(c=0; c<d ;c++)
    {
        coordinate = (double) ((double)(v1[c] - v2[c]))*((double)(v1[c] - v2[c]));
        distance += coordinate;
    }
    return distance;
}

static int findClusterOfClosestCentroid(double* point, cluster* clusters, int k, int vic_size) {
    double minVal = DBL_MAX;    /* defined in <float.h> */
    int clusterOfMinVal = -1;
    int j;
    double *centroidOfClusterJ;
    double dist;
    /* loop over all clusters: */
    for (j = 0 ; j < k ; j++) {
        centroidOfClusterJ = clusters[j].centroid;
        dist = distance_calculator(centroidOfClusterJ, point, vic_size);
        if (dist < minVal) {
            minVal = dist;
            clusterOfMinVal = j;
        }
    }
    return clusterOfMinVal;
}

static int xfindCluster(cluster* clusters, double* points, int k, int N, int vic_size) {
    int i;
    int j;
    double *point;
    int closestCluster;
    point = (double*) calloc (vic_size, sizeof(double));
    /* points is an array of N doubles - all points in the input file
    clusters is an array of k clusters - each cluster is of type cluster
     must first zero numberOfPointsInThisCluster and sumOfAllPointsInCluster in each cluster: */
    if (point == NULL){
        printf("An Error Has Occurred");
        return 0;
    }
    for (i = 0; i < k; i++) {
        for (j = 0; j < vic_size; j++)
        {
            clusters[i].sumOfAllPointsInCluster[j] = 0.0 ;
        }
        clusters[i].numberOfPointsInThisCluster = 0;
    }

    for (i=0; i < N; i++){
        /* warning your next bug may be this */
        for (j = 0; j < vic_size; j++){
            point[j] = points[i*vic_size + j];
        }

        /* find cluster of centroid that is closest to point */
        closestCluster = findClusterOfClosestCentroid(point, clusters, k, vic_size);
        /* "add" point to cluster closestCluster - actually we don't need to add it,
         just add it to the sum, and increase number of points in this cluster: */
        clusters[closestCluster].numberOfPointsInThisCluster += 1 ;
        for (j = 0; j < vic_size; j++){
            clusters[closestCluster].sumOfAllPointsInCluster[j] += point[j] ;
        }

    }
    free(point);
    return 1;
}

static void updateCentroids(cluster* clusters, int k, int vic_size) {
    int i;
    int j;
    for (i=0; i < k; i++){
        for (j=0 ; j < vic_size; j++)
        {
            clusters[i].centroid[j] = (clusters[i].sumOfAllPointsInCluster[j] / clusters[i].numberOfPointsInThisCluster);
        }
    }
}

static int Conver(cluster *clusters, double * prev_centroids, int k, int vic_len, double epsilon){
    int i;
    int j;
    double val;
    for (i = 0; i < k; i++){
        for (j = 0; j < vic_len; j ++){
            val = fabs(prev_centroids[i * vic_len + j] - clusters[i].centroid[j]);
            if (val > epsilon){
                return 0;
            }
        }
    }
/* val <epsilon --> 1= there's convergence */
    return 1;
}

static void free_2D_double_array(double **arr, int len){
    int i;
    for (i=0;i<len;i++){
        free(arr[i]);
    }
}

static double **kmeanspp(double *points, double **centroids, int line_cnt, int vic_size, int k, int maxiter, double epsilon)
{
    
    int Conv;
    int ox;
    int i;
    int j;
    int iter_num;
    double *oldPoints;
    cluster *clusters;
    Conv = 0;
    iter_num = 0;
    int m;
    oldPoints = (double*) calloc(k * vic_size, sizeof(double));
    if(oldPoints==NULL){
        printf("An Error Has Occurred");
        return NULL;
    }

    clusters = (cluster*) calloc(k, sizeof(cluster));
    if(clusters==NULL){
        printf("An Error Has Occurred");
        return NULL;
    }

    for (i = 0; i < k; i++) {
        clusters[i].centroid = (double*) calloc(vic_size, sizeof(double));
        if(clusters[i].centroid==NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
        clusters[i].sumOfAllPointsInCluster = (double*) calloc(vic_size, sizeof(double));
        if(clusters[i].sumOfAllPointsInCluster==NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
        for (j = 0; j < vic_size; j ++){
            clusters[i].centroid[j] = centroids[i][j]; /* initialize the centroid to the first k data points */
            clusters[i].sumOfAllPointsInCluster[j] = centroids[i][j];
        }
        clusters[i].numberOfPointsInThisCluster = 1;
    }


    while (iter_num < maxiter && Conv == 0)
    {
        ox = xfindCluster(clusters, points, k, line_cnt , vic_size);
        if (ox == 0){
            return NULL;
        }

        for (i = 0; i < k; i ++){
            for (j = 0; j < vic_size; j++){
                oldPoints[i * vic_size + j] = clusters[i].centroid[j];
            }
        }

        updateCentroids(clusters, k, vic_size);

        Conv = Conver(clusters, oldPoints, k, vic_size, epsilon);
       
        iter_num++;
    }

    for(i=0; i < k; i++)
    {
        for(j = 0; j < vic_size; j++)
        {
            centroids[i][j]=clusters[i].centroid[j];
        }
    }
    for(m=0;m<k;m++){
        free(clusters[m].centroid);
        free(clusters[m].sumOfAllPointsInCluster);
    }
    free(clusters);
    free(oldPoints);

    return centroids;
}

double** pyMatToCMAt(PyObject *pythonMat, int rowsCnt, int colsCnt)
{
    int i, j;
    double** res;
    PyObject *coordinate, *rowMat;
    /* need to replaced with Create2dARR*/
    res = create_2d_mat(rowsCnt, colsCnt);
    if (res == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;

    }
    
    for (i = 0; i < rowsCnt; i++) 
    {
        rowMat = PyList_GetItem(pythonMat, i);
        if (!PyList_Check(rowMat)){
            continue;
        }
        for (j = 0; j < colsCnt; j++) {
            coordinate = PyList_GetItem(rowMat, j);
            res[i][j] = PyFloat_AsDouble(coordinate); /** convert a Python float to C double **/
            /** check casting from python float to C double **/
            if (PyErr_Occurred() && res[i][j]  == -1.0){
                return NULL;
            }
        }
    }
    return res;
}

static PyObject* kmeansppCalc(PyObject* self, PyObject *args)
{
    int numOfClusters, vecCnt, vecDim, max_iter;
    double epsilon;
    double *points, **centroids, **final_centroids;
    int casting_failed;
    PyObject *pythonCentroids, *pythonVectors, *pythonFinalCentroids;
    casting_failed = 0;

    /** validate we got all parameters **/
    if (!PyArg_ParseTuple(args, "OOiiiid", &pythonCentroids, &pythonVectors, &max_iter, &vecCnt, &vecDim, &numOfClusters, &epsilon))
        return NULL;
    if (!PyList_Check(pythonCentroids))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    /** fill the points**/
    points = pyObjToC1DArr(pythonVectors, vecCnt, vecDim);
    
    /** fill the centroids**/
    centroids = pyMatToCMAt(pythonCentroids, numOfClusters, vecDim);

    if(points == NULL || centroids == NULL)
    {
        casting_failed = 1;
    }

    final_centroids = NULL;

    if (casting_failed == 0){
        final_centroids = kmeanspp(points, centroids, vecCnt, vecDim, numOfClusters, max_iter, epsilon);
        if (final_centroids == NULL)
        {
            return NULL;
        }
    }

    pythonFinalCentroids = cMatToPyMat(final_centroids, numOfClusters, vecDim);
    free(points);
    free_2D_double_array(centroids, numOfClusters);
    free(centroids);
    
    return pythonFinalCentroids;
}

static PyObject* calcWAM(PyObject* self, PyObject *args)
{
    int vecCnt, vecDim, casting_failed;
    PyObject *pythonVectors, *res;
    double** wam;
    double* points;
    casting_failed = 0;
    if (!PyArg_ParseTuple(args, "Oii", &pythonVectors, &vecCnt, &vecDim))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;
    
    points = pyObjToC1DArr(pythonVectors, vecCnt, vecDim);
    wam = waj_calc(points, vecCnt, vecDim);
    res = cMatToPyMat(wam, vecCnt, vecCnt);
    free(points);
    free_2dMat(wam, vecCnt);
    return res;
}

static PyObject* calcDDG(PyObject* self, PyObject *args)
{
    int vecCnt, vecDim, casting_failed;
    PyObject *pythonVectors, *res;
    double** input;
    double** ddg;
    casting_failed = 0;
    if (!PyArg_ParseTuple(args, "Oii", &pythonVectors, &vecCnt, &vecDim))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;
    input = pyMatToCMAt(pythonVectors, vecCnt, vecDim);
    ddg = ddg_calc(input, vecCnt, vecDim);
    res = cMatToPyMat(ddg, vecCnt, vecDim);
    free_2dMat(input, vecCnt);
    free_2dMat(ddg, vecCnt);
    return res;
}

static PyObject* calcLnorm(PyObject* self, PyObject *args)
{
    int vecCnt, vecDim, casting_failed;
    double *points;
    double** waj;
    double** ddg;
    double** lnorm;

    PyObject *pythonVectors, *res;
    casting_failed = 0;

    if (!PyArg_ParseTuple(args, "Oii", &pythonVectors, &vecCnt, &vecDim))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;

    points = pyObjToC1DArr(pythonVectors, vecCnt, vecDim);
    waj = waj_calc(points, vecCnt, vecDim);
    ddg = ddg_calc(waj, vecCnt, vecCnt);
    lnorm = lnorm_calc(ddg, waj, vecCnt);
    res = cMatToPyMat(lnorm, vecCnt, vecCnt);
    free(points);
    free_2dMat(waj, vecCnt);
    free_2dMat(ddg, vecCnt);
    free_2dMat(lnorm, vecCnt);
    return res;
}

static PyObject* calcJacobi(PyObject* self, PyObject *args)
{
    int vecCnt, vecDim, casting_failed;
    PyObject *pythonVectors, *res;
    double** input;
    double ** eigenValuesMAt;
    double ** eigenVectorsMAT;
    double*** jacobiAlgOutput;
    casting_failed = 0;
    if (!PyArg_ParseTuple(args, "Oii", &pythonVectors, &vecCnt, &vecDim))
        return NULL;
    if (!PyList_Check(pythonVectors))
        return NULL;
    input = pyMatToCMAt(pythonVectors, vecCnt, vecDim);
    jacobiAlgOutput = jacobi_AlgVer2(input, vecCnt);
    eigenValuesMAt = jacobiAlgOutput[0];
    eigenVectorsMAT = jacobiAlgOutput[1];
    res = jacobiOutputToPy3DMat(jacobiAlgOutput, vecCnt, vecCnt);
    free(jacobiAlgOutput);
    free_2dMat(eigenValuesMAt, vecCnt);
    free_2dMat(eigenVectorsMAT, vecCnt);
    return res; 
}

static PyMethodDef spkmeansMethods[] = {
        {"kmeansppCalc", (PyCFunction) kmeansppCalc, METH_VARARGS, PyDoc_STR("C kmeanspp algorithm")},
        {"calcTMatrix", (PyCFunction) calcTMarix, METH_VARARGS, PyDoc_STR("C spkmeans algorithm")},
        {"calcWAM", (PyCFunction) calcWAM, METH_VARARGS, PyDoc_STR("C spkmeans algorithm")},
        {"calcDDG", (PyCFunction) calcDDG, METH_VARARGS, PyDoc_STR("C spkmeans algorithm")},
        {"calcLnorm", (PyCFunction) calcLnorm, METH_VARARGS, PyDoc_STR("C spkmeans algorithm")},
        {"calcJacobi", (PyCFunction) calcJacobi, METH_VARARGS, PyDoc_STR("C spkmeans algorithm")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT, "myspkmeans", NULL, -1, spkmeansMethods
};

PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}