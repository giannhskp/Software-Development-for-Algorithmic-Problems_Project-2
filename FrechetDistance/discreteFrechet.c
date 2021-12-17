#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../BinaryTree/binaryTree.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MAX3(a, b, c) MAX((MAX((a),(b))),(c))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MIN3(a, b, c) MIN((MIN((a),(b))),(c))

// structure that used to save the optimal path indexes
struct index_node {
  int i;
  int j;
};
typedef struct index_node Index;

static double l2_metric(double x1,double y1,double x2,double y2){
  // calculate the Euclidean distance (or L2) between the given vectors and return it
  double temp = sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1)));
  return temp;
}

double discreteFrechet(Vector v1,Vector v2){
  // calculate the Discrete Frechet distance between the 2 given timeseries.
  // these two timeseries may have different dimensions, so the dynamic 2-d array
  // has as many rows as the first given timeseries and as many columns as the second timeseries has.
  // We find the Discrete Frechet distance between them by applying Dynamic Programming,
  // the requested distance will be appeared in the lower right element of the dynamic array.

  int i_dim = getDim(v1);
  int j_dim = getDim(v2);
  double **dynamicArray=malloc(i_dim*sizeof(double*));
  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  double *coordsTime1 = getTime(v1);
  double *coordsTime2 = getTime(v2);


  for(int i=0;i<i_dim;i++){
    dynamicArray[i]=malloc(j_dim*sizeof(double));
  }

  // initialize the dynamic array
  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime1[0],coords2[0],coordsTime2[0]);

  // for the first row of the dynamic array -> i=0, j>0
  for(int j=1;j<j_dim;j++){
    double p1qj = l2_metric(coords1[0],coordsTime1[0],coords2[j],coordsTime2[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
  }

  // for the first column of the dynamic array -> i>0, j=0
  for(int i=1;i<i_dim;i++){
    double piq1 = l2_metric(coords1[i],coordsTime1[i],coords2[0],coordsTime2[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
  }

  // for the remaining elements of the dynamic array - > i>0, j>0
  for(int i=1;i<i_dim;i++){
    for(int j=1;j<j_dim;j++){
      double piqj = l2_metric(coords1[i],coordsTime1[i],coords2[j],coordsTime2[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }

  // Discrete Distance has been calculated
  double finalDistance = dynamicArray[i_dim-1][j_dim-1];

  for(int i=0;i<i_dim;i++){
    free(dynamicArray[i]);
  }
  free(dynamicArray);
  return finalDistance;
}

Index *discreteFrechet_optimalPath(Vector v1,Vector v2,int *pathLength){
  // find the optimal path between the 2 given timeseries that comes up through the calculation of Discrete Frechet distance.
  // these two timeseries may have different dimensions, so the dynamic 2-d array
  // has as many rows as the first given timeseries and as many columns as the second timeseries has.
  // Αs the dynamic array is completed, the backtracking path is completed by the same time.
  // Therefore, function returns the optimal path that has been calculated with backtracking.

  int i_dim = getDim(v1);
  int j_dim = getDim(v2);


  double **dynamicArray=malloc(i_dim*sizeof(double*));

  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  double *coordsTime1 = getTime(v1);
  double *coordsTime2 = getTime(v2);

  Index **backtrackingPath=malloc(i_dim*sizeof(Index*));
  Index *optimalPath = malloc((i_dim+j_dim)*sizeof(Index));


  for(int i=0;i<i_dim;i++){
    dynamicArray[i]=malloc(j_dim*sizeof(double));
    backtrackingPath[i]=malloc(j_dim*sizeof(Index));
  }

  // initialize the dynamic array
  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime1[0],coords2[0],coordsTime2[0]);
  // initialize the backtracking array
  backtrackingPath[0][0].i=-1;
  backtrackingPath[0][0].j=-1;

  // for the first row of the dynamic array -> i=0, j>0
  for(int j=1;j<j_dim;j++){
    double p1qj = l2_metric(coords1[0],coordsTime1[0],coords2[j],coordsTime2[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
    // and the corresponding indexes saved at the backtracking array
    backtrackingPath[0][j].i=0;
    backtrackingPath[0][j].j=j-1;
  }

  // for the first column of the dynamic array -> i>0, j=0
  for(int i=1;i<i_dim;i++){
    double piq1 = l2_metric(coords1[i],coordsTime1[i],coords2[0],coordsTime2[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
    // and the corresponding indexes saved at the backtracking array
    backtrackingPath[i][0].i=i-1;
    backtrackingPath[i][0].j=0;
  }

  // for the remaining elements of the dynamic array - > i>0, j>0
  for(int i=1;i<i_dim;i++){
    for(int j=1;j<j_dim;j++){
      double piqj = l2_metric(coords1[i],coordsTime1[i],coords2[j],coordsTime2[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      // for every dynamic array's element find the indexes of the point that optimal path came up
      // so find the minimum element among the three 3 above [i-1][j] , [i-1][j-1] , [i][j-1]
      // and save the corresponding indexes at the backtracking array
      if(minC==dynamicArray[i-1][j]){
        backtrackingPath[i][j].i=i-1;
        backtrackingPath[i][j].j=j;
      }else if(minC==dynamicArray[i-1][j-1]){
        backtrackingPath[i][j].i=i-1;
        backtrackingPath[i][j].j=j-1;
      }else{
        backtrackingPath[i][j].i=i;
        backtrackingPath[i][j].j=j-1;
      }
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }
  // initialize the optimal path
  int optimalPathLenth = 0;
  optimalPath[optimalPathLenth].i = i_dim-1;
  optimalPath[optimalPathLenth++].j = j_dim-1;


  // let's form the optimal path with backtracking
  // we begin from the lower right element of the backtracking array [i-1][j-1]
  // and we are following the correspodings indexes until [0][0] element reached.
  // This path from [i-1][j-1] to [0][0] is the optimal path
  Index backtrackIndex = backtrackingPath[i_dim-1][j_dim-1];
  optimalPath[optimalPathLenth++] = backtrackIndex;
  while(1){

    backtrackIndex = backtrackingPath[backtrackIndex.i][backtrackIndex.j];

    if(optimalPathLenth>=(i_dim+j_dim)){  // never actually happens
      break;
    }
    if(backtrackIndex.i==-1 && backtrackIndex.j==-1){
      break;
    }
    optimalPath[optimalPathLenth++] = backtrackIndex;
  }

  (*pathLength) = optimalPathLenth;

  for(int i=0;i<i_dim;i++){
    free(dynamicArray[i]);
    free(backtrackingPath[i]);
  }
  free(dynamicArray);
  free(backtrackingPath);

  // finally return the optimal path
  return optimalPath;
}


Vector meanCurveBetween2Curves(Vector v1,Vector v2){
  // find and return the mean curve between the two given timeseries/curves.
  // more especially to find the mean curve, optimal path will be calculated
  // from the function call of discreteFrechet_optimalPath(...)
  // optimal path is a sequence of indexes (i,j),
  // i corresponds to the first timeseries and j corresponds to the second one.
  // for every tuple (i,j) of the optimal path the correspoding point (x,y) of the mean curve
  // comes up by the following formula:
  // (x,y) = ( (vector1_x[i]+vector2_x[j])/2 , (vector1_y[i]+vector2_y[j])/2)
  // thus, the mean curve is calculated and finally formed.

  int pathLength = 0;
  Index *optimalPath = discreteFrechet_optimalPath(v1,v2,&pathLength);
  double *meanCurveCoord = malloc(pathLength*sizeof(double));
  double *meanCurveTime = malloc(pathLength*sizeof(double));
  for(int i=0;i<pathLength;i++){
    Index index = optimalPath[i];
    double t1 = getCoords(v1)[index.i];
    double t2 = getCoords(v2)[index.j];
    double meanCoord = (t1+t2)/2;
    double t3 = getTime(v1)[index.i];
    double t4 = getTime(v2)[index.j];
    double meanTime = (t3+t4)/2;
    meanCurveCoord[pathLength-1-i] = meanCoord;
    meanCurveTime[pathLength-1-i] = meanTime;
  }
  Vector meanCurve = initTimeSeries(meanCurveCoord,meanCurveTime,"meanCurve",pathLength);
  free(meanCurveCoord);
  free(meanCurveTime);
  free(optimalPath);
  return meanCurve;
}

Vector computeFrechetMeanCurve(List list,int count){
  Tree tree = createTreeFromList(list,count);
  Vector newCenter = treeFindMeanCurve(tree);
  destroyTree(tree);
  return newCenter;
}

Vector computeFrechetMeanCurveLSH(HashTable ht,int count){
  Tree tree = createTreeFromHt(ht,count);
  Vector newCenter = treeFindMeanCurve(tree);
  destroyTree(tree);
  return newCenter;
}
