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

static double l2_metric(double x1,double y1,double x2,double y2){
  // calculate the Euclidean distance (or L2) between the given vectors and return it
  return sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1)));
}

double discreteFrechet(Vector v1,Vector v2){
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
  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime1[0],coords2[0],coordsTime2[0]);
  for(int j=1;j<j_dim;j++){
    double p1qj = l2_metric(coords1[0],coordsTime1[0],coords2[j],coordsTime2[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
  }

  for(int i=1;i<i_dim;i++){
    double piq1 = l2_metric(coords1[i],coordsTime1[i],coords2[0],coordsTime2[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
  }

  for(int i=1;i<i_dim;i++){
    for(int j=1;j<j_dim;j++){
      double piqj = l2_metric(coords1[i],coordsTime1[i],coords2[j],coordsTime2[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }

  double finalDistance = dynamicArray[i_dim-1][j_dim-1];

  for(int i=0;i<i_dim;i++){
    free(dynamicArray[i]);
  }
  free(dynamicArray);
  return finalDistance;
}

int *discreteFrechet_optimalPath(Vector v1,Vector v2,int *pathLength){
  int i_dim = getDim(v1);
  int j_dim = getDim(v2);
  double **dynamicArray=malloc(i_dim*sizeof(double*));
  int **backtrackingPath=malloc(i_dim*sizeof(int*));
  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  double *coordsTime1 = getTime(v1);
  double *coordsTime2 = getTime(v2);
  int *optimalPath = malloc((i_dim+j_dim)*sizeof(int));

  for(int i=0;i<i_dim;i++){
    dynamicArray[i]=malloc(j_dim*sizeof(double));
    backtrackingPath[i]=malloc(j_dim*sizeof(int));
  }

  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime1[0],coords2[0],coordsTime2[0]);
  backtrackingPath[0][0] = -1;
  for(int j=1;j<j_dim;j++){
    double p1qj = l2_metric(coords1[0],coordsTime1[0],coords2[j],coordsTime2[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
    backtrackingPath[0][j] =  j-1;
  }

  for(int i=1;i<i_dim;i++){
    double piq1 = l2_metric(coords1[i],coordsTime1[i],coords2[0],coordsTime2[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
    backtrackingPath[i][0] =  (i-1)*i_dim;
  }

  for(int i=1;i<i_dim;i++){
    for(int j=1;j<j_dim;j++){
      double piqj = l2_metric(coords1[i],coordsTime1[i],coords2[j],coordsTime2[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      if(minC==dynamicArray[i-1][j]){
        backtrackingPath[i][j] = (i-1)*i_dim+(j);
      }else if(minC==dynamicArray[i-1][j-1]){
        backtrackingPath[i][j] = (i-1)*i_dim+(j-1);
      }else{
        backtrackingPath[i][j] = (i)*i_dim+(j-1);
      }
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }
  int optimalPathLenth = 0;
  int backtrackIndex = backtrackingPath[i_dim-1][j_dim-1];
  optimalPath[optimalPathLenth++] = backtrackIndex;

  while(1){
    if(backtrackIndex==-1){
      optimalPath[optimalPathLenth++] = -1;
      break;
    }
    int x = backtrackIndex / i_dim;
    int y = backtrackIndex % i_dim;
    backtrackIndex = backtrackingPath[x][y];
    if(optimalPathLenth>=(i_dim+j_dim)){
      printf("!!! optimalPath OUT OF RANGE!\n");
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
  return optimalPath;
}

Vector meanCurveBetween2Curves(Vector v1,Vector v2){
  int pathLength = 0;
  int i_dim = getDim(v1);
  int *optimalPath = discreteFrechet_optimalPath(v1,v2,&pathLength);
  double *meanCurveCoord = malloc(pathLength*sizeof(double));
  double *meanCurveTime = malloc(pathLength*sizeof(double));
  for(int i=0;i<pathLength;i++){
    int index = optimalPath[i];
    int x = index / i_dim;
    int y = index % i_dim;
    double meanCoord = (getCoords(v1)[x] + getCoords(v2)[y])/2;
    double meanTime = (getTime(v1)[x] + getTime(v2)[y])/2;
    meanCurveCoord[i] = meanCoord;
    meanCurveTime[i] = meanTime;
  }
  Vector meanCurve = initTimeSeries(meanCurveCoord,meanCurveTime,"meanCurve",pathLength);
  //free(optimalPath)
  // free(meanCurveCoord)
  // free(meanCurveTime)
  return meanCurve;


}

Vector computeFrechetMeanCurve(List list,int count){
  Tree tree = createTreeFromList(list,count);
  printf("createTreeFromList OK!!!!!!!!\n");
  Vector newCenter = treeFindMeanCurve(tree);
  printf("treeFindMeanCurve OK!!!\n" );
  fflush(stdout);
  // destroyTree(tree);
  return newCenter;
}
