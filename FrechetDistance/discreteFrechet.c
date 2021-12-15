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

Index *discreteFrechet_optimalPath(Vector v1,Vector v2,int *pathLength){
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

  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime1[0],coords2[0],coordsTime2[0]);
  backtrackingPath[0][0].i=-1;
  backtrackingPath[0][0].j=-1;
  for(int j=1;j<j_dim;j++){
    double p1qj = l2_metric(coords1[0],coordsTime1[0],coords2[j],coordsTime2[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
    backtrackingPath[0][j].i=0;
    backtrackingPath[0][j].j=j-1;
  }

  for(int i=1;i<i_dim;i++){
    double piq1 = l2_metric(coords1[i],coordsTime1[i],coords2[0],coordsTime2[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
    backtrackingPath[i][0].i=i-1;
    backtrackingPath[i][0].j=0;
  }

  for(int i=1;i<i_dim;i++){
    for(int j=1;j<j_dim;j++){
      double piqj = l2_metric(coords1[i],coordsTime1[i],coords2[j],coordsTime2[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
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
  int optimalPathLenth = 0;
  optimalPath[optimalPathLenth].i = i_dim-1;
  optimalPath[optimalPathLenth++].j = j_dim-1;

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
  return optimalPath;
}


Vector meanCurveBetween2Curves(Vector v1,Vector v2){
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
