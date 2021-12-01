#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "../Vector/vector.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MAX3(a, b, c) MAX((MAX((a),(b))),(c))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MIN3(a, b, c) MIN((MIN((a),(b))),(c))

static double l2_metric(double x1,double y1,double x2,double y2){
  // calculate the Euclidean distance (or L2) between the given vectors and return it
  return sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1)));
}

double discreteFrechet(Vector v1,Vector v2,Vector time,int d){
  double **dynamicArray=malloc(d*sizeof(double*));
  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  double *coordsTime = getCoords(time);

  for(int i=0;i<d;i++){
    dynamicArray[i]=malloc(d*sizeof(double));
  }

  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime[0],coords2[0],coordsTime[0]);
  for(int i=1;i<d;i++){
    double p1qj = l2_metric(coords1[0],coordsTime[0],coords2[i],coordsTime[i]);
    dynamicArray[0][i] = MAX( dynamicArray[0][i-1], p1qj);
  }

  for(int j=1;j<d;j++){
    double piq1 = l2_metric(coords1[j],coordsTime[j],coords2[0],coordsTime[0]);
    dynamicArray[j][0] = MAX(dynamicArray[j-1][0],piq1);
  }

  for(int i=1;i<d;i++){
    for(int j=1;j<d;j++){
      double piqj = l2_metric(coords1[i],coordsTime[i],coords2[j],coordsTime[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }

  double finalDistance = dynamicArray[d-1][d-1];

  for(int i=0;i<d;i++){
    free(dynamicArray[i]);
  }
  free(dynamicArray);
  return finalDistance;
}

int *discreteFrechet_optimalPath(Vector v1,Vector v2,Vector time,int d,int *pathLength){
  double **dynamicArray=malloc(d*sizeof(double*));
  int **backtrackingPath=malloc(d*sizeof(int*));
  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  double *coordsTime = getCoords(time);
  int *optimalPath = malloc(2*d*sizeof(int*));

  for(int i=0;i<d;i++){
    dynamicArray[i]=malloc(d*sizeof(double));
    backtrackingPath[i]=malloc(d*sizeof(int));
  }

  dynamicArray[0][0] = l2_metric(coords1[0],coordsTime[0],coords2[0],coordsTime[0]);
  backtrackingPath[0][0] = -1;
  for(int j=1;j<d;j++){
    double p1qj = l2_metric(coords1[0],coordsTime[0],coords2[j],coordsTime[j]);
    dynamicArray[0][j] = MAX( dynamicArray[0][j-1], p1qj);
    backtrackingPath[0][j] =  j-1;
  }

  for(int i=1;i<d;i++){
    double piq1 = l2_metric(coords1[i],coordsTime[i],coords2[0],coordsTime[0]);
    dynamicArray[i][0] = MAX(dynamicArray[i-1][0],piq1);
    backtrackingPath[i][0] =  (i-1)*d;
  }

  for(int i=1;i<d;i++){
    for(int j=1;j<d;j++){
      double piqj = l2_metric(coords1[i],coordsTime[i],coords2[j],coordsTime[j]);
      double minC = MIN3(dynamicArray[i-1][j], dynamicArray[i-1][j-1], dynamicArray[i][j-1]);
      if(minC==dynamicArray[i-1][j]){
        backtrackingPath[i][j] = (i-1)*d+(j);
      }else if(minC==dynamicArray[i-1][j-1]){
        backtrackingPath[i][j] = (i-1)*d+(j-1);
      }else{
        backtrackingPath[i][j] = (i)*d+(j-1);
      }
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }
  int optimalPathLenth = 0;
  int backtrackIndex = backtrackingPath[d-1][d-1];
  optimalPath[optimalPathLenth++] = backtrackIndex;

  while(1){
    if(backtrackIndex==-1){
      optimalPath[optimalPathLenth++] = -1;
      break;
    }
    int y = backtrackIndex / d;
    int x = backtrackIndex % d;
    backtrackIndex = backtrackingPath[x][y];
    optimalPath[optimalPathLenth++] = backtrackIndex;
  }


  (*pathLength) = optimalPathLenth;

  for(int i=0;i<d;i++){
    free(dynamicArray[i]);
    free(backtrackingPath[i]);
  }
  free(dynamicArray);
  free(backtrackingPath);
  return optimalPath;
}
