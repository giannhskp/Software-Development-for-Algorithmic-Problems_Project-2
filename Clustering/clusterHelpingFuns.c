#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../hashTable/hashTable.h"
#include "../LSH/lsh.h"

#define SQUARE(x) ((x)*(x))

#define CONVERGENCE 5
#define TRUE 1
#define FALSE 0

extern int d;

int existsInArray(int *array,int check,int arraySize){
  for(int i=0;i<arraySize;i++){
    if(array[i]==-1){
      return FALSE;
    }
    if(array[i]==check){
      return TRUE;
    }
  }
  return FALSE;
}


void minDistToCentroids(Vector v,Vector* vecs,Vector *clusters,int numOfClusters,double *minDistance){
  for(int i=0;i<numOfClusters;i++){
    if(clusters[i]==NULL){
      break;
    }
    double tempDist = distance_metric(clusters[i],v,d);
    if(tempDist<(*minDistance)){
        (*minDistance) = tempDist;
    }
  }
}

void minDistbetweenCentroids(Vector *centroids,int numOfClusters,double *minDistance){
  for(int i=0;i<numOfClusters;i++){
    if(centroids[i]==NULL){
      break;
    }
    for(int j=i+1;j<numOfClusters;j++){
        if(centroids[j]==NULL){
          break;
        }
        double tempDist = distance_metric(centroids[i],centroids[j],d);
        if(tempDist<(*minDistance)){
          (*minDistance) = tempDist;
        }
    }
  }
}

int centroidsConverge(Vector *new,Vector *old,int numOfClusters,int d){
  if(old==NULL) return FALSE;
  for(int i=0;i<numOfClusters;i++){
    if(distance_metric(new[i],old[i],d)>CONVERGENCE){
      return FALSE;
    }
  }
  return TRUE;
}

int findClosestCentroid(Vector v,Vector *clusters,int numOfClusters){
  int minDistIndex = -1;
  double minDist = DBL_MAX;
  for(int i=0;i<numOfClusters;i++){
    double tempDist = distance_metric(v,clusters[i],d);
    if(tempDist<minDist){
      minDistIndex = i;
      minDist = tempDist;
    }
  }
  return minDistIndex;
}
