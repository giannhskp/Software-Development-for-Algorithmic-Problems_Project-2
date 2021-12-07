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
  if(isinf(temp)){
    printf("BKLJBVSDKLBDVKLBJVSDLJKBDVJKLSDVBJKLVDSBVJKDBSDVKLBSDVJKL\n");
    printf(" x1 = %f | y1 = %f | x2 = %f | y2 = %f \n",x1,y1,x2,y2);
    exit(0);
  }
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
  // int *optimalPath = malloc((i_dim+j_dim)*sizeof(int));
  Index **backtrackingPath=malloc(i_dim*sizeof(Index*));
  Index *optimalPath = malloc((i_dim+j_dim)*sizeof(Index));


  // printf("V1 = [ ");
  // for(int i=0;i<i_dim;i++){
  //     printf(" (%f,%f) ",coords1[i],coordsTime1[i]);
  // }
  // printf(" ]\n");
  //
  // printf("V2 = [ ");
  // for(int i=0;i<j_dim;i++){
  //     printf(" (%f,%f) ",coords2[i],coordsTime2[i]);
  // }
  // printf(" ]\n");

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
        // printf("1 ([%d,%d]) backtrackingPath = %d\n",i,j, (i-1)*(i_dim-1)+(j));
      }else if(minC==dynamicArray[i-1][j-1]){
        backtrackingPath[i][j].i=i-1;
        backtrackingPath[i][j].j=j-1;
        // printf("2 ([%d,%d]) backtrackingPath = %d\n",i,j, (i-1)*(i_dim-1)+(j-1));
      }else{
        backtrackingPath[i][j].i=i;
        backtrackingPath[i][j].j=j-1;
        // printf("3 ([%d,%d]) backtrackingPath = %d\n",i,j, (i)*(i_dim-1)+(j-1));
      }
      dynamicArray[i][j] = MAX(minC,piqj);
    }
  }
  int optimalPathLenth = 0;
  optimalPath[optimalPathLenth].i = i_dim-1;
  optimalPath[optimalPathLenth++].j = j_dim-1;
  // printf("** [%d,%d]\n",optimalPath[optimalPathLenth-1].i,optimalPath[optimalPathLenth-1].j );

  Index backtrackIndex = backtrackingPath[i_dim-1][j_dim-1];
  // printf("** [%d,%d]\n",backtrackIndex.i,backtrackIndex.j );
  optimalPath[optimalPathLenth++] = backtrackIndex;


  // printf("OPTIMAL PATH = [ ");
  // printf("-----------------------------------------\n");
  while(1){

    // printf(" backtrackIndex=%d | x = %d | y = %d    |*| [%d,%d]\n",backtrackIndex,x,y,i_dim,j_dim);
    backtrackIndex = backtrackingPath[backtrackIndex.i][backtrackIndex.j];
    if(optimalPathLenth>=(i_dim+j_dim)){
      printf("!!! optimalPath OUT OF RANGE!\n");
      break;
    }
    if(backtrackIndex.i==-1 && backtrackIndex.j==-1){
      // optimalPath[optimalPathLenth++] = -1;
      break;
    }
    optimalPath[optimalPathLenth++] = backtrackIndex;
    // printf("** [%d,%d]\n",optimalPath[optimalPathLenth-1].i,optimalPath[optimalPathLenth-1].j );
    // printf("%d ",backtrackIndex);
  }
  // printf("-----------------------------------------\n");
  // printf(" ]\n");

  // printf("OPTIMAL PATH 1 = [ ");
  // for(int i=0;i<optimalPathLenth;i++){
  //   printf("(%d,%d) ",optimalPath[i].i,optimalPath[i].j);
  // }
  // printf(" ]\n");


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
  Index *optimalPath = discreteFrechet_optimalPath(v1,v2,&pathLength);
  // printf(" PATH LENGTH = %d\n",pathLength);
  double *meanCurveCoord = malloc(pathLength*sizeof(double));
  double *meanCurveTime = malloc(pathLength*sizeof(double));
  // printf("-----------------------------------------------------\n");
  // printf("OPTIMAL PATH 2 = [ ");
  for(int i=0;i<pathLength;i++){
    Index index = optimalPath[i];
    // printf("(%d,%d) ",optimalPath[i].i,optimalPath[i].j);
    // int x = index / (i_dim-1);
    // int y = index % (i_dim-1);
    // if(x>=getDim(v1))
    //   printf("OUT OF BOUNDS V1\n");
    // if(y>=getDim(v2))
    //   printf("OUT OF BOUNDS V2\n");
    // printf("{%d} [index = %d]    V1 SIZE = %d | x = %d  |*| V2 SIZE = %d | y = %d\n",pathLength-2-i,index,getDim(v1),x,getDim(v2),y);
    double t1 = getCoords(v1)[index.i];
    double t2 = getCoords(v2)[index.j];
    double meanCoord = (t1+t2)/2;
    double t3 = getTime(v1)[index.i];
    double t4 = getTime(v2)[index.j];
    double meanTime = (t3+t4)/2;
    // double meanCoord = (getCoords(v1)[x] + getCoords(v2)[y])/2;
    // double meanTime = (getTime(v1)[x] + getTime(v2)[y])/2;
    meanCurveCoord[pathLength-1-i] = meanCoord;
    meanCurveTime[pathLength-1-i] = meanTime;
  }
  // printf(" ]\n");
  // getchar();
  // printf("-----------------------------------------------------\n");
  Vector meanCurve = initTimeSeries(meanCurveCoord,meanCurveTime,"meanCurve",pathLength);
  //free(optimalPath)
  // free(meanCurveCoord)
  // free(meanCurveTime)
  return meanCurve;


}

Vector computeFrechetMeanCurve(List list,int count){
  Tree tree = createTreeFromList(list,count);
  Vector newCenter = treeFindMeanCurve(tree);
  fflush(stdout);
  // destroyTree(tree);
  return newCenter;
  // return NULL;
}

Vector computeFrechetMeanCurveLSH(HashTable ht,int count){
  Tree tree = createTreeFromHt(ht,count);
  Vector newCenter = treeFindMeanCurve(tree);
  fflush(stdout);
  // destroyTree(tree);
  return newCenter;
  // return NULL;
}
