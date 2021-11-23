#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define TRUE 1
#define FALSE 0
extern int d;

typedef struct extra_info_node{
  int assignedCluster; // index of cluster that the vector assigned
  int iterationAssigned; // iteration that the vector assigned
  double assignedAtRadius; // radius that the vector assigned
}extraInfoNode;
typedef extraInfoNode *extraInfo;


typedef struct vec_node{
  char *vec_id; // to save the corresponding id from the given file
  double* coords; // an array to save the coordinates of the vector
  extraInfo clusterInfo; // extra info for the vector tha used at clustering (reverseAssignment with LSH and reverseAssignment with Hypercube)
}vec;
typedef vec *Vector;

double* getCoords(Vector v){
  return v->coords;
}

int assignedToCluster(Vector v){
  return (v->clusterInfo->assignedCluster==-1) ? FALSE : TRUE;
}

int getAssignedCluster(Vector v){
  return v->clusterInfo->assignedCluster;
}

void setAssignedCluster(Vector v,int cluster){
   v->clusterInfo->assignedCluster = cluster;
}

int getAssignedIteration(Vector v){
  return v->clusterInfo->iterationAssigned;
}

void setAssignedIteration(Vector v,int iter){
   v->clusterInfo->iterationAssigned = iter;
}

double getAssignedAtRadius(Vector v){
  return v->clusterInfo->assignedAtRadius;
}

void setAssignedAtRadius(Vector v,double radius){
   v->clusterInfo->assignedAtRadius = radius;
}

Vector initVector(double *vec, char id[]){
  Vector v=malloc(sizeof(struct vec_node));
  v->coords = malloc(d*sizeof(double));
  for(int i=0;i<d;i++){
    (v->coords)[i] = vec[i];
  }
  v->vec_id = malloc((strlen(id)+1)*sizeof(char));
  strcpy(v->vec_id,id);
  v->clusterInfo = NULL;
  return v;
}

void initializeClusterInfo(Vector v){
  v->clusterInfo = malloc(sizeof(extraInfoNode));
  v->clusterInfo->assignedCluster = -1;
  v->clusterInfo->iterationAssigned = -1;
  v->clusterInfo->assignedAtRadius = -1;
}


Vector copyVector(Vector vec){
  double *coords = getCoords(vec);
  Vector v=malloc(sizeof(struct vec_node));
  v->coords = malloc(d*sizeof(double));
  for(int i=0;i<d;i++){
    (v->coords)[i] = coords[i];
  }
  v->vec_id = malloc((strlen(vec->vec_id)+1)*sizeof(char));
  strcpy(v->vec_id,vec->vec_id);
  v->clusterInfo=NULL;
  return v;
}


void deleteVector(Vector v){
  free(v->coords);
  free(v->vec_id);
  if(v->clusterInfo!=NULL){
    free(v->clusterInfo);
  }
  free(v);
}


void printVector(Vector v){
  if(v==NULL)
    return;
  printf("\n[");
  for(int i=0;i<d;i++){
    printf(" %f",v->coords[i]);
  }
  printf(" ]\n");
}

void printVectorId(Vector v){
  if(v==NULL)
    return;
  printf("%s\n",v->vec_id);
}

void printVectorIdInFileNoNewline(Vector v,FILE *fptr){
  if(v==NULL)
    return;
  fprintf(fptr,"%s",v->vec_id);
}

void printVectorIdInFile(Vector v,FILE *fptr){
  if(v==NULL)
    return;
  fprintf(fptr,"%s\n",v->vec_id);
}

void printVectorInFile(Vector v,FILE *fptr){
  if(v==NULL)
    return;
  fprintf(fptr,"\n[");
  for(int i=0;i<d;i++){
    fprintf(fptr," %f",v->coords[i]);
  }
  fprintf(fptr," ]\n");
}

int compareVectors(Vector v1,Vector v2){
  for(int i=0;i<d;i++){
    if(v1->coords[i]!=v2->coords[i])
      return 0;
  }
  return 1;
}
