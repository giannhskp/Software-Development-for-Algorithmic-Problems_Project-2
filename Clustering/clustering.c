#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../hashTable/hashTable.h"
#include "../LSH/lsh.h"
#include "../Hypercube/hypercube.h"
#include "./clusterHelpingFuns.h"
#include "./kmeansPlusPlus.h"
#include "../FrechetDistance/discreteFrechet.h"

#define SQUARE(x) ((x)*(x))

#define TRUE 1
#define FALSE 0
#define MAX_RECENTER_ITERATIONS 10
#define W_DIVIDER_LSH 60
#define W_DIVIDER_CUBE 20
#define METHOD_VECTOR 2
#define METHOD_FRECHET 3

extern int numOfVecs;
// extern int d;
extern int hashTableSize;
extern int silhouette;
extern int w;
char *distanceMetric;
// Vector timeVector;


int wValueCalculation(List list,int numberOfVectorsInFile,int dim){
  long double sumDist = 0.0;
  int count=0;
  double persentageToCheck;
  if(numberOfVectorsInFile<=1000){
    persentageToCheck = 0.1;
  }else if(numberOfVectorsInFile<=10000){
    persentageToCheck = 0.001;
  }else if (numberOfVectorsInFile<=100000){
    persentageToCheck = 0.0001;
  }else{
    persentageToCheck = 0.000001;
  }
  persentageToCheck = 0.0001; // TODO: REMOVE
  int stopBound = persentageToCheck*numberOfVectorsInFile*numberOfVectorsInFile;
  while(list!=NULL){
    List nested = list;
    while(nested!=NULL){
      if(count>stopBound){
        return floor(sumDist/count);
      }
      sumDist += distance_metric(getVector(list),getVector(nested));
      count++;
      nested = getNext(nested);
    }
    list=getNext(list);
  }
  return floor(sumDist/count);
}



void lloyds(Vector* clusters,Vector *oldClusters,Vector* vectors,List* clustersList,int numberOfVectors,int numOfClusters,int *vectorCount,int *firstTime,int dim) {
  // lloyds Algorithm
  if(*firstTime) // skip it for the first time (original centroids from the kmeansPlusPlus)
    for(int i=0;i<numOfClusters;i++){
      Vector newCenter;
      if(clustersList[i]!=NULL){ // check if each cluster has been formed (has vectors)
        // ok then find the new centroid for this cluster
        if(strcmp(distanceMetric,"discreteFrechet")==0){
          newCenter = computeFrechetMeanCurve(clustersList[i],vectorCount[i]);
        }else{
          newCenter=listMeanOfCluster(clustersList[i],dim);
        }
      }
      else{
        // this cluster hasn't been formed, let as centroid the previous one
        newCenter=copyVector(oldClusters[i]);
      }
      // finally delete each cluster in order to form a new one based to the new centroid
      listDelete(clustersList[i],0);
      clustersList[i] = NULL;
      // save the new centroid
      clusters[i]=newCenter;
      vectorCount[i] = 0;
    }

  for(int i=0;i<numberOfVectors;i++){ // for every vector
    // find the closest centroid with the euclidean distance
    int closestCentroid = findClosestCentroid(vectors[i],clusters,numOfClusters);
    // printf("Vector[%d] | closestCentroid= %d\n",i,closestCentroid);
    // and assign this vector to the corresponding cluster
    vectorCount[closestCentroid] += 1;
    clustersList[closestCentroid] = listInsert(clustersList[closestCentroid],vectors[i],dim);
  }
  *firstTime=1;
}

double *silhouetteLloyds(List *clustersList,Vector *clusters,int numOfClusters,int *vectorCount,double *stotal,int dim){
  // used to find the silhouettes of each cluster in Lloyds Algorithm
  double *silhouettes = calloc(sizeof(double),numOfClusters);
  for(int i=0;i<numOfClusters;i++){
    silhouettes[i] = silhouetteofClusterLloyds(clustersList,clusters,i,numOfClusters,vectorCount[i],dim,stotal);
  }
  return silhouettes;
}

void clusteringLloyds(List vecList,int numOfClusters,FILE* fptr,int dim){
  Vector *vectors;
  Vector *clusters;
  Vector *oldClusters = NULL;
  double *props;
  int *vectorCount;

  vectors = transformListToArray(vecList,numOfVecs);
  clusters = malloc(numOfClusters*sizeof(Vector)); // used to store the new centroids
  oldClusters = malloc(numOfClusters*sizeof(Vector)); // used to store the old centroids

  for(int i=0;i<numOfClusters;i++){
    clusters[i] = NULL;
    oldClusters[i] = NULL;
  }
  props = calloc(numOfVecs,sizeof(double));
  vectorCount = calloc(numOfVecs,sizeof(int));

  clock_t cluster_start = clock();

  // find the original centroids with the kmeans++ Algorithm
  kmeansplusplus(vectors,numOfClusters,clusters,props);

  int firstIter = TRUE;

  List *clustersList=malloc(numOfClusters*sizeof(List *)); // array of lists, each list represents a cluster that the vectors will be stored
  for(int i=0;i<numOfClusters;i++){
    clustersList[i]=initializeList();
  }
  int count=0;
  int firstTime=0;
  // lloyds Algorithm runs until convergence between the old cluster centroids and the new ones is achieved
  while((count<2) || !centroidsConverge(clusters,oldClusters,numOfClusters)){ // check for convergence after the second one iteration
    count++;
    if(!firstIter){
      Vector *temp = oldClusters;
      oldClusters=clusters;
      clusters = temp;
      for(int i=0;i<numOfClusters;i++){
        if(clusters[i]!=NULL){
          deleteVector(clusters[i]);
          clusters[i] = NULL;
        }
        // vectorCount[i] = 0;
      }
    }
    // lloyds Algorithm
    lloyds(clusters,oldClusters,vectors,clustersList,numOfVecs,numOfClusters,vectorCount,&firstTime,dim);

    firstIter=FALSE;
  }

  clock_t cluster_end = clock();
  double cluster_time = (double)(cluster_end - cluster_start) / CLOCKS_PER_SEC;

  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"CLUSTER-%d {size: %d",i+1,vectorCount[i]);
    printVectorInFile(clusters[i],fptr);
    fprintf(fptr,"}\n");
  }

  fprintf(fptr, "clustering_time: %f seconds\n",cluster_time);
  fflush(fptr);

  printf("- COMPUTING SILHOUETTES FOR CLUSTERS\n");
  double stotal = 0.0;
  double * silhouettes = silhouetteLloyds(clustersList,clusters,numOfClusters,vectorCount,&stotal,dim);
  printf("- FINISHED COMPUTING SILHOUETTES\n");
  fprintf(fptr, "Silhouette: [ ");
  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"s%d = %f ,",i+1,silhouettes[i]);
  }
  fprintf(fptr,"stotal = %f ]\n\n",stotal/numOfVecs);

  if(silhouette){
    for(int i=0;i<numOfClusters;i++){
      fprintf(fptr,"CLUSTER-%d {",i+1);
      printVectorInFile(clusters[i],fptr);
      listPrintClusteringInFile(clustersList[i],fptr);
      fprintf(fptr,"}\n\n");
    }
  }

  // free the allocated space
  for(int i=0;i<numOfClusters;i++){
    listDelete(clustersList[i],0);
    deleteVector(oldClusters[i]);
    deleteVector(clusters[i]);
  }

  free(silhouettes);
  free(clusters);
  free(props);
  free(vectors);
  free(clustersList);
  free(oldClusters);
  free(vectorCount);
}

void reverseAssignmentLSH(LSH lsh,Vector *vectors,Vector *clusters,Vector *oldClusters,HashTable *clustersHt,int numOfClusters,int iteration,int *firstTime,int dim,int method,Grids grids,double delta){
  if(*firstTime) // skip it for the first time (original centroids from the kmeansPlusPlus)
    for(int i=0;i<numOfClusters;i++){
      Vector newCenter = NULL;
      if(getNumberOfVectors(clustersHt[i])<=0){ // this cluster hasn't been formed, let as centroid the previous one
        newCenter=copyVector(oldClusters[i]);
      }else if(method == METHOD_VECTOR){
        newCenter = htMeanOfCluster(clustersHt[i],dim); // find the new centroid for every cluster
      }else if(method == METHOD_FRECHET){
        newCenter = computeFrechetMeanCurveLSH(clustersHt[i],getNumberOfVectors(clustersHt[i]));
      }

      // finally delete each cluster in order to form a new one based to the new centroid
      htDelete(clustersHt[i],0);
      clustersHt[i] = htInitialize(numOfVecs/(4*numOfClusters));
      // save the new centroid
      clusters[i]=newCenter;
    }
  double radius=DBL_MAX;
  // find the min distance between the centroids in order to initialize the radius for the range search
  minDistbetweenCentroids(clusters,numOfClusters,&radius);
  radius/=2;
  int assignCounter = 0;
  int previousAssigns = -1;
  int loopCounter = 0;
  while((double)assignCounter<(0.8*numOfVecs) && loopCounter<MAX_RECENTER_ITERATIONS){ // stop when the 80% of vectors have been assigned in to the clusters
    if(assignCounter==previousAssigns && assignCounter!=0 && loopCounter>5){ // or stop when the assigns number stays some with the previous one for more than 5 iterations
      break;
    }
    previousAssigns = assignCounter;
    List confList=initializeList(); // list that used to store the vectors that in range search assigned at more than one cluster
    // assign each vector to the corresponding cluster with the help of range search
    for(int i=0;i<numOfClusters;i++){
      if(method == METHOD_VECTOR){
        radiusNeigborsClustering(lsh,clusters[i],radius,clustersHt[i],i,&confList,&assignCounter,iteration);
      }else if(method == METHOD_FRECHET){
        radiusNeigborsClusteringTimeSeries(lsh,clusters[i],radius,clustersHt[i],i,&confList,&assignCounter,iteration,grids,delta);
      }
    }
    // manage the vectors that presenting conflict
    listSolveRangeConflicts(confList,clustersHt,clusters,numOfClusters,dim,iteration);
    listDelete(confList,0);
    radius*=2; // doubled the radius for the next range search
    loopCounter++;
  }
  int remainderCounter = 0;
  // finally one big percentage of vectors has been assigned into clusters
  // the remaining vectors will be assigned based on the nearest centroid at the corresponding cluster
  if(assignCounter<numOfVecs){
    for(int i=0;i<numOfVecs;i++){
      if(assignedToCluster(vectors[i]) && (getAssignedIteration(vectors[i])==iteration)){
        continue;
      }
      int closestCentroid = findClosestCentroid(vectors[i],clusters,numOfClusters);
      htRangeInsert(clustersHt[closestCentroid],vectors[i],-1,dim);
      setAssignedCluster(vectors[i],closestCentroid);
      setAssignedIteration(vectors[i],iteration);
      setAssignedAtRadius(vectors[i],radius);
      remainderCounter++;
    }
  }
  *firstTime=1;
}

double *silhouetteLSH_Hypercube(HashTable *clustersHt,Vector *clusters,int numOfClusters,double *stotal,int dim){
    // used to find the silhouettes of each cluster in reverseAssignmentHypercube
  double *silhouettes = calloc(sizeof(double),numOfClusters);
  for(int i=0;i<numOfClusters;i++){
    silhouettes[i] = silhouetteofClusterLSH(clustersHt,clusters,i,numOfClusters,dim,stotal);
  }
  return silhouettes;
}

void clusteringLSH(List vecList,int numOfClusters,int l,FILE* fptr,int dim,int method,int delta){
  Vector *vectors;
  Vector *clusters;
  Vector *oldClusters = NULL;
  double *props;
  HashTable *clustersHt=malloc(numOfClusters*sizeof(HashTable *)); // array of hash tables each hash table represents a cluster that the vectors will be stored
  for(int i=0;i<numOfClusters;i++){
    clustersHt[i]= htInitialize(numOfVecs/(4*numOfClusters));
  }
  vectors = transformListToArray(vecList,numOfVecs);
  clusters = malloc(numOfClusters*sizeof(Vector)); // used to store the new centroids
  oldClusters = malloc(numOfClusters*sizeof(Vector)); // used to store the old centroids
  for(int i=0;i<numOfClusters;i++){
    clusters[i]=NULL;
    oldClusters[i]=NULL;
  }
  props = calloc(numOfVecs,sizeof(double));

  // allocate and initialize the LSH with the vectors tha will be inserted into clusters
  if(numOfVecs>1000){
    hashTableSize = (int)(numOfVecs*0.005);
  }else{
    hashTableSize=numOfVecs/32;
  }

  clock_t begin = clock();
  // w = wValueCalculation(vecList,numOfVecs,dim);
  // w /= W_DIVIDER_LSH;
  w = 6;
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  // printf("Found value of w in %f seconds, w = %d\n",time_spent,w );



  begin = clock();
  int vector_v_dim = -1;
  if(method == METHOD_VECTOR){
    vector_v_dim = dim;
  }else if(method == METHOD_FRECHET){
    vector_v_dim = 2*dim*numOfVecs; // TODO: ADD AN UPPER BOUND
  }
  LSH lsh = initializeLSH(l,vector_v_dim);
  Grids grids;
  if(method == METHOD_VECTOR){
    grids = NULL;
  }else if(method == METHOD_FRECHET){
    grids = initializeGrids(delta,l);
  }
  for(int i=0;i<numOfVecs;i++){
    initializeClusterInfo(vectors[i]);
    if(method == METHOD_VECTOR){
      insertToLSH(lsh,vectors[i]);
    }else if(method == METHOD_FRECHET){
      insertTimeSeriesToLSH(lsh,grids,delta,vectors[i]);
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Created LSH in : %f seconds\n",time_spent);


  clock_t cluster_start = clock();
  // find the original centroids with the kmeans++ Algorithm
  kmeansplusplus(vectors,numOfClusters,clusters,props);

  printf("Initialized clusters with Kmeans++\n");

  int firstIterLSH = TRUE;
  int countLSH=0;
  int firstTime=0;
  // reverseAssignmentLSH Algorithm runs until convergence between the old cluster centroids and the new ones is achieved
  while((countLSH<2) || !centroidsConverge(clusters,oldClusters,numOfClusters)){ // check for convergence after the second one iteration
    if(countLSH==MAX_RECENTER_ITERATIONS)
      break;
    countLSH++;
    if(!firstIterLSH){
      Vector *temp = oldClusters;
      oldClusters=clusters;
      clusters = temp;
      for(int i=0;i<numOfClusters;i++){
        if(clusters[i]!=NULL){
          deleteVector(clusters[i]);
          clusters[i] = NULL;
        }
      }
    }

    // reverseAssignmentLSH Algorithm
    reverseAssignmentLSH(lsh,vectors,clusters,oldClusters,clustersHt,numOfClusters,countLSH,&firstTime,dim,method,grids,delta);

    firstIterLSH=FALSE;

  }

  clock_t cluster_end = clock();
  double cluster_time = (double)(cluster_end - cluster_start) / CLOCKS_PER_SEC;

  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"CLUSTER-%d {size: %d",i+1,getNumberOfVectors(clustersHt[i]));
    printVectorInFile(clusters[i],fptr);
    fprintf(fptr,"}\n");
  }

  fprintf(fptr, "clustering_time: %f seconds\n",cluster_time);
  fflush(fptr);
  fflush(stdout);

  printf("- COMPUTING SILHOUETTES FOR CLUSTERS\n");
  double stotal = 0.0;
  double * silhouettes = silhouetteLSH_Hypercube(clustersHt,clusters,numOfClusters,&stotal,dim);
  printf("- FINISHED COMPUTING SILHOUETTES\n");
  fprintf(fptr, "Silhouette: [ ");
  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"s%d = %f ,",i+1,silhouettes[i]);
  }
  fprintf(fptr,"stotal = %f ]\n\n",stotal/numOfVecs);

  if(silhouette){
    for(int i=0;i<numOfClusters;i++){
      fprintf(fptr,"CLUSTER-%d {",i+1);
      printVectorInFile(clusters[i],fptr);
      htPrintClustering(clustersHt[i],fptr);
      fprintf(fptr,"}\n\n");
    }
  }

  // free the allocated space
  for(int i=0;i<numOfClusters;i++){
    if(oldClusters[i]!=NULL)
      deleteVector(oldClusters[i]);
    if(clusters[i]!=NULL)
      deleteVector(clusters[i]);
    htDelete(clustersHt[i],0);
  }

  free(silhouettes);
  free(props);
  free(vectors);
  free(oldClusters);
  free(clustersHt);
  free(clusters);
  destroyLSH(lsh);
}

void reverseAssignmentHypercube(HyperCube cube,Vector *vectors,Vector *clusters,Vector *oldClusters,HashTable *clustersHt,int numOfClusters,int iteration,int m,int probes,int *firstTime,int dim){
  if(*firstTime) // skip it for the first time (original centroids from the kmeansPlusPlus)
    for(int i=0;i<numOfClusters;i++){

      Vector newCenter = htMeanOfCluster(clustersHt[i],dim); // find the new centroid for every cluster
      if(newCenter==NULL){ // this cluster hasn't been formed, let as centroid the previous one
        newCenter=copyVector(oldClusters[i]);
      }
        // finally delete each cluster in order to form a new one based to the new centroid
      htDelete(clustersHt[i],0);
      clustersHt[i] = htInitialize(numOfVecs/(4*numOfClusters));
      // save the new centroid
      clusters[i]=newCenter;
    }
  double radius=DBL_MAX;
  // find the min distance between the centroids in order to initialize the radius for the range search
  minDistbetweenCentroids(clusters,numOfClusters,&radius);
  radius/=2;
  int assignCounter = 0;
  int previousAssigns = -1;
  int loopCounter = 0;
  while((double)assignCounter<(0.8*numOfVecs) && loopCounter<MAX_RECENTER_ITERATIONS){ // stop when the 80% of vectors have been assigned in to the clusters
    if(assignCounter==previousAssigns && assignCounter!=0){ // or stop when the assigns number stays some with the previous one for more than 5 iterations
      break;
    }

    previousAssigns = assignCounter;
    List confList=initializeList(); // list that used to store the vectors that in range search assigned at more than one cluster
    // assign each vector to the corresponding cluster with the help of range search
    for(int i=0;i<numOfClusters;i++){
      radiusNeigborHypercubeClustering(cube,clusters[i],clustersHt[i],radius,probes,m,i,&confList,&assignCounter,iteration);
    }
    // manage the vectors that presenting conflict
    listSolveRangeConflicts(confList,clustersHt,clusters,numOfClusters,dim,iteration);
    listDelete(confList,0);
    radius*=2; // doubled the radius for the next range search
    loopCounter++;
  }
  printf("ASSIGN COUNTER = %d\n",assignCounter);
  int remainderCounter = 0;
  // finally one big percentage of vectors has been assigned into clusters
  // the remaining vectors will be assigned based on the nearest centroid at the corresponding cluster
  if(assignCounter<numOfVecs){
    for(int i=0;i<numOfVecs;i++){
      if(assignedToCluster(vectors[i]) && (getAssignedIteration(vectors[i])==iteration)){
        continue;
      }
      int closestCentroid = findClosestCentroid(vectors[i],clusters,numOfClusters);
      htRangeInsert(clustersHt[closestCentroid],vectors[i],-1,dim);
      setAssignedCluster(vectors[i],closestCentroid);
      setAssignedIteration(vectors[i],iteration);
      setAssignedAtRadius(vectors[i],radius);
      remainderCounter++;
    }
  }
  printf("remainderCounter COUNTER = %d\n",remainderCounter);
  *firstTime=1;
}

void clusteringHypercube(List vecList,int numOfClusters,int m,int probes,FILE* fptr,int dim){
  Vector *vectors;
  Vector *clusters;
  Vector *oldClusters = NULL;
  double *props;
  HashTable *clustersHt=malloc(numOfClusters*sizeof(HashTable *)); // array of hash tables each hash table represents a cluster that the vectors will be stored
  for(int i=0;i<numOfClusters;i++){
    clustersHt[i]= htInitialize(numOfVecs/(4*numOfClusters));
  }
  vectors = transformListToArray(vecList,numOfVecs);
  clusters = malloc(numOfClusters*sizeof(Vector)); // used to store the new centroids
  oldClusters = malloc(numOfClusters*sizeof(Vector)); // used to store the old centroids
  for(int i=0;i<numOfClusters;i++){
    clusters[i]=NULL;
    oldClusters[i]=NULL;
  }
  props = calloc(numOfVecs,sizeof(double));

  clock_t begin = clock();
  // w = wValueCalculation(vecList,numOfVecs,dim);
  // w /= W_DIVIDER_CUBE;
  w = 6;
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  // printf("Found value of w in %f seconds, w = %d\n",time_spent,w );

  // allocate and initialize the Hypercube with the vectors tha will be inserted into clusters
  hashTableSize=numOfVecs/16;
  begin = clock();
  HyperCube cube = initializeHyperCube(dim);
  for(int i=0;i<numOfVecs;i++){
    initializeClusterInfo(vectors[i]);
    insertToHyperCube(cube,vectors[i]);
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Created Hypercube in : %f seconds\n",time_spent);

  clock_t cluster_start = clock();
  // find the original centroids with the kmeans++ Algorithm
  kmeansplusplus(vectors,numOfClusters,clusters,props);

  int firstIterLSH = TRUE;
  int countLSH=0;
  int firstTime=0;
  // reverseAssignmentHypercube Algorithm runs until convergence between the old cluster centroids and the new ones is achieved
  while((countLSH<2) || !centroidsConverge(clusters,oldClusters,numOfClusters)){
    if(countLSH==MAX_RECENTER_ITERATIONS)
      break;
    countLSH++;
    if(!firstIterLSH){
      Vector *temp = oldClusters;
      oldClusters=clusters;
      clusters = temp;
      for(int i=0;i<numOfClusters;i++){
        if(clusters[i]!=NULL){
          deleteVector(clusters[i]);
          clusters[i] = NULL;
        }
      }
    }

    // reverseAssignmentHypercube Algorithm
    reverseAssignmentHypercube(cube,vectors,clusters,oldClusters,clustersHt,numOfClusters,countLSH,m,probes,&firstTime,dim);

    firstIterLSH=FALSE;

  }

  clock_t cluster_end = clock();
  double cluster_time = (double)(cluster_end - cluster_start) / CLOCKS_PER_SEC;

  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"CLUSTER-%d {size: %d",i+1,getNumberOfVectors(clustersHt[i]));
    printVectorInFile(clusters[i],fptr);
    fprintf(fptr,"}\n");
  }

  fprintf(fptr, "clustering_time: %f seconds\n",cluster_time);
  fflush(fptr);
  fflush(stdout);

  printf("- COMPUTING SILHOUETTES FOR CLUSTERS\n");
  double stotal = 0.0;
  double * silhouettes = silhouetteLSH_Hypercube(clustersHt,clusters,numOfClusters,&stotal,dim);
  printf("- FINISHED COMPUTING SILHOUETTES\n");
  fprintf(fptr, "Silhouette: [ ");
  for(int i=0;i<numOfClusters;i++){
    fprintf(fptr,"s%d = %f ,",i+1,silhouettes[i]);
  }
  fprintf(fptr,"stotal = %f ]\n\n",stotal/numOfVecs);

  if(silhouette){
    for(int i=0;i<numOfClusters;i++){
      fprintf(fptr,"CLUSTER-%d {",i+1);
      printVectorInFile(clusters[i],fptr);
      htPrintClustering(clustersHt[i],fptr);
      fprintf(fptr,"}\n\n");
    }
  }
  // free the allocated space
  for(int i=0;i<numOfClusters;i++){
    if(oldClusters[i]!=NULL)
      deleteVector(oldClusters[i]);
    if(clusters[i]!=NULL)
      deleteVector(clusters[i]);
    htDelete(clustersHt[i],0);
  }

  free(silhouettes);
  free(props);
  free(vectors);
  free(oldClusters);
  free(clustersHt);
  free(clusters);
  deleteHyperCube(cube);
}

void clustering(List vecList,FILE* fptr,char* assignment,char *update,int numOfClusters,int l,int mHyper,int probes,int dim,int delta){
  if(strcmp(assignment,"Classic")==0){
    if(strcmp(update,"Mean Vector")==0){
      fprintf(fptr,"Algorithm: Lloyds for Vectors\n");
      distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));
      strcpy(distanceMetric,"l2");
    }else if(strcmp(update,"Mean Frechet")==0){
      fprintf(fptr,"Algorithm: Lloyds for Curves\n");
      distanceMetric=malloc(sizeof(char)*(strlen("discreteFrechet")+1));
      strcpy(distanceMetric,"discreteFrechet");
    }else{
      printf("Wrong update method!\n");
      exit(-1);
    }
    fprintf(fptr,"Algorithm: Lloyds\n");
    clusteringLloyds(vecList,numOfClusters,fptr,dim);
  }
  else if(strcmp(assignment,"LSH")==0){
    int method;
    if(strcmp(update,"Mean Vector")==0){
      fprintf(fptr,"Algorithm: Range Search LSH for Vectors\n");
      distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));
      strcpy(distanceMetric,"l2");
      method = METHOD_VECTOR;
    }else if(strcmp(update,"Mean Frechet")==0){
      fprintf(fptr,"Algorithm: Range Search LSH for Curves\n");
      distanceMetric=malloc(sizeof(char)*(strlen("discreteFrechet")+1));
      strcpy(distanceMetric,"discreteFrechet");
      method = METHOD_FRECHET;
    }else{
      printf("Wrong update method!\n");
      exit(-1);
    }

    clusteringLSH(vecList,numOfClusters,l,fptr,dim,method,delta);
  }
  else if(strcmp(assignment,"Hypercube")==0){
    fprintf(fptr,"Algorithm: Range Search Hypercube for vectors\n");
    distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));
    strcpy(distanceMetric,"l2");
    clusteringHypercube(vecList,numOfClusters,mHyper,probes,fptr,dim);
  }
  else{
    printf("%s\n",assignment );
    printf("INVALID assignment NAME!\n");
    exit(EXIT_FAILURE);
  }

}
