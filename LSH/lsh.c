#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTable.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "./helperFunctions.h"
// #include "../Fred-master/include/frechet.hpp"
// #include "../Fred-master/include/point.hpp"
// #include "../Fred-master/include/types.hpp"


#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

extern int w;
extern int d;
extern int k_LSH;
extern int hashTableSize;

#define PADDING_M 300


typedef struct hfunc{
  double *v; // v ∼ N (0, 1)^d
  double t; // t variable uniformly ∈ R [0, w)
}h_function;

typedef struct gfunc{
  h_function *h_functions;
  int *r;
  unsigned int m;
}g_function;

typedef struct lsh_n{
  int l;
  g_function *g_fun; // g function is a random combination of hi's, every g function has k_LSH h functions
  HashTable *hts; // LSH has l hash tables (g functions used like hash functions at the corresponding hash tables )
}lshNode;
typedef lshNode *LSH;

typedef struct grid_n{
  double *t_array;
}gridNode;
typedef gridNode *Grids;

Grids initializeGrids(double delta,int l){
  Grids grids = malloc(sizeof(gridNode));
  grids->t_array = malloc(l*sizeof(double));
  for(int i=0;i<l;i++){
    double temp = uniform_distribution(0,delta);
    grids->t_array[i] = temp;
  }
  return grids;
}

void deleteGrids(Grids grids){
  free(grids->t_array);
  free(grids);
}

double getTofGrid(Grids grids,int index){
  return grids->t_array[index];
}


int getL(LSH lsh){
  return lsh->l;
}
HashTable *getHts(LSH lsh){
  return lsh->hts;
}
g_function *getGfuns(LSH lsh){
  return lsh->g_fun;
}

Vector timeSeriesSnapping(Vector v,Vector time,int d,double gridDelta,double t){
  double *coordsVector = getCoords(v);
  double *coordsTime = getCoords(time);
  double snappedVector[d];
  double snappedTime[d];
  double snappedFinal[2*d];

  int index=0;
  int indexFinal=0;
  for(int i=0;i<d;i++){
    double temp=coordsVector[i] + t;  // displacement
    double tempTime=coordsTime[i] + t;   // displacement
    int found=0;
    double keepX;
    double keepY;

    // x
    tempTime = tempTime+(0.5);
    tempTime = tempTime/gridDelta;
    tempTime = floor(tempTime);
    tempTime = tempTime * gridDelta;
    keepX = tempTime;
    // y
    temp = temp+(0.5);
    temp = temp/gridDelta;
    temp = floor(temp);
    temp = temp * gridDelta;
    keepY = temp;


    for(int j=0;j<index;j++){
      if(snappedTime[j]==keepX && snappedVector[j]==keepY){
        found=1;
        break;
      }
    }

    if(!found){
      snappedFinal[indexFinal++]=keepX;
      snappedFinal[indexFinal++]=keepY;
      snappedTime[index]=keepX;
      snappedVector[index++]=keepY;
    }
  }
  // TODO: padding
  ////////////////////////////
  // for(int i=index;i<d;i++){
  //   snappedVector[i]=PADDING_M;
  //   snappedTime[i]=PADDING_M;
  // }
  for(int i=index;i<2*d;i++){
    snappedFinal[i]=PADDING_M;
    // snappedTime[i]=PADDING_M;
  }
  ////////////////////////////
  // printf("%d\n",d);
  d*=2;
  // printf("%d\n",d);
  Vector vecTmp=initVector(snappedFinal,getID(v));
  d/=2;
  // getchar();
  // printVector(vecTmp);
  // getchar();
  // printVector(snappedTime);
  // printVector(snappedVector);
  return vecTmp;
}

Vector continuousTimeSeriesSnapping(Vector v,int d,double gridDelta){
  double *coordsVector = getCoords(v);
  double snappedVector[d];

  int index=0;
  for(int i=0;i<d;i++){
    double temp=coordsVector[i];  // displacement
    if(temp==-1){
      continue;
    }
    double keepY;
    // y
    // temp = temp + gridDelta/2;
    // temp = temp - fmod(temp,gridDelta);
    // keepY=temp;

    temp = temp+(0.5);
    temp = temp/gridDelta;
    temp = floor(temp);
    temp = temp * gridDelta;
    keepY = temp;


    snappedVector[index++]=keepY;
  }
  // TODO: padding
  ////////////////////////////
  for(int i=index;i<d;i++){
    snappedVector[i]=PADDING_M;
  }
  ////////////////////////////
  Vector vecTmp=initVector(snappedVector,getID(v));
  return vecTmp;
}

Vector minima_maxima(Vector v,int d){
  double *coordsVector = getCoords(v);
  double keyVector[d];
  int newIndex = 0;
  int prev = 0;
  for(int i=1;i<(d-1);i++){
    if(coordsVector[i+1]==PADDING_M){
      break;
    }
    double minimum = MIN(coordsVector[prev],coordsVector[i+1]);
    double maximum = MAX(coordsVector[prev],coordsVector[i+1]);
    if(coordsVector[i]>=minimum && coordsVector[i]<=maximum){
      continue;
    }else{
      keyVector[newIndex++] = coordsVector[i];
      prev=i;
    }
  }
  for(int i=newIndex;i<d;i++){
    keyVector[i] = PADDING_M;
  }
  Vector tempVec = initVector(keyVector,"keyVector");
  return tempVec;
}

Vector filtering(Vector v,double epsilon){
  double *filteredCoords = malloc(d*(sizeof(double)));
  double *originalCoords = getCoords(v);
  filteredCoords[0]=originalCoords[0];
  double previous = originalCoords[0];
  for(int i=1;i<(d-1);i++){
    if((fabs(previous-originalCoords[i])<epsilon) && (fabs(originalCoords[i]-originalCoords[i+1])<epsilon)){
      filteredCoords[i]=-1;
    }else{
      filteredCoords[i]=originalCoords[i];
      previous = originalCoords[i];
    }
  }
  filteredCoords[d-1]=originalCoords[d-1];
  Vector tempVec = initVector(filteredCoords,getID(v));
  free(filteredCoords);
  return tempVec;
}


/* H FUNCTIONS*/

void generateH_LSH(h_function *hfun){
  // generate v vector coordinates of h function, v ∼ N (0, 1)^d
  hfun->v=malloc(d*sizeof(double));
  for(int i=0;i<d;i++){
    hfun->v[i] = normalRandom();
  }
  // pick t variable uniformly ∈ R [0, w)
  double temp = uniform_distribution(0,w);
  hfun->t=temp;
}

void destroyH_LSH(h_function h){
  free(h.v);
}

int computeH_LSH(h_function hfun,Vector vector){
  // compute the dot product of the given vector with the v vector of h function (p · v)
  double pv = dot_product(hfun.v,getCoords(vector),d);
  // finally calculate the value of h function
  double temp = (double) (pv+hfun.t)/(double)w;
  return floor(temp);
}

/* G FUNCTIONS*/

void generateG(g_function *gfun){
  // allocate and generate the h functions tha will be used at the calculation of G function, k_LSH (number of h functions) has been given from the command line
  // g function is a random combination of hi's, every g function has k_LSH h functions
  gfun->h_functions = malloc(k_LSH*sizeof(h_function));

  for(int i=0;i<k_LSH;i++){
     generateH_LSH(&gfun->h_functions[i]);
  }
  // allocate and generate as many variables r (r[i] are a random int ≤ 32 bits) as the functions h
  gfun->r = malloc(k_LSH*sizeof(int));
  for(int i=0;i<k_LSH;i++){
    int temp = (rand()%(10)) + 1;
    gfun->r[i]=temp;
  }
  gfun->m=(INT_MAX-4);
}

void destroyG(g_function g){
  for(int i=0;i<k_LSH;i++){
    destroyH_LSH(g.h_functions[i]);
  }
  free(g.h_functions);
  free(g.r);
}

int computeG(g_function gfun,Vector p,unsigned int *id){
  long long int sum = 0;
  // g(p) = [(r1h1(p) + r2h2(p) + · · · + rkhk (p)) mod M] mod TableSize
  // g is a random combination of hi's, every g function has k_LSH h functions
  // compute g function and at the same time compute all the h functions that make it up
  // with this property of modulus operator overflow is avoided ((a _ b) mod m = ((a mod m)(b mod m)) mod m)
  for(int i=0;i<k_LSH;i++){
    int h = computeH_LSH(gfun.h_functions[i],p);
    sum += mod_LLI_UI(gfun.r[i]*h,gfun.m);
  }

  int temp_ID = mod_LLI_I(sum,gfun.m);


  // Store object ID along with pointer to object (Querying trick), for all bucket elements to avoid to compute the Euclidean distance for all vectors p in bucket
  // do it only for p which: ID(p) = ID(q)
  (*id) = temp_ID;
  return mod_Int_Int(temp_ID,hashTableSize);
}


/* LSH IMPLEMENTATION*/

LSH initializeLSH(int l){
  LSH tempLSH = malloc(sizeof(lshNode));
  // allocate as many G functions as the number of hash tables (g functions used like hash functions at the corresponding  hash tables )
  tempLSH->g_fun = malloc(l*sizeof(g_function));
  // allocate the hash tables that the LSH need
  tempLSH->hts = malloc(l*sizeof(HashTable));
  printf("-HASHSIZE = %d\n",hashTableSize);
  // generate G functions and initialize the correspodings hash tables
  for(int i=0;i<l;i++){
     generateG(&(tempLSH->g_fun[i]));
     tempLSH->hts[i] = htInitialize(hashTableSize);
  }
  // save l (the number of hash tables)
  tempLSH->l=l;
  // and finally return the LSH node
  return tempLSH;
}

void insertToLSH(LSH lsh,Vector v){
  // insert the given vector in all LSΗ hash tables
  // the bucket of the hash table that the vector will be inserted depends from the corresponding g function of the specific hash Table (hash function)
  // at the new node tha will be inserted at the hash Tables save the id (Querying trick)
  int l = lsh->l;
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    unsigned int id;
    int index = computeG(lsh->g_fun[i],v,&id); // compute the value of the g function for the given vector that will be inserted
    // finally insert the vector at the corresponding bucket of the current hash table
    htInsert(lsh->hts[i],v,index,id);
  }
}

void insertTimeSeriesToLSH(LSH lsh,Grids grids,Vector timeVector,double delta,Vector v){
  // insert the given vector in all LSΗ hash tables
  // the bucket of the hash table that the vector will be inserted depends from the corresponding g function of the specific hash Table (hash function)
  // at the new node tha will be inserted at the hash Tables save the id (Querying trick)
  int l = lsh->l;
  // printf("*** FOR VECTOR WITH ID ");
  // printVectorId(v);
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    // printf("----- L = %d -------\n",i);
    double t_of_grid = getTofGrid(grids,i);
    // printf("*-*  T = %f\n",t_of_grid);
    // printf("ORIGINAL = ");
    // printVector(v);
    Vector snappedToGrid = timeSeriesSnapping(v,timeVector,d,delta,t_of_grid);
    // d*=2;
    // printf("SNAPPED = ");
    // printVector(snappedToGrid);
    unsigned int id;
    int index = computeG(lsh->g_fun[i],snappedToGrid,&id); // compute the value of the g function for the given vector that will be inserted
    // finally insert the vector at the corresponding bucket of the current hash table
    // d/=2;
    htInsert(lsh->hts[i],v,index,id);
    d*=2;
    deleteVector(snappedToGrid);
    d/=2;
  }
  // getchar();
}


void insertContinuousTimeSeriesToLSH(LSH lsh,Grids grids,Vector timeVector,double delta,Vector v,double epsilon){
  // insert the given vector in all LSΗ hash tables
  // the bucket of the hash table that the vector will be inserted depends from the corresponding g function of the specific hash Table (hash function)
  // at the new node tha will be inserted at the hash Tables save the id (Querying trick)
  // printf("*** FOR VECTOR WITH ID ");
  // printVectorId(v);
  // printf("ORIGINAL = ");
  // printVector(v);
  Vector v2 = filtering(v,epsilon);
  // printf("FILTERED = ");
  // printVector(v2);
  // printf("DELTA = %f\n",delta);
  Vector v3 = continuousTimeSeriesSnapping(v2,d,delta);
  // printf("SNAPPED = ");
  // printVector(v3);
  Vector v4 = minima_maxima(v3,d);
  // printf("KEY_VECTOR = ");
  // printVector(v4);


  int l = lsh->l;
  if(l!=1){ // TODO: REMOVE
    printf("L!=1 ON CONTINUOUS LSH\n");

  }
  unsigned int id;
  int index = computeG(lsh->g_fun[0],v4,&id); // compute the value of the g function for the given vector that will be inserted
  // finally insert the vector at the corresponding bucket of the current hash table
  htInsert(lsh->hts[0],v,index,id);

  // TODO: FREE VECTORS
}

void insertFromListToLSH(List list,LSH lsh){
  // insert every vector of the list at the corresponding LSH
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
      insertToLSH(lsh,getVector(temp));
      temp=getNext(temp);
  }
}

void insertTimeSeriesFromListToLSH(List list,LSH lsh,Grids grids,Vector timeVector,double delta){
  // insert every vector of the list at the corresponding LSH
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
      insertTimeSeriesToLSH(lsh,grids,timeVector,delta,getVector(temp));
      temp=getNext(temp);
  }
}

void insertContinuousTimeSeriesFromListToLSH(List list,LSH lsh,Grids grids,Vector timeVector,double delta,double epsilon){
  // insert every vector of the list at the corresponding LSH
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
      insertContinuousTimeSeriesToLSH(lsh,grids,timeVector,delta,getVector(temp),epsilon);
      temp=getNext(temp);
  }
}

void printLSH(LSH lsh){
  // just print the LSH
  int l = lsh->l;
  for(int i=0;i<l;i++){
    printf("-------- HASH %d --------\n",i+1);
    htPrint(lsh->hts[i]);
    printf("--------------------------\n");
  }
}

void destroyLSH(LSH lsh){
  for(int i=0;i<lsh->l;i++){
     destroyG(lsh->g_fun[i]);
     htDelete(lsh->hts[i],!i);
  }
  free(lsh->g_fun);
  free(lsh->hts);
  free(lsh);
}

void nearestNeigborLSH(LSH lsh,Vector q,Vector *nNearest,double *trueDist,FILE *fptr){
  // find the nearest neighbor of the given vector q with the help of LSH
  Vector nearest=NULL;
  double nearestDist=-1;
  int l = getL(lsh);
  HashTable *hts = getHts(lsh);
  g_function *gfuns = getGfuns(lsh);
  // to find the nearest neighbor of the given vector q, euclidean distance must be applied between the vector q and the vectors of the corresponding bucket of every LSH hash table
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    unsigned int q_ID;
    int q_index = computeG(gfuns[i],q,&q_ID); // compute the value of the g function for the given vector
    // and go to the corresponding bucket of the current hash table to find the  nearest neighbor for the given query vector
    htFindNearestNeighbor(hts[i],q_index,q,&nearest,&nearestDist,d,q_ID);
  }
  // check if nearest neighbor of the given vector q found or not
  if(nearestDist>=0 && nearest!=NULL){
    fprintf(fptr,"Approximate Nearest neighbor: ");
    printVectorIdInFile(nearest,fptr);
    fprintf(fptr,"True Nearest neighbor: ");
    printVectorIdInFile(*nNearest,fptr);
    fprintf(fptr,"distanceApproximate: %f\n",nearestDist);
    fprintf(fptr,"distanceTrue: %f\n", *trueDist);
  }else{
    fprintf(fptr,"- DID NOT FIND NEAREST NEIGHBOR\n");
  }
}

void nearestNeigborLSH_DiscreteFrechet(LSH lsh,Vector q,Vector *nNearest,double *trueDist,FILE *fptr,Grids grids,Vector timeVector,double delta){
  // find the nearest neighbor of the given vector q with the help of LSH
  Vector nearest=NULL;
  double nearestDist=-1;
  int l = getL(lsh);
  HashTable *hts = getHts(lsh);
  g_function *gfuns = getGfuns(lsh);
  // to find the nearest neighbor of the given vector q, euclidean distance must be applied between the vector q and the vectors of the corresponding bucket of every LSH hash table
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    double t_of_grid = getTofGrid(grids,i);
    Vector snappedToGrid = timeSeriesSnapping(q,timeVector,d,delta,t_of_grid);
    unsigned int q_ID;
    // d*=2;
    int q_index = computeG(gfuns[i],snappedToGrid,&q_ID); // compute the value of the g function for the given vector
    // d/=2;
    // and go to the corresponding bucket of the current hash table to find the  nearest neighbor for the given query vector
    htFindNearestNeighbor(hts[i],q_index,q,&nearest,&nearestDist,d,q_ID);
    deleteVector(snappedToGrid);
  }
  // check if nearest neighbor of the given vector q found or not
  if(nearestDist>=0 && nearest!=NULL){
    fprintf(fptr,"Approximate Nearest neighbor: ");
    printVectorIdInFile(nearest,fptr);
    fprintf(fptr,"True Nearest neighbor: ");
    printVectorIdInFile(*nNearest,fptr);
    fprintf(fptr,"distanceApproximate: %f\n",nearestDist);
    fprintf(fptr,"distanceTrue: %f\n", *trueDist);
  }else{
    fprintf(fptr,"- DID NOT FIND NEAREST NEIGHBOR\n");
  }
}

void nearestNeigborLSH_ContinuousFrechet(LSH lsh,Vector q,Vector *nNearest,double *trueDist,FILE *fptr,Vector timeVector,double delta,double epsilon){
  // find the nearest neighbor of the given vector q with the help of LSH
  Vector nearest=NULL;
  double nearestDist=-1;
  int l = getL(lsh);
  if(l!=1){printf("L!=1 ON CONTINUOUS LSH\n");} // TODO: REMOVE
  HashTable *hts = getHts(lsh);

  Vector v2 = filtering(q,epsilon);
  Vector v3 = continuousTimeSeriesSnapping(v2,d,delta);
  Vector v4 = minima_maxima(v3,d);

  unsigned int q_id;
  int q_index = computeG(lsh->g_fun[0],v4,&q_id); // compute the value of the g function for the given vector that will be inserted

  htFindNearestNeighbor(hts[0],q_index,q,&nearest,&nearestDist,d,q_id);

  // check if nearest neighbor of the given vector q found or not
  if(nearestDist>=0 && nearest!=NULL){
    fprintf(fptr,"Approximate Nearest neighbor: ");
    printVectorIdInFile(nearest,fptr);
    fprintf(fptr,"True Nearest neighbor: ");
    printVectorIdInFile(*nNearest,fptr);
    fprintf(fptr,"distanceApproximate: %f\n",nearestDist);
    fprintf(fptr,"distanceTrue: %f\n", *trueDist);
  }else{
    fprintf(fptr,"- DID NOT FIND NEAREST NEIGHBOR\n");
  }
}

void kNearestNeighborsLSH(LSH lsh,Vector q,int knn,double *knearestTrueDists,FILE* fptr){
  // find the k nearest neighbours of the given vector q with the help of LSH
  Vector nearest[knn]; // an array to save the  k nearest neighbours (vectors)
  double knearestDists[knn]; // an array to save the distances of k nearest neighbours
  for (int i = 0; i < knn; i++){
    knearestDists[i]=-1;
    nearest[i]=NULL;
  }
  // to find the k nearest neighbours of the given vector q, euclidean distance must be applied between the vector q and the vectors of the corresponding bucket of every LSH hash table
  int l = getL(lsh);
  HashTable *hts = getHts(lsh);
  g_function *gfuns = getGfuns(lsh);
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    unsigned int q_ID;
    int q_index = computeG(gfuns[i],q,&q_ID); // compute the value of the g function for the given vector
    // and go to the corresponding bucket of the current hash table to find the k nearest neighbors for the given query vector
    htKFindNearestNeighbors(hts[i], q_index, q, nearest, knearestDists, d,knn,q_ID);
  }
  int flag=1;
  for (int i = knn-1; i >= 0; i--){
    // check if k nearest neighbor of the given vector q found or not
    if (knearestDists[i] >= 0 && nearest[i] != NULL){
      fprintf(fptr,"Nearest neighbor-%d: ",knn-i);
      printVectorIdInFile(nearest[i],fptr);
      fprintf(fptr,"distanceLSH: %f\n", knearestDists[i]);
      fprintf(fptr,"distanceTrue: %f\n", knearestTrueDists[i]);
      flag=0;
    }
  }
  if(flag){
    fprintf(fptr,"- DID NOT FIND NEAREST NEIGHBOURS\n");
  }
}

void radiusNeigborsLSH(LSH lsh,Vector q,double radius,FILE *fptr){
  // find the neighbours of the given vector q inside the given radius with the help of LSH
  // store adjacent vectors in a hash table
  int vecsInRadius_size=1;
  int l = getL(lsh);
  HashTable *hts = getHts(lsh);
  g_function *gfuns = getGfuns(lsh);
  if(l>0){
    vecsInRadius_size = getNumberOfVectors(hts[0])/8;
  }
  HashTable vecsInRadius = htInitialize(vecsInRadius_size);  // hash table to store the adjacent vectors of the given vector q
  // to find the neighbours of the given vector q, euclidean distance must be applied between the vector q and the vectors of the corresponding bucket that are inside the radius of every LSH hash table
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    unsigned int q_ID;
    int q_index = computeG(gfuns[i],q,&q_ID); // compute the value of the g function for the given vector
    // and go to the corresponding bucket of the current hash table to do the range search (to find the neighbors of the query vector inside the given radius)
    htFindNeighborsInRadius(hts[i],q_index,vecsInRadius,q,d,q_ID,radius);
  }
  htRangePrint(vecsInRadius,q,d,fptr);

  htDelete(vecsInRadius,0);
}
void radiusNeigborsClustering(LSH lsh,Vector q,double radius,HashTable vecsInRadius,int centroidIndex,List* confList,int *assignCounter,int iteration){
  // based on the given centroids find the clusters that the given vectors belong with the help of LSH (this function used for the "reverseAssignmentLSH")
  // the clusters are represented by hash tables
  int l = getL(lsh);
  HashTable *hts = getHts(lsh);
  g_function *gfuns = getGfuns(lsh);
  for(int i=0;i<l;i++){ // go at every hash table of lsh
    unsigned int q_ID;
    int q_index = computeG(gfuns[i],q,&q_ID); // compute the value of the g function for the given vector
    // and go to the corresponding bucket of the current hash table to do the range search (to find the vectors that belong to the corresponding cluster)
    htFindNeighborsInRadiusClustering(hts[i],q_index,centroidIndex,confList,vecsInRadius,q,d,q_ID,radius,assignCounter,iteration);
  }
}
