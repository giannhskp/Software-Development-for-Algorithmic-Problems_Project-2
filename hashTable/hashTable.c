#include <stdio.h>
#include <stdlib.h>
#include "../Vector/vector.h"
#include "./hashTableList/hashTableList.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))



struct hashtable_node{
  int key;
  List head;
};
typedef struct hashtable_node *hashtable_nodePtr;


struct hashtable_head{
  hashtable_nodePtr table;
  int buckets;
  int numberOfVectors;
};
typedef struct hashtable_head *HashTable;


int getNumberOfVectors(const HashTable ht){
  return ht->numberOfVectors;
}



HashTable htInitialize(int buckets) {
  HashTable ht=malloc(sizeof(struct hashtable_head));
  ht->table=malloc(sizeof(struct hashtable_node)*buckets);
  if(ht->table==NULL){
    printf("MALLOC ERROR : NOT ENOUGH MEMORY TO STORE THE HASH TABLE\n");
    fflush(stdout);
    exit(0);
  }
  ht->buckets=buckets;
  for (int i=0;i<buckets;i++){
    ht->table[i].key=i;
    ht->table[i].head=initializeList();
  }
  ht->numberOfVectors=0; // nodes counter
  return ht;
}



int hashFunction(const HashTable ht,Vector v,int d){ /* Hash function used for hash tables in Radius Search and in clustering (reverseAssignmentLSH and reverseAssignmentHypercube) */
  /* Not used in LSH*/
  /* Based on djb2 hash function*/
  double *coords = getCoords(v);
  double hash = 5381;
  int c;
  for(int i=0;i<d;i++){
    c = coords[i];
    hash = (int)(((hash * 33) + hash) + c)% ht->buckets; /* hash * 33 + c */
  }
  return (((int)hash) % ht->buckets);
}

int htInsert(HashTable ht, Vector v,int index,int id){
  fflush(stdout);
  ht->table[index].head=listInsert(ht->table[index].head,v,id);
  ht->numberOfVectors++;
  return 1;
}

void htRangeInsert(HashTable ht, Vector v,int id,int d){
  int index=hashFunction(ht,v,d);
  ht->table[index].head=listUniqueInsert(ht->table[index].head,v,id);
  ht->numberOfVectors++;
}

void htRangeDelete(HashTable ht, Vector v,int id,int d){
  int index=hashFunction(ht,v,d);
  ht->table[index].head=listDeleteItem(ht->table[index].head,v,id);
  ht->numberOfVectors--;
}


void htPrint(const HashTable ht){
  for (int i=0;i<ht->buckets;i++){
    printf("\n>>Bucket %d: ",i+1);
    listPrint(ht->table[i].head);
  }
  printf("\n\n");
}

void htPrintClustering(const HashTable ht,FILE* fptr){
  for (int i=0;i<ht->buckets;i++){
    listPrintClusteringInFile(ht->table[i].head,fptr);
    if(i!=ht->buckets-1 && (ht->table[i].head!=NULL))
      fprintf(fptr,", ");
  }
  printf("\n\n");
}

void htRangePrint(const HashTable ht,Vector q,int d,FILE *fptr){
  fprintf(fptr,"R-near neighbors:\n");
  int counter=1;
  for (int i=0;i<ht->buckets;i++){
    listRangePrint(ht->table[i].head,q,d,&counter,fptr);
  }
  if(counter==1){
      fprintf(fptr," Did not find any neighbor inside the given range\n");
  }
}


HashTable htDelete(HashTable ht,int freeVectors){
  // delete whole hash table
  for (int i=0;i<ht->buckets;i++){
    ht->table[i].head=listDelete(ht->table[i].head,freeVectors);
  }
  free(ht->table);
  free(ht);
  return NULL;
}

void htFindNearestNeighbor(HashTable ht,int index,Vector q,Vector *nearest,double *nearestDist,int d,int id){
  listFindNearestNeighbor(ht->table[index].head,q,nearest,nearestDist,d,id);
}
void htKFindNearestNeighbors(HashTable ht,int index,Vector q,Vector *nearest,double *nearestDist,int d,int k,int id){
  listFindKNearestNeighbors(ht->table[index].head, q, nearest, nearestDist, d,k,id);
}

void htFindNeighborsInRadius(HashTable ht,int index,HashTable storeNeighbors,Vector q,int d,int id,int radius){
  listFindNeighborsInRadius(ht->table[index].head,storeNeighbors,q,d,id,radius);
}
void htFindNeighborsInRadiusClustering(HashTable ht,int index,int centroidIndex,List* confList,HashTable storeNeighbors,Vector q,int d,int id,int radius,int *assignCounter,int iteration){
  listFindNeighborsInRadiusClustering(ht->table[index].head,centroidIndex,confList,storeNeighbors,q,d,id,radius,assignCounter,iteration);
}

void htFindNeighborsInRadiusClusteringCube(HashTable ht,int index,int centroidIndex,List* confList,HashTable storeNeighbors,Vector q,int d,double radius,int *numOfSearched,int maxToSearch,int *assignCounter,int iteration){
  listFindNeighborsInRadiusClusteringCube(ht->table[index].head,centroidIndex,confList,storeNeighbors,q,d,radius,numOfSearched,maxToSearch,assignCounter,iteration);
}


void htFindNearestNeighborCube(HashTable ht,int index,Vector q,Vector *nearest,double *nearestDist,int d,int *numOfSearched,int maxToSearch){
  listFindNearestNeighborCube(ht->table[index].head,q,nearest,nearestDist,d,numOfSearched,maxToSearch);
}

void htKFindNearestNeighborsCube(HashTable ht,int index,Vector q,Vector *nearest,double *nearestDist,int d,int k,int *numOfSearched,int maxToSearch){
  listFindKNearestNeighborsCube(ht->table[index].head, q, nearest, nearestDist, d,k,numOfSearched,maxToSearch);
}

void htFindNeighborsInRadiusCube(HashTable ht,int index,HashTable storeNeighbors,Vector q,int d,int radius,int *numOfSearched,int maxToSearch){
  listFindNeighborsInRadiusCube(ht->table[index].head,storeNeighbors,q,d,radius,numOfSearched,maxToSearch);
}

Vector htMeanOfCluster(HashTable ht,int d){
  // used to find the centroid of the cluster
  // cluster represented by the hash table
  int count = 0;
  double *sumDims=calloc(d,sizeof(double));
  for (int i=0;i<ht->buckets;i++){
    double *tempSum = listSumOfVectors(ht->table[i].head,d,&count);
    if(tempSum!=NULL){
      for(int i=0;i<d;i++){
        sumDims[i]+=tempSum[i]; // add up all the coordinates of the cluster vectors
      }
      free(tempSum);
    }
  }
  if(count==0){ free(sumDims); return NULL;}
  for(int i=0;i<d;i++){
    sumDims[i]/=(double)count; // divide the coordinates with their number
  }
  // finally create the centroid
  Vector newCentroid  = initVector(sumDims,"temp");
  free(sumDims);
  // and return it
  return newCentroid;
}

double htFindAverageDistanceOfVectorInCluster(HashTable ht,Vector v,int d){
  int count=0;
  double sumOfDists = 0.0;
  for (int i=0;i<ht->buckets;i++){
    sumOfDists += listFindSumOfDistancesOfVector(ht->table[i].head,v,&count,d);
  }
  if(count==0) return 0.0;
  return sumOfDists/count;
}

double silhouetteofClusterLSH(HashTable *clustersHt,Vector *clusters,int current_cluster,int numOfClusters,int d,double *stotal){
  HashTable ht = clustersHt[current_cluster];
  List allBuckets[ht->buckets];
  for (int i=0;i<ht->buckets;i++){
    allBuckets[i]=ht->table[i].head;
  }
  // a(i) = average distance of i to objects in same cluster
  // b(i) = average distance of i to objects in next best (neighbor) cluster
  double *a = calloc(sizeof(double),ht->numberOfVectors);
  double *b = calloc(sizeof(double),ht->numberOfVectors);
  int count = 0;
  listComputeAverageDistOfEveryPointOfCluster(allBuckets,ht->buckets,clustersHt,current_cluster,clusters,numOfClusters,a,b,&count,d);
  double sumOfS_i = 0.0;
  for(int i=0;i<ht->numberOfVectors;i++){
    double s_i = (b[i]-a[i])/ MAX(b[i],a[i]);
    (*stotal) += s_i;
    sumOfS_i += s_i;
  }
  free(a);
  free(b);
  return sumOfS_i/ht->numberOfVectors;
}















//
