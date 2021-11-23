#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "../../Vector/vector.h"
#include "../hashTable.h"


#define SQUARE(x) ((x)*(x))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


double distance_metric(Vector v1,Vector v2,int d){
  // calculate the Euclidean distance (or L2) between the given vectors and return it
  double sum = 0.0;
  double *coords1 = getCoords(v1);
  double *coords2 = getCoords(v2);
  for(int i=0;i<d;i++){
    sum += SQUARE(coords1[i]-coords2[i]);
  }
  return sqrt(sum);
}


struct listNode {
  Vector v;
  int vector_ID; // ID used to LSH for th Querying trick, (at the other implementations ID initialized with -1 )
  List next;
};
typedef struct listNode *List;

Vector getVector(List l ){
  return l->v;
}
List getNext(List l){
  return l->next;
}


List initializeList(){
  return NULL;
}


List allocateListNode(Vector v,int id){
  List node=malloc(sizeof(struct listNode));
  node->v=v;
  node->vector_ID=id;
  node->next=NULL;
  return node;
}



List listInsert(List list,Vector v,int id){
  if (list==NULL){
    // list is empty
    List node=allocateListNode(v,id);
    return node;
  }
  List newnode = allocateListNode(v,id);
  newnode->next = list;

  return newnode;
}

List listUniqueInsert(List list,Vector v,int id){
  if (list==NULL){
    // list is empty
    List node=allocateListNode(v,id);
    return node;
  }
  List origin = list;
  while(list!=NULL){
    if(compareVectors(v, list->v)){
      return origin;
    }
    list=list->next;
  }

  List newnode = allocateListNode(v,id);
  newnode->next = origin;

  return newnode;
}

List listDeleteItem(List list,Vector v,int id){
  if (list==NULL){
    // list is empty
    return list;
  }
  List origin = list;
  List prev = NULL;
  while(list!=NULL){
    if(compareVectors(v, list->v)){
      if(prev==NULL){
        origin = list->next;
      }else{
        prev->next = list->next;
      }
      free(list);
      return origin;
    }
    prev = list;
    list=list->next;
  }

  return origin;
}


void listPrint(List list){
  int count =0;
    if(list==NULL){ printf("List Empty!\n");  return;}
    List temp=list;
    while(temp!=NULL){
        // printVectorId(temp->v);
        count++;
        temp=temp->next;
    }
    printf(" %d\n",count);
}

void listPrintClusteringInFile(List list,FILE* fptr){
    if(list==NULL){ return;}
    List temp=list;
    while(temp!=NULL){
        printVectorIdInFileNoNewline(temp->v,fptr);
        if(temp->next!=NULL)
          fprintf(fptr,", ");
        temp=temp->next;
    }
}


void listRangePrint(List list,Vector q,int d,int *counter,FILE *fptr){
    if(list==NULL){ return;}
    List temp=list;
    while(temp!=NULL){
        fprintf(fptr,"%d: ",(*counter)++);
        fprintf(fptr,"Found vector: ");
        printVectorIdInFile(temp->v,fptr);
        temp=temp->next;
    }
}


Vector *transformListToArray(List list,int size){
  // used to convert the given list to an array of vectors and return it
  Vector *array = malloc(size*sizeof(Vector));
  List temp = list;
  int index=0;
  while(temp!=NULL){
      if(index>=size){
        printf("!! OUT OF RANGE (transformListToArray)\n");
        break;
      }
      array[index++] = temp->v;
      temp=temp->next;
  }
  return array;
}


List listDelete(List list,int freeVectors){
    // delete whole list
    if(list==NULL) return NULL;

    List current = list;
    List next;

    while(current!=NULL){
        next=current->next;
        if(freeVectors){
          deleteVector(current->v);
        }
        free(current);
        current=next;
    }
    list=NULL;
    return list;
}

void listFindNearestNeighbor(List list,Vector q,Vector *nearest,double *nearestDist,int d,int id){
  //  used to LSH to find the nearest neighbor and the corresponding distance for the given vector q at a list
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
    if(id==(temp->vector_ID)){ // (Querying trick, to avoid to compute Euclidean distance for all vectors in bucket)
      double dist = distance_metric(temp->v,q,d);
      if(dist<(*nearestDist) || (*nearestDist)<0){
        (*nearestDist) = dist;
        (*nearest) = temp->v;
      }
    }
    temp=temp->next;
  }
}

void listFindNearestNeighborCube(List list,Vector q,Vector *nearest,double *nearestDist,int d,int *numOfSearched,int maxToSearch){
  // used to Hypercube to find the nearest neighbor and the corresponding distance for the given vector q at a list
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
    if((*numOfSearched)>=maxToSearch){
      return;
    }
    double dist = distance_metric(temp->v,q,d);
    if(dist<(*nearestDist) || (*nearestDist)<0){
      (*nearestDist) = dist;
      (*nearest) = temp->v;
    }
    (*numOfSearched)+=1; // increase the number of the searched vectors
    temp=temp->next;
  }
}

void swapDoubles(double *xp, double *yp)
{
  double temp = *xp;
  *xp = *yp;
  *yp = temp;
}
void swapVectors(Vector *xp, Vector *yp)
{
  Vector temp = *xp;
  *xp = *yp;
  *yp = temp;
}

int partition(double arr[], Vector *nearest, int low, int high){
  // used at quickSort to find the partition position
  double pivot = arr[high];
  int i = (low - 1);
  for (int j = low; j < high; j++) {
    if (arr[j] >= pivot) {

      i++;

      swapDoubles(&arr[i], &arr[j]);
      swapVectors(&nearest[i], &nearest[j]);
    }
  }

  swapDoubles(&arr[i + 1], &arr[high]);
  swapVectors(&nearest[i + 1], &nearest[high]);

  return (i + 1);
}

void quickSort(double arr[], Vector *nearest, int low, int high){
  if (low < high) {

    int pi = partition(arr, nearest, low, high);

    quickSort(arr,nearest, low, pi - 1);

    quickSort(arr, nearest, pi + 1, high);
  }
}

int binarySearch(double arr[], int l, int r, double x){
  if (r >= l){
    int mid = l + (r - l) / 2;
    if (arr[mid] == x)
      return mid;

    if (arr[mid] < x)
      return binarySearch(arr, l, mid - 1, x);

    return binarySearch(arr, mid + 1, r, x);
  }

  return -1;
}

void listFindKNearestNeighbors(List list,Vector q,Vector *nearest,double *nearestDist,int d,int k,int id){
  // used to LSH to find the k nearest neighbors and the corresponding distance for the given vector q at a list
  if(list==NULL){ return;}
  List temp=list;
  int filledFlag=0; // the array of k nearest neighbors initialize with NULL (nearest[i]) and "-1" value (nearestDist[i]) respectively
  // this flag used to be sure that these two array will be take/initialize by real distances and vectors (then the filledFlag becomes "1")
  while(temp!=NULL){
      if(id!=(temp->vector_ID)){ // (Querying trick, to avoid to compute Euclidean distance for all vectors in bucket)
        temp = temp->next;
        continue;
      }
      int flag = 1;   // this flag also used to be sure that these two array will be take/initialize by real distances and vectors (when this flag does not change value from 1 to 0, it means that the corresponding arrays are initialized)
      int added=0;  // used to check if a vector added at the corresponding arrays to do quickSort
      double dist = distance_metric(temp->v,q,d);


      // check if the vector is alreast exists at the array of k nearest neighbors
      int index=binarySearch(nearestDist,0,k-1,dist);
      if(index!=-1 && nearest[index] != NULL && compareVectors(nearest[index], temp->v)){
        // vector is alreast exists, go to the next one
        temp = temp->next;
        continue;
      }
      if(!filledFlag) // check if the array still has its initial values
        for (int i = 0; i < k; i++){
          if(nearestDist[i]<0){
            flag=0;
            added=1;
            nearestDist[i]=dist;
            nearest[i]=temp->v;
            break;
          }
        }

      if (flag){
        filledFlag=1;
        if(dist<nearestDist[0]){ // array is sorted in descending order, compare the new distance with the first one distance of the array only
          // replace it
          added = 1;
          nearestDist[0] = dist;
          nearest[0] = temp->v;
        }
      }

      if(added) // check if a vector added at the corresponding arrays to do quickSort
        quickSort(nearestDist, nearest,0, k-1);
      temp = temp->next;
  }
}

void listFindKNearestNeighborsCube(List list,Vector q,Vector *nearest,double *nearestDist,int d,int k,int *numOfSearched,int maxToSearch){
    // used to Hypercube to find the k nearest neighbors and the corresponding distance for the given vector q at a list
  if(list==NULL){ return;}
  List temp=list;
  int filledFlag=0;// the array of k nearest neighbors initialize with NULL (nearest[i]) and "-1" value (nearestDist[i]) respectively
  // this flag used to be sure that these two array will be take/initialize by real distances and vectors (then the filledFlag becomes "1")
  while(temp!=NULL){
      if((*numOfSearched)>=maxToSearch){ // number of vector searches reached, stop searching
        return;
      }
      int flag = 1; // this flag also used to be sure that these two array will be take/initialize by real distances and vectors (when this flag does not change value from 1 to 0, it means that the corresponding arrays are initialized)
      int added=0; // used to check if a vectoradded at the corresponding arrays to do quickSort
      double dist = distance_metric(temp->v,q,d);
      (*numOfSearched) += 1;

      // check if the vector is alreast exists at the array of k nearest neighbors
      int index=binarySearch(nearestDist,0,k-1,dist);
      if(index!=-1 && nearest[index] != NULL && compareVectors(nearest[index], temp->v)){
        // vector is alreast exists, go to the next one
        temp = temp->next;
        continue;
      }
      if(!filledFlag) // check if the array still has its initial values
        for (int i = 0; i < k; i++){
          if(nearestDist[i]<0){
            flag=0;
            added=1;
            nearestDist[i]=dist;
            nearest[i]=temp->v;
            break;
          }
        }

      if (flag){
        filledFlag=1;
        if(dist<nearestDist[0]){ // array is sorted in descending order, compare the new distance with the first one distance of the array only
          // replace it
          added = 1;
          nearestDist[0] = dist;
          nearest[0] = temp->v;
        }
      }
      if(added) // check if a vector added at the corresponding arrays to do quickSort
        quickSort(nearestDist, nearest,0, k-1);
      temp = temp->next;
  }
}



void listFindNeighborsInRadius(List list,HashTable storeNeighbors,Vector q,int d,int id,int radius){
  // (LSH) used to find the neighbors of the given vector q inside at the given radius
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
    if(id==(temp->vector_ID)){ // (Querying trick, to avoid to compute Euclidean distance for all vectors in bucket)
      double dist = distance_metric(temp->v,q,d);
      if(dist<=radius){
        htRangeInsert(storeNeighbors,temp->v,temp->vector_ID,d);
      }
    }
    temp=temp->next;
  }
}


void listFindNeighborsInRadiusClustering(List list,int centroidIndex,List* confList,HashTable storeNeighbors,Vector q,int d,int id,int radius,int *assignCounter,int iteration){
  // used at reverseAssignmentLSH to assign the vectors in to the clusters (clusters represented with hash tables)
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
    if(id==(temp->vector_ID)){ // (Querying trick, to avoid to compute Euclidean distance for all vectors in bucket)
        // check if vector has already been assigned at the same iteration in one cluster
        if(assignedToCluster(temp->v) && (getAssignedIteration(temp->v)==iteration)){
          int assignedCluster = getAssignedCluster(temp->v);
          // then
          // check if vector has already been assigned at the some cluster (check the centroid index)
          // or if vector has already been assigned previously in cluster at a search with different radius
          if(assignedCluster==centroidIndex || (((int)getAssignedAtRadius(temp->v))!=radius)){
            // ok, then skip it
            temp=temp->next;
            continue;
          }
          else{
            double dist = distance_metric(temp->v,q,d);
            if(dist<=radius){ // check if for the given radius the current vector lies in ≥ 2 clusters, add it into the confList
              *confList=listInsert(*confList,temp->v,-1);
            }
            temp=temp->next;
            continue;
          }
        }
        else{
          double dist = distance_metric(temp->v,q,d);
          if(dist<=radius){
            // insert the vector at the corresponding cluster (hash table)
            htRangeInsert(storeNeighbors,temp->v,temp->vector_ID,d);
            setAssignedCluster(temp->v,centroidIndex);
            setAssignedIteration(temp->v,iteration);
            setAssignedAtRadius(temp->v,radius);
            (*assignCounter)++;
          }
        }

    }
    temp=temp->next;
  }
}

void listFindNeighborsInRadiusClusteringCube(List list,int centroidIndex,List* confList,HashTable storeNeighbors,Vector q,int d,double radius,int *numOfSearched,int maxToSearch,int *assignCounter,int iteration){
  // used at reverseAssignmentHypercube to assign the vectors in to the clusters (clusters represented with hash tables)
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
      if((*numOfSearched)>=maxToSearch){
        return;
      }
      (*numOfSearched)++;
      // check if vector has already been assigned at the same iteration in one cluster
      if(assignedToCluster(temp->v) && (getAssignedIteration(temp->v)==iteration)){
        int assignedCluster = getAssignedCluster(temp->v);
        // then
        // check if vector has already been assigned at the some cluster (check the centroid index)
        // or if vector has already been assigned previously in cluster at a search with different radius
        if(assignedCluster==centroidIndex || (getAssignedAtRadius(temp->v)!=radius )){
          // ok, then skip it
          temp=temp->next;
          continue;
        }
        else{
          double dist = distance_metric(temp->v,q,d);
          if(dist<=radius){ // check if for the given radius the current vector lies in ≥ 2 clusters, add it into the confList
            *confList=listInsert(*confList,temp->v,-1);
          }
          temp=temp->next;
          continue;
        }
      }
      else{
        double dist = distance_metric(temp->v,q,d);
        if(dist<=radius){
          // insert the vector at the corresponding cluster (hash table)
          htRangeInsert(storeNeighbors,temp->v,temp->vector_ID,d);
          setAssignedCluster(temp->v,centroidIndex);
          setAssignedIteration(temp->v,iteration);
          setAssignedAtRadius(temp->v,radius);
          (*assignCounter)++;
        }
      }
    temp=temp->next;
    }
}

void listSolveRangeConflicts(List conflictList,HashTable *clustersHt,Vector *clusters,int numOfClusters,int d,int iteration){
  // used to solve the conflicts for all vectors that has been inserted in the given list
  // more specifically these vectors has been assigned in more than one clusters previously at range search,
  // so let's find the cluster that really belongs these vectors with the help of Euclidean Distance
  if(conflictList==NULL){ return;}
  List temp=conflictList;
  int count = 0;
  int countReal =0;
  while(temp!=NULL){
    count++;
    // check if the conflict vector has already been "solved" (at this case negative radius has been assigned)
    if(getAssignedAtRadius(temp->v)<-1.0){
      // then skip it
      temp=temp->next;
      continue;
    }
    // let's solve the conflict and assign the vector at the "correct" cluster based on euclidean distance
    countReal++;
    // first delete it from the already assigned cluster
    htRangeDelete(clustersHt[getAssignedCluster(temp->v)],temp->v,temp->vector_ID,d);
    double minDist = DBL_MAX;
    int closestCentroid = -1;
    // and find the closest centroid based on euclidean distance
    for(int i=0;i<numOfClusters;i++){
      double dist = distance_metric(temp->v,clusters[i],d);
      if(dist<minDist){
        minDist = dist;
        closestCentroid = i;
      }
    }
    // closest centroid found, insert the vector at the corresponding cluster
    htRangeInsert(clustersHt[closestCentroid],temp->v,temp->vector_ID,d);
    setAssignedCluster(temp->v,closestCentroid);
    setAssignedAtRadius(temp->v,-3.0);
    setAssignedIteration(temp->v,iteration);
    temp=temp->next;
    // go to the next one
  }

}

void listFindNeighborsInRadiusCube(List list,HashTable storeNeighbors,Vector q,int d,int radius,int *numOfSearched,int maxToSearch){
  if(list==NULL){ return;}
  List temp=list;
  while(temp!=NULL){
    if((*numOfSearched)>=maxToSearch){
      return;
    }
    double dist = distance_metric(temp->v,q,d);
    (*numOfSearched) += 1;
    if(dist<=radius){
      htRangeInsert(storeNeighbors,temp->v,temp->vector_ID,d);
    }
    temp=temp->next;
  }
}


Vector listMeanOfCluster(List list,int d){
  // used to find/calculate the centroid of the cluster (cluster is the given list)
  if(list==NULL) return NULL;
  List temp = list;
  int count=0;
  double *sumDims=calloc(d,sizeof(double));
  while(temp!=NULL){
    for(int i=0;i<d;i++){
      sumDims[i]+=getCoords(temp->v)[i]; // add up all the coordinates of the cluster vectors
    }
    count++;
    temp=temp->next;
  }
  for(int i=0;i<d;i++){
    sumDims[i]/=(double)count; // divide the coordinates with their number
  }
  // finally create the centroid
  Vector newCentroid  = initVector(sumDims,"tempCentroid");
  free(sumDims);
  // and return it
  return newCentroid;
}

double *listSumOfVectors(List list,int d,int *count){
  // used to calculate the sum of vectors at one cluster (cluster is the given list)
  if(list==NULL) return NULL;
  List temp = list;
  double *sumDims=calloc(d,sizeof(double));
  while(temp!=NULL){
    for(int i=0;i<d;i++){
      sumDims[i]+=getCoords(temp->v)[i];
    }
    (*count)++;
    temp=temp->next;
  }
  return sumDims;
}


int closestCentroid(Vector v,Vector *clusters,int numOfClusters,int d,int except){
  // used to find for the given vector the closest centroid based on euclidean distance
  int minDistIndex = -1;
  double minDist = DBL_MAX;
  for(int i=0;i<numOfClusters;i++){
    if(i==except) continue;
    double tempDist = distance_metric(v,clusters[i],d);
    if(tempDist<minDist){
      minDistIndex = i;
      minDist = tempDist;
    }
  }
  return minDistIndex;
}

double listFindSumOfDistancesOfVector(List list,Vector v,int *count,int d){
  // used to find the sum of distances from the given vector to vectors of one cluster (cluster is the given list)
  if(list==NULL) return 0.0;
  List temp = list;
  double tempSum = 0.0;
  while(temp!=NULL){
    double dist = distance_metric(temp->v,v,d);
    tempSum += dist;
    (*count)++;
    temp=temp->next;
  }
  return tempSum;
}

void listComputeAverageDistOfEveryPointOfCluster(List *listArray,int arraySize,HashTable *clustersHt,int current_cluster,Vector *clusters,int numOfClusters,double *a,double *b,int *count,int d){
  // used to compute a(i) and b(i) for the silhouette of each cluster at reverseAssignmentLSH and reverseAssignmentHypercube Algorithms
  for(int i=0;i<arraySize;i++){
    List temp = listArray[i];
    while(temp!=NULL){
      // a(i) = average distance of i to objects in same cluster
      a[(*count)]=htFindAverageDistanceOfVectorInCluster(clustersHt[current_cluster],temp->v,d);
      int second_closest_centroid = closestCentroid(temp->v,clusters,numOfClusters,d,current_cluster);
      // b(i) = average distance of i to objects in next best (neighbor) cluster
      b[(*count)]=htFindAverageDistanceOfVectorInCluster(clustersHt[second_closest_centroid],temp->v,d);
      (*count)++;
      temp=temp->next;
    }
  }
}


double listFindAverageDistOfVector(List list,Vector v,int d){
  // used to compute average distance of given vector in the given cluster (list)
  if(list==NULL) return 0.0;
  List temp = list;
  double tempSum = 0.0;
  int count = 0;
  while(temp!=NULL){
    double dist = distance_metric(temp->v,v,d);
    tempSum += dist;
    (count)++;
    temp=temp->next;
  }
  if(count==0) return 0.0;
  return tempSum/count;
}

void listComputeAverageDistOfEveryPointOfClusterLloyds(List list,List *clustersLists,int current_cluster,Vector *clusters,int numOfClusters,double *a,double *b,int *count,int d){
  // used to compute compute a(i) and b(i) for the silhouette of each cluster at Lloyds Algorithm
    List temp = list;
    while(temp!=NULL){
      // a(i) = average distance of i to objects in same cluster
      a[(*count)]=listFindAverageDistOfVector(clustersLists[current_cluster],temp->v,d);
      int second_closest_centroid = closestCentroid(temp->v,clusters,numOfClusters,d,current_cluster);
      // b(i) = average distance of i to objects in next best (neighbor) cluster
      b[(*count)]=listFindAverageDistOfVector(clustersLists[second_closest_centroid],temp->v,d);
      (*count)++;
      temp=temp->next;
    }
}


double silhouetteofClusterLloyds(List *clustersLists,Vector *clusters,int current_cluster,int numOfClusters,int numOfVectorsInCluster,int d,double *stotal){
  // used to find the silhouettes of each clusters (clusters given in a array of lists)
  // a(i) = average distance of i to objects in same cluster
  // b(i) = average distance of i to objects in next best (neighbor) cluster
  double *a = calloc(sizeof(double),numOfVectorsInCluster);
  double *b = calloc(sizeof(double),numOfVectorsInCluster);
  int count = 0;
  listComputeAverageDistOfEveryPointOfClusterLloyds(clustersLists[current_cluster],clustersLists,current_cluster,clusters,numOfClusters,a,b,&count,d);
  double sumOfS_i = 0.0;
  for(int i=0;i<numOfVectorsInCluster;i++){
    double s_i = (b[i]-a[i])/ MAX(b[i],a[i]);
    (*stotal) += s_i;
    sumOfS_i += s_i;
  }
  free(a);
  free(b);
  return sumOfS_i/numOfVectorsInCluster;
}






//
