#ifndef LIST_H
#define LIST_H



typedef struct listNode *List;
typedef struct hashtable_head *HashTable;

Vector getVector(List );
List getNext(List );

double distance_metric(Vector ,Vector ,int );

List initializeList();
List allocateListNode(Vector ,int );
List listInsert(List, Vector,int );
List listUniqueInsert(List ,Vector ,int );
List listSearchId(List ,int );
void listPrint(List );
void listRangePrint(List ,Vector ,int,int *,FILE* );
List listDeleteItem(List ,Vector ,int );
List listDelete(List ,int );
Vector *transformListToArray(List ,int );
void listPrintClusteringInFile(List ,FILE* );

void listFindNearestNeighbor(List ,Vector ,Vector *,double *,int ,int );
void listFindKNearestNeighbors(List, Vector, Vector *, double *, int,int ,int );
void listFindNeighborsInRadius(List ,HashTable ,Vector ,int ,int ,int );
void listFindNeighborsInRadiusClustering(List ,int ,List* ,HashTable ,Vector ,int ,int ,int ,int *,int );
void listFindNeighborsInRadiusClusteringCube(List ,int ,List* ,HashTable ,Vector ,int ,double ,int *,int ,int *,int );

void listFindNearestNeighborCube(List ,Vector ,Vector *,double *,int ,int *,int );
void listFindKNearestNeighborsCube(List ,Vector ,Vector *,double *,int ,int ,int *,int );
void listFindNeighborsInRadiusCube(List ,HashTable ,Vector ,int ,int ,int *,int );

Vector listMeanOfCluster(List ,int );
void listSolveRangeConflicts(List ,HashTable *,Vector *,int ,int ,int );
double *listSumOfVectors(List ,int ,int *);
double listFindSumOfDistancesOfVector(List ,Vector ,int *,int );
void listComputeAverageDistOfEveryPointOfCluster(List *,int ,HashTable *,int ,Vector *,int ,double *,double *,int *,int );
double silhouetteofClusterLloyds(List *,Vector *,int ,int ,int ,int ,double *);
#endif
