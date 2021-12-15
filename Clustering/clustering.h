#ifndef CLUSTERING_H
#define CLUSTERING_H

int existsInArray(int *,int );
void minDistToCentroids(Vector ,Vector* ,int *,int ,double *,int);
int* kmeansplusplus(Vector* ,int ,Vector *,double *,int);
void clustering(List ,FILE* ,char* ,char *,int ,int ,int ,int ,int ,double );


#endif
