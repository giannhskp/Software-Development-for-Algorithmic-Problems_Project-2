#ifndef CLUSTERING_H
#define CLUSTERING_H

int existsInArray(int *,int ,int );
void minDistToCentroids(Vector ,Vector* ,int *,int ,double *);
int* kmeansplusplus(Vector* ,int ,Vector *,double *);
void clustering(List ,FILE*,char * ,int ,int ,int ,int );


#endif
