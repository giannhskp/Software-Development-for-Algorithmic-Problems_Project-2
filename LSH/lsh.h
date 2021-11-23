#ifndef LSH_H
#define LSH_H

typedef struct lsh_n * LSH;

typedef struct listNode *List;
LSH initializeLSH(int );
void insertToLSH(LSH ,Vector );
void insertFromListToLSH(List ,LSH );
void printLSH(LSH );
void destroyLSH(LSH );

void nearestNeigborLSH(LSH ,Vector,double *,FILE* );
void kNearestNeighborsLSH(LSH, Vector,int,double *,FILE*);
void radiusNeigborsLSH(LSH ,Vector ,double,FILE* );
void radiusNeigborsClustering(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int );
#endif
