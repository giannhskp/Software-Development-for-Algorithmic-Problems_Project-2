#ifndef LSH_H
#define LSH_H

typedef struct lsh_n * LSH;
typedef struct listNode *List;
typedef struct grid_n *Grids;

LSH initializeLSH(int );
void insertToLSH(LSH ,Vector );
void insertFromListToLSH(List ,LSH );
void insertTimeSeriesFromListToLSH(List ,LSH ,Grids ,Vector ,double );
void printLSH(LSH );
void destroyLSH(LSH );

void nearestNeigborLSH(LSH ,Vector,double *,FILE* );
void kNearestNeighborsLSH(LSH, Vector,int,double *,FILE*);
void radiusNeigborsLSH(LSH ,Vector ,double,FILE* );
void radiusNeigborsClustering(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int );



Vector timeSeriesSnapping(Vector,Vector ,int ,double ,double );
Grids initializeGrids(double ,int );
void deleteGrids(Grids );
double getTofGrid(Grids ,int );
void nearestNeigborLSH_DiscreteFrechet(LSH ,Vector ,double *,FILE *,Grids ,Vector ,double );
#endif
