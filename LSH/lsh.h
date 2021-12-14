#ifndef LSH_H
#define LSH_H

typedef struct lsh_n * LSH;
typedef struct listNode *List;
typedef struct grid_n *Grids;

LSH initializeLSH(int,int );
void insertToLSH(LSH ,Vector );
void insertTimeSeriesToLSH(LSH ,Grids ,double ,Vector );
void insertFromListToLSH(List ,LSH );
void insertTimeSeriesFromListToLSH(List ,LSH ,Grids ,double );
void insertContinuousTimeSeriesFromListToLSH(List ,LSH ,double ,double ,Grids );
void printLSH(LSH );
void destroyLSH(LSH );

void nearestNeigborLSH(LSH ,Vector,Vector *,double *,FILE* );
void kNearestNeighborsLSH(LSH, Vector,int,double *,FILE*);
void radiusNeigborsLSH(LSH ,Vector ,double,FILE* );
void radiusNeigborsClustering(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int );
void radiusNeigborsClusteringTimeSeries(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int ,Grids ,double );



Vector timeSeriesSnapping(Vector,Vector ,int ,double ,double );
Grids initializeGrids(double ,int );
void deleteGrids(Grids );
double getTofGrid(Grids ,int );
void nearestNeigborLSH_DiscreteFrechet(LSH ,Vector ,Vector *,double *,FILE *,Grids ,double );
void nearestNeigborLSH_ContinuousFrechet(LSH ,Vector ,Vector *,double *,FILE * ,double ,double ,Grids );

int getValueOfFirstGFun(LSH ,Vector ,unsigned int * );
#endif
