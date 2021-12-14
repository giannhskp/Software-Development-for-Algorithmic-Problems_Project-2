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

void nearestNeigborLSH(LSH ,Vector,Vector *,double *,FILE* ,double *,int *);
void kNearestNeighborsLSH(LSH, Vector,int,double *,FILE*);
void radiusNeigborsLSH(LSH ,Vector ,double,FILE* );
void radiusNeigborsClustering(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int );
void radiusNeigborsClusteringTimeSeries(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int ,Grids ,double ,int );


Grids initializeGrids(double ,int ,int );
void deleteGrids(Grids ,int );
double getTofGrid(Grids ,int ,int);
void nearestNeigborLSH_DiscreteFrechet(LSH ,Vector ,Vector *,double *,FILE *,Grids ,double ,double *,int *);
void nearestNeigborLSH_ContinuousFrechet(LSH ,Vector ,Vector *,double *,FILE * ,double ,double ,Grids ,double *,int *);

int getValueOfFirstGFun(LSH ,Vector ,unsigned int * );
#endif
