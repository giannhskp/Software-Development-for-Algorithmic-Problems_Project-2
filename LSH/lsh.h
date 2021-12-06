#ifndef LSH_H
#define LSH_H

typedef struct lsh_n * LSH;
typedef struct listNode *List;
typedef struct grid_n *Grids;

LSH initializeLSH(int,int );
void insertToLSH(LSH ,Vector );
void insertFromListToLSH(List ,LSH );
void insertTimeSeriesFromListToLSH(List ,LSH ,Grids ,double );
void insertContinuousTimeSeriesFromListToLSH(List ,LSH ,double ,double );
void printLSH(LSH );
void destroyLSH(LSH );

void nearestNeigborLSH(LSH ,Vector,Vector *,double *,FILE* );
void kNearestNeighborsLSH(LSH, Vector,int,double *,FILE*);
void radiusNeigborsLSH(LSH ,Vector ,double,FILE* );
void radiusNeigborsClustering(LSH ,Vector ,double ,HashTable ,int ,List* ,int *,int );



Vector timeSeriesSnapping(Vector,Vector ,int ,double ,double );
Grids initializeGrids(double ,int );
void deleteGrids(Grids );
double getTofGrid(Grids ,int );
void nearestNeigborLSH_DiscreteFrechet(LSH ,Vector ,Vector *,double *,FILE *,Grids ,double );
void nearestNeigborLSH_ContinuousFrechet(LSH ,Vector ,Vector *,double *,FILE * ,double ,double );

// #ifdef __cplusplus
// #define EXTERNC extern "C"
// #else
// #define EXTERNC
// #endif
//
// typedef void* string;
//
// EXTERNC string str();
// EXTERNC void mylibrary_mytype_destroy(mylibrary_mytype_t mytype);
// EXTERNC void mylibrary_mytype_doit(mylibrary_mytype_t self, int param);
//
// #undef EXTERNC

#endif
