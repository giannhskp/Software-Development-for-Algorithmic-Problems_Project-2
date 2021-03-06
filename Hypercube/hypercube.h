#ifndef HYPERCUBE_H
#define HYPERCUBE_H

typedef struct hc_n *HyperCube;

HyperCube initializeHyperCube(int);

void insertToHyperCube(HyperCube ,Vector );
void insertFromListToHyperCube(List ,HyperCube );

void printHyperCube(HyperCube );

void deleteHyperCube(HyperCube );

void nearestNeigborHypercube(HyperCube ,Vector ,Vector *,int ,int,double *,FILE* ,double *,int *,int );

void kNearestNeigborsHypercube(HyperCube ,Vector ,int ,int ,int,double *,FILE* );

void radiusNeigborsHypercube(HyperCube ,Vector ,double ,int ,int,FILE* );
void radiusNeigborHypercubeClustering(HyperCube ,Vector ,HashTable ,double ,int ,int ,int ,List* ,int *,int );

#endif
