#ifndef CLUSTERHELPFUNS_H
#define CLUSTERHELPFUNS_H

int existsInArray(int *,int ,int );

void minDistToCentroids(Vector ,Vector* ,Vector *,int ,double *);

void minDistbetweenCentroids(Vector *,int ,double *);

int centroidsConverge(Vector *,Vector *,int ,int );

int findClosestCentroid(Vector ,Vector *,int );

#endif
