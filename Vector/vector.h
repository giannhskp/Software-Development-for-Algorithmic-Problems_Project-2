#ifndef VECTOR_H
#define VECTOR_H

typedef struct vec_node *Vector;

// Vector initVector(double *, char []);
Vector initVector(double *, char [],int);
// Vector initTimeSeries(double *,double *, char []);
Vector initTimeSeries(double *,double *, char [],int);

Vector copyVector(Vector );

void deleteVector(Vector);

void printVector(Vector );
void printTimes(Vector );
void printVectorId(Vector );

int getDim(Vector );
void printVectorInFile(Vector,FILE* );
void printVectorIdInFile(Vector,FILE* );
void printVectorIdInFileNoNewline(Vector,FILE* );

int compareVectors(Vector , Vector );

int assignedToCluster(Vector );

int getAssignedCluster(Vector );

void setAssignedCluster(Vector ,int );

double *getCoords(Vector);
double* getTime(Vector );
char* getID(Vector );

void initializeClusterInfo(Vector );
int getAssignedIteration(Vector );
void setAssignedIteration(Vector ,int );
double getAssignedAtRadius(Vector );
void setAssignedAtRadius(Vector ,double );

Vector shiftVector(Vector ,double );

int compareTimeSeries(Vector ,Vector);

#endif
