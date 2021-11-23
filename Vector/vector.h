#ifndef VECTOR_H
#define VECTOR_H

typedef struct vec_node *Vector;

Vector initVector(double *, char []);

Vector copyVector(Vector );

void deleteVector(Vector);

void printVector(Vector );
void printVectorId(Vector );

void printVectorInFile(Vector,FILE* );
void printVectorIdInFile(Vector,FILE* );
void printVectorIdInFileNoNewline(Vector,FILE* );

int compareVectors(Vector , Vector );

int assignedToCluster(Vector );

int getAssignedCluster(Vector );

void setAssignedCluster(Vector ,int );

double *getCoords(Vector);

void initializeClusterInfo(Vector );
int getAssignedIteration(Vector );
void setAssignedIteration(Vector ,int );
double getAssignedAtRadius(Vector );
void setAssignedAtRadius(Vector ,double );

#endif
