#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTable.h"
#include "../LSH/lsh.h"
#include "../hashTable/hashTableList/hashTableList.h"

#define OUT    0
#define IN    1
#define FALSE    0
#define TRUE    1

#define MAX_INPUT_LENGTH 10240

// returns number of words in str
static int countWords(char *str){
    char * token = strtok(str, "	");
    token = strtok(NULL, "	");
     // loop through the string to extract all other tokens
     int counter = 0;
     while( token != NULL ) {
        counter++;
        token = strtok(NULL, "	");
     }
    return counter;
}

int findDimLSH(char* fileName){
  FILE *file = fopen(fileName, "r"); // read mode

  if (file == NULL){
     perror("Error while opening the file.\n");
     exit(-1);
  }

 if (feof(file)){ // empty file, return
   return -1;
 }
 char buffer[MAX_INPUT_LENGTH];
 fflush(stdin);  // clear stdin buffer
 if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
   perror("Error while reading the file.\n");
   exit(-1);
 }
 int dims = countWords(buffer);
 fclose(file);
 return dims-1;
}


void readFileLSH(char* fileName,List *inputs,int *vectorCount,int keepTimes,double *times,int dim){

   FILE *file = fopen(fileName, "r"); // read mode

   if (file == NULL){
      perror("Error while opening the file.\n");
      exit(-1);
   }

  if (feof(file)){ // empty file, return
    return;
  }

  char buffer[MAX_INPUT_LENGTH];

  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }

    double vec[dim];
    char * token = strtok(buffer, "	 ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL, "	");
     // loop through the string to extract all other tokens
     int counter = 0;
     while( token != NULL ) {
        vec[counter++]=atof(token);
        token = strtok(NULL, "	");
     }
     Vector vecTmp;
     if(keepTimes){
       vecTmp=initTimeSeries(vec,times,name,dim);
     }else{
       vecTmp=initVector(vec,name,dim);
     }
     (*inputs) = listInsert((*inputs),vecTmp,-1);
     (*vectorCount)++;
  }

  fclose(file);
}


void readQueryFileLSH(char* queryFile,char* outputFile,LSH lsh,List inputs,int dim,int distanceTrueOff){

   FILE *file = fopen(queryFile, "r"); // read mode

   if (file == NULL){
      perror("Error while opening the file.\n");
      exit(-1);
   }

  if (feof(file)){ // empty file, return
    return;
  }
  FILE* fptr;
  fptr = fopen(outputFile, "w");
  if(fptr == NULL){
    /* File not created hence exit */
    printf("Unable to create file.\n");
    exit(EXIT_FAILURE);
  }

  int n =1 ; // 1 nearest neighbor
  char buffer[MAX_INPUT_LENGTH];
  Vector nNearest[n]; // here store the true k nearest neighbors
  double knearestDists[n]; // here store the true distances from the k nearest neighbors
  double vec[dim];

  int query_count =0;
  double total_lsh_time = 0.0;
  double total_true_time = 0.0;
  double max_aproximation_factor = -1;
  double min_aproximation_factor = DBL_MAX;
  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }


    char * token = strtok(buffer, "	 ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL, "	");
    int counter = 0;
    while( token != NULL ) {
      vec[counter++]=atof(token);
      token = strtok(NULL, "	");
    }
    Vector vecTmp=initVector(vec,name,dim);

    fprintf(fptr, "Query %s:\n",name);

    for (int i = 0; i < n; i++){
      nNearest[i]=NULL;
      knearestDists[i]=-1;
    }

    clock_t begin_true = clock(); // time calculation for the k nearest neighbors with brute force method
    // find with the brute force method the k nearest neighbors for the corresponding query vector
    if(distanceTrueOff!=1)
      listFindNearestNeighbor(inputs,vecTmp,nNearest,knearestDists,dim,-1);

    clock_t end_true = clock();
    double time_spent_true = (double)(end_true - begin_true) / CLOCKS_PER_SEC;

    clock_t begin_lsh = clock(); // time calculation for the k nearest neighbors with LSH
    // find with the help of LSH the k nearest neighbor for the corresponding query vector
    double approximation_factor=-1;
    int found_neighbor = -1;
    nearestNeigborLSH(lsh,vecTmp,nNearest,knearestDists,fptr,&approximation_factor,&found_neighbor,distanceTrueOff);
    if(approximation_factor>0){
      if(approximation_factor>max_aproximation_factor)
        max_aproximation_factor=approximation_factor;
      if(approximation_factor<min_aproximation_factor)
        min_aproximation_factor=approximation_factor;
    }

    clock_t end_lsh = clock();
    double time_spent_lsh = (double)(end_lsh - begin_lsh) / CLOCKS_PER_SEC;
    if(found_neighbor>0){
      total_lsh_time += time_spent_lsh;
      total_true_time += time_spent_true;
      query_count++;
    }
    fprintf(fptr, "\n");
    fflush(fptr);
    deleteVector(vecTmp);
  }
  if(query_count>0){
    fprintf(fptr, "tApproximateAverage: %f seconds\n",total_lsh_time/query_count);
    if(distanceTrueOff!=1){
      fprintf(fptr, "tTrueAverage: %f seconds\n",total_true_time/query_count);
      fprintf(fptr, "Max Approximation Factor: %f\n",max_aproximation_factor);
      fprintf(fptr, "Min Approximation Factor: %f\n",min_aproximation_factor);
    }else{
      fprintf(fptr, "tTrueAverage: True Distances were not computed\n");
      fprintf(fptr, "Min/Max Approximation Factor are not defined without distanceTrue\n");
    }
  }
  fclose(fptr);
  fclose(file);
}

void readQueryFileLSH_DiscreteFrechet(char* queryFile,char* outputFile,LSH lsh,List inputs,Grids grids,double delta,double *time,int dim,int distanceTrueOff){

   FILE *file = fopen(queryFile, "r"); // read mode

   if (file == NULL){
      perror("Error while opening the file.\n");
      exit(-1);
   }

  if (feof(file)){ // empty file, return
    return;
  }
  FILE* fptr;
  fptr = fopen(outputFile, "w");
  if(fptr == NULL){
    /* File not created hence exit */
    printf("Unable to create file.\n");
    exit(EXIT_FAILURE);
  }

  int n =1 ; // 1 nearest neighbor
  char buffer[MAX_INPUT_LENGTH];
  Vector nNearest[n]; // here store the true k nearest neighbors
  double knearestDists[n]; // here store the true distances from the k nearest neighbors
  double vec[dim];

  int query_count =0;
  double total_lsh_time = 0.0;
  double total_true_time = 0.0;
  double max_aproximation_factor = -1;
  double min_aproximation_factor = DBL_MAX;
  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }


    char * token = strtok(buffer, "	 ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL, "	");
    int counter = 0;
    while( token != NULL ) {
      vec[counter++]=atof(token);
      token = strtok(NULL, "	");
    }
    Vector vecTmp=initTimeSeries(vec,time,name,dim);

    fprintf(fptr, "Query %s:\n",name);

    for (int i = 0; i < n; i++){
      nNearest[i]=NULL;
      knearestDists[i]=-1;
    }

    clock_t begin_true = clock(); // time calculation for the k nearest neighbors with brute force method
    // find with the brute force method the k nearest neighbors for the corresponding query vector
    if(distanceTrueOff!=1)
      listFindNearestNeighbor(inputs,vecTmp,nNearest,knearestDists,dim,-1);

    clock_t end_true = clock();
    double time_spent_true = (double)(end_true - begin_true) / CLOCKS_PER_SEC;

    clock_t begin_lsh = clock(); // time calculation for the k nearest neighbors with LSH
    // find with the help of LSH the k nearest neighbor for the corresponding query vector
    double approximation_factor=-1;
    int found_neighbor=-1;
    nearestNeigborLSH_DiscreteFrechet(lsh,vecTmp,nNearest,knearestDists,fptr,grids,delta,&approximation_factor,&found_neighbor,distanceTrueOff);
    if(approximation_factor>0){
      if(approximation_factor>max_aproximation_factor)
        max_aproximation_factor=approximation_factor;
      if(approximation_factor<min_aproximation_factor)
        min_aproximation_factor=approximation_factor;
    }

    clock_t end_lsh = clock();
    double time_spent_lsh = (double)(end_lsh - begin_lsh) / CLOCKS_PER_SEC;
    if(found_neighbor>0){
      total_lsh_time += time_spent_lsh;
      total_true_time += time_spent_true;
      query_count++;
    }
    // fprintf(fptr, "tApproximateAverage: %f seconds\n",time_spent_lsh);
    // fprintf(fptr, "tTrueAverage: %f seconds\n",time_spent_true);
    fprintf(fptr, "\n");
    fflush(fptr);
    deleteVector(vecTmp);
  }
  fprintf(fptr, "tApproximateAverage: %f seconds\n",total_lsh_time/query_count);
  if(distanceTrueOff!=1){
    fprintf(fptr, "tTrueAverage: %f seconds\n",total_true_time/query_count);
    fprintf(fptr, "Max Approximation Factor: %f\n",max_aproximation_factor);
    fprintf(fptr, "Min Approximation Factor: %f\n",min_aproximation_factor);
  }else{
    fprintf(fptr, "tTrueAverage: True Distances were not computed\n");
    fprintf(fptr, "Min/Max Approximation Factor are not defined without distanceTrue\n");
  }
  fclose(fptr);
  fclose(file);
}

void readQueryFileLSH_ContinuousFrechet(char* queryFile,char* outputFile,LSH lsh,List inputs,double delta,double epsilon,int dim,Grids grid,int distanceTrueOff){

   FILE *file = fopen(queryFile, "r"); // read mode

   if (file == NULL){
      perror("Error while opening the file.\n");
      exit(-1);
   }

  if (feof(file)){ // empty file, return
    return;
  }
  FILE* fptr;
  fptr = fopen(outputFile, "w");
  if(fptr == NULL){
    /* File not created hence exit */
    printf("Unable to create file.\n");
    exit(EXIT_FAILURE);
  }

  int n =1 ; // 1 nearest neighbor
  char buffer[MAX_INPUT_LENGTH];
  Vector nNearest[n]; // here store the true k nearest neighbors
  double knearestDists[n]; // here store the true distances from the k nearest neighbors
  double vec[dim];

  int query_count =0;
  double total_lsh_time = 0.0;
  double total_true_time = 0.0;
  double max_aproximation_factor = -1;
  double min_aproximation_factor = DBL_MAX;
  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }


    char * token = strtok(buffer, "	 ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL, "	");
    int counter = 0;
    while( token != NULL ) {
      vec[counter++]=atof(token);
      token = strtok(NULL, "	");
    }
    Vector vecTmp=initVector(vec,name,dim);

    fprintf(fptr, "Query %s:\n",name);

    for (int i = 0; i < n; i++){
      nNearest[i]=NULL;
      knearestDists[i]=-1;
    }

    clock_t begin_true = clock(); // time calculation for the k nearest neighbors with brute force method
    // find with the brute force method the k nearest neighbors for the corresponding query vector
    if(distanceTrueOff!=1)
      listFindNearestNeighbor(inputs,vecTmp,nNearest,knearestDists,dim,-1);

    clock_t end_true = clock();
    double time_spent_true = (double)(end_true - begin_true) / CLOCKS_PER_SEC;

    clock_t begin_lsh = clock(); // time calculation for the k nearest neighbors with LSH
    double approximation_factor=-1;
    int found_neighbor=-1;
    // find with the help of LSH the k nearest neighbor for the corresponding query vector
    nearestNeigborLSH_ContinuousFrechet(lsh,vecTmp,nNearest,knearestDists,fptr,delta,epsilon,grid,&approximation_factor,&found_neighbor,distanceTrueOff);
    if(approximation_factor>0){
      if(approximation_factor>max_aproximation_factor)
        max_aproximation_factor=approximation_factor;
      if(approximation_factor<min_aproximation_factor)
        min_aproximation_factor=approximation_factor;
    }

    clock_t end_lsh = clock();
    double time_spent_lsh = (double)(end_lsh - begin_lsh) / CLOCKS_PER_SEC;
    if(found_neighbor>0){
      total_lsh_time += time_spent_lsh;
      total_true_time += time_spent_true;
      query_count++;
    }
    // fprintf(fptr, "tApproximateAverage: %f seconds\n",time_spent_lsh);
    // fprintf(fptr, "tTrueAverage: %f seconds\n",time_spent_true);
    fprintf(fptr, "\n");
    fflush(fptr);
    deleteVector(vecTmp);
  }
  if(query_count>0){
    fprintf(fptr, "tApproximateAverage: %f seconds\n",total_lsh_time/query_count);
    if(distanceTrueOff!=1){
      fprintf(fptr, "tTrueAverage: %f seconds\n",total_true_time/query_count);
      fprintf(fptr, "Max Approximation Factor: %f\n",max_aproximation_factor);
      fprintf(fptr, "Min Approximation Factor: %f\n",min_aproximation_factor);
    }else{
      fprintf(fptr, "tTrueAverage: True Distances were not computed\n");
      fprintf(fptr, "Min/Max Approximation Factor are not defined without distanceTrue\n");
    }
  }

  fclose(fptr);
  fclose(file);
}
