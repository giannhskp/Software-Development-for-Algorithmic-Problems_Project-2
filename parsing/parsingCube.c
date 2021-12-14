#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTable.h"
#include "../Hypercube/hypercube.h"
#include "../hashTable/hashTableList/hashTableList.h"

#define OUT    0
#define IN    1
#define FALSE    0
#define TRUE    1

#define MAX_INPUT_LENGTH 10240

// extern int d;

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


int findDimCube(char* fileName){
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
 int dims = countWords(buffer)-1;
 fclose(file);
 return dims;
}


void readFileCube(char* fileName,List *inputs,int *vectorCount,int dim){

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
    int read_result = fscanf(file,"%[^\n]\n",buffer);
    if(read_result<0){ // read a line from the file
      continue;
    }
    double vec[dim];
    char * token = strtok(buffer, "	 ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    name[strlen(name)]='\0';
    token = strtok(NULL, "	");
    int counter = 0;
    while( token != NULL ) {
      vec[counter++]=atof(token);
      token = strtok(NULL, "	");
     }
     Vector vecTmp=initVector(vec,name,dim);
     (*inputs) = listInsert((*inputs),vecTmp,-1);
     (*vectorCount)++;
  }

  fclose(file);
}



void readQueryFileCube(char* queryFile,char* outputFile,HyperCube hc,List inputs,int hammingDist,int m,int dim){

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

  char buffer[MAX_INPUT_LENGTH];
  int n = 1;  // 1 nearest neighbor
  Vector nNearest[n]; // here store the true k nearest neighbors
  double knearestDists[n]; // here store the true distances from the k nearest neighbors
  double vec[dim];

  int query_count =0;
  double total_cube_time = 0.0;
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

    clock_t begin_true = clock(); // time calculation for the brute force method
    // find with the brute force method the k nearest neighbors for the corresponding query vector
    listFindNearestNeighbor(inputs,vecTmp,nNearest,knearestDists,dim,-1);
    clock_t end_true = clock();
    double time_spent_true = (double)(end_true - begin_true) / CLOCKS_PER_SEC;

    clock_t begin_cube = clock(); // time calculation for the Hypercube

    double approximation_factor=-1;
    int found_neighbor = -1;
    // find with the help of Hypercube the k nearest neighbors corresponding query vector
    nearestNeigborHypercube(hc,vecTmp,nNearest,hammingDist,m,knearestDists,fptr,&approximation_factor,&found_neighbor);
    if(approximation_factor>0){
      if(approximation_factor>max_aproximation_factor)
        max_aproximation_factor=approximation_factor;
      if(approximation_factor<min_aproximation_factor)
        min_aproximation_factor=approximation_factor;
    }
    clock_t end_cube = clock();
    double time_spent_cube = (double)(end_cube - begin_cube) / CLOCKS_PER_SEC;
    if(found_neighbor){
      total_cube_time += time_spent_cube;
      total_true_time += time_spent_true;
      query_count++;
    }
    // fprintf(fptr, "tCube: %f seconds\n",time_spent_cube);
    // fprintf(fptr, "tTrue: %f seconds\n",time_spent_true);
    fprintf(fptr, "\n");
    fflush(fptr);
    deleteVector(vecTmp);
  }
  if(query_count>0){
    fprintf(fptr, "tApproximateAverage: %f seconds\n",total_cube_time/query_count);
    fprintf(fptr, "tTrueAverage: %f seconds\n",total_true_time/query_count);
    fprintf(fptr, "Max Approximation Factor: %f\n",max_aproximation_factor);
    fprintf(fptr, "Min Approximation Factor: %f\n",min_aproximation_factor);
  }
  fclose(fptr);
  fclose(file);
}
