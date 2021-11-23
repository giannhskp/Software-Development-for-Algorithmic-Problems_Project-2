#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTable.h"
#include "../LSH/lsh.h"
#include "../hashTable/hashTableList/hashTableList.h"

#define OUT    0
#define IN    1
#define FALSE    0
#define TRUE    1

#define MAX_INPUT_LENGTH 1024

extern int d;

// returns number of words in str
int countWords(char *str){
    int state = OUT;
    int wc = 0;  // word count

    // Scan all characters one by one
    while (*str)
    {
        // If next character is a separator, set the
        // state as OUT
        if (*str == ' ' || *str == '\n' || *str == '\t' || *str == '\0')
            state = OUT;

        // If next character is not a word separator and
        // state is OUT, then set the state as IN and
        // increment word count
        else if (state == OUT)
        {
            state = IN;
            ++wc;
        }

        // Move to next character
        ++str;
    }

    return wc;
}

// returns number of words in str
int countWords2(char *str){
    char * token = strtok(str, " ");
    token = strtok(NULL, " ");
     // loop through the string to extract all other tokens
     int counter = 0;
     while( token != NULL ) {
        counter++;
        token = strtok(NULL, " ");
     }

    return counter;
}

int countLines(FILE *fp){
  // count the lines of the given file
  int count=0;
  if(fp==NULL){ // error checking
    perror("Error");
    exit(-1);
  }
  while(!feof(fp)){
    if(fgetc(fp)=='\n'){
      count++;
    }
  }
  return count;
}

int findDim(char* fileName){
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
 int dims = countWords2(buffer);
 fclose(file);
 return dims-1;
}


void readFile(char* fileName,List *inputs,int *vectorCount){

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

    double vec[d];
    char * token = strtok(buffer, " ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL, " ");
     // loop through the string to extract all other tokens
     int counter = 0;
     while( token != NULL ) {
        vec[counter++]=atof(token);
        token = strtok(NULL, " ");
     }
     Vector vecTmp=initVector(vec,name);
     (*inputs) = listInsert((*inputs),vecTmp,-1);
     (*vectorCount)++;


  }

  fclose(file);
}


void readQueryFile(char* queryFile,char* outputFile,LSH lsh,List inputs,int n,double radius){

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
  Vector nNearest[n]; // here store the true k nearest neighbors
  double knearestDists[n]; // here store the true distances from the k nearest neighbors
  double vec[d];


  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }

    for (int i = 0; i < n; i++){
      nNearest[i]=NULL;
      knearestDists[i]=-1;
    }

    int id;
    char * token = strtok(buffer, " ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    id=atoi(token);
    token = strtok(NULL, " ");
    int counter = 0;
    while( token != NULL ) {
      vec[counter++]=atof(token);
      token = strtok(NULL, " ");
    }
    Vector vecTmp=initVector(vec,name);

    fprintf(fptr, "Query %d:\n",id);

    clock_t begin_true = clock(); // time calculation for the k nearest neighbors with brute force method
    // find with the brute force method the k nearest neighbors for the corresponding query vector
    if(n==1)
      listFindNearestNeighbor(inputs,vecTmp,nNearest,knearestDists,d,-1);
    else
      listFindKNearestNeighbors(inputs,vecTmp,nNearest,knearestDists,d,n,-1);

    clock_t end_true = clock();
    double time_spent_true = (double)(end_true - begin_true) / CLOCKS_PER_SEC;

    clock_t begin_lsh = clock(); // time calculation for the k nearest neighbors with LSH
    // find with the help of LSH the k nearest neighbor for the corresponding query vector
    if(n==1)
      nearestNeigborLSH(lsh,vecTmp,knearestDists,fptr);
    else
      kNearestNeighborsLSH(lsh,vecTmp,n,knearestDists,fptr);

    clock_t end_lsh = clock();
    double time_spent_lsh = (double)(end_lsh - begin_lsh) / CLOCKS_PER_SEC;
    fprintf(fptr, "tLSH: %f seconds\n",time_spent_lsh);
    fprintf(fptr, "tTrue: %f seconds\n",time_spent_true);


    clock_t begin_radius = clock(); // time calculation for the range search with LSH
    // neighbors range search of LSH for the corresponding query vector
    radiusNeigborsLSH(lsh,vecTmp,radius,fptr);
    clock_t end_radius = clock();
    double time_spent_radius = (double)(end_radius - begin_radius) / CLOCKS_PER_SEC;
    fprintf(fptr, "tRadiusSearch: %f seconds\n\n\n",time_spent_radius);
    deleteVector(vecTmp);
  }
  fclose(fptr);
  fclose(file);
}
