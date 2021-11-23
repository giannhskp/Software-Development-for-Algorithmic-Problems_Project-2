#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../Hypercube/hypercube.h"

#define OUT    0
#define IN    1
#define FALSE    0
#define TRUE    1

#define MAX_INPUT_LENGTH 1024

extern int d;
extern int k_LSH;
extern int new_dimension;

// returns the number of words in string
int countWords(char *str){
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
 int dims = countWords(buffer)-1;
 fclose(file);
 return dims;
}


void readConfFile(char* fileName,int * numOfClusters,int *l,int *mHyper,int *probes){
  FILE *file = fopen(fileName, "r"); // read mode

  if (file == NULL){
     perror("Error while opening the file.\n");
     exit(-1);
  }

 if (feof(file)){ // empty file, return
   return;
 }

  char buffer[MAX_INPUT_LENGTH];
  char command[100];
  char temp[100];
  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }
    if(strstr(buffer, "number_of_clusters:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      *numOfClusters=atoi(temp);
      printf("number_of_clusters: %d\n",*numOfClusters);
      continue;
    }
    else if(strstr(buffer, "number_of_vector_hash_tables:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      *l=atoi(temp);
      printf("number_of_vector_hash_tables: %d\n",*l);
      continue;
    }
    else if(strstr(buffer, "number_of_vector_hash_functions:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      k_LSH=atoi(temp);
      printf("number_of_vector_hash_functions: %d\n",k_LSH);
      continue;
    }
    else if(strstr(buffer, "max_number_M_hypercube:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      *mHyper=atoi(temp);
      printf("max_number_M_hypercube: %d\n",*mHyper);
      continue;
    }
    else if(strstr(buffer, "number_of_hypercube_dimensions:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      new_dimension=atoi(temp);
      printf("number_of_hypercube_dimensions: %d\n",new_dimension);
      continue;
    }
    else if(strstr(buffer, "number_of_probes:") != NULL) {
      sscanf(buffer,"%s %s\n",command,temp);
      *probes=atoi(temp);
      printf("number_of_probes: %d\n",*probes);
      continue;
    }
  }
  fclose(file);
}


void readFile(char* fileName,List *list,int *numOfVecs){

   FILE *file = fopen(fileName, "r"); // read mode

   if (file == NULL){
      perror("Error while opening the file.\n");
      exit(-1);
   }

  if (feof(file)){ // empty file, return
    return;
  }
  (*numOfVecs)=0;


  char buffer[MAX_INPUT_LENGTH];


  while(!feof(file)){
    fflush(stdin);  // clear stdin buffer
    if(fscanf(file,"%[^\n]\n",buffer)<0){ // read a line from the file
      continue;
    }

    double vec[d];
    char * token = strtok(buffer," ");
    char name[MAX_INPUT_LENGTH];
    strcpy(name,token);
    token = strtok(NULL," ");
     int counter = 0;
     while( token != NULL ) {
        vec[counter++]=atof(token);
        token = strtok(NULL," ");
     }
     Vector vecTmp=initVector(vec,name);
     (*list) = listInsert((*list),vecTmp,-1);
     (*numOfVecs)++;

  }


  fclose(file);


}
