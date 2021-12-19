#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "Vector/vector.h"
#include "./hashTable/hashTable.h"
#include "./hashTable/hashTableList/hashTableList.h"
#include "Hypercube/hypercube.h"
#include "./parsing/parsingCluster.h"
#include "./Clustering/clustering.h"
#include "./BinaryTree/binaryTree.h"

#define W_VALUE 4

int new_dimension;
int m;
int probes;
int w;
int numOfVecs;
int hashTableSize;
int silhouette;
int k_LSH;
int complete;

void printOptions(){
  printf("_________________Options____________________\n\n");
		printf("1. /repeat <new_input_file> <configuration file> <method: Classic OR LSH or Hypercube> <output file> \n\n");
	printf("2. /exit\n\n");
	printf("_____________________________________\n\n");
}



int main(int argc, char *argv[]) {
  srand(time(NULL));
  silhouette=0;
  complete=0;
  char str[200];
  char inputFile[200];
  int inputflag=0;
  char confFile[200];
  int confflag=0;
  char outputFile[200];
  int outputflag=0;
  char update[200];
  strcpy(update,"Mean Vector"); // default
  char assignment[200];
  strcpy(assignment,"Classic"); // default

  double delta = 1;


  for(int i = 1 ; i < argc ; i++){
    strcpy(str,argv[i]);
    if(strcmp(str,"-i")==0 && (argc > i+1)){
      inputflag++;
      strcpy(inputFile,argv[i+1]);
      printf("Given input File : %s\n", inputFile);
    }
    else if(strcmp(str,"-c")==0 && (argc > i+1)){
      confflag++;
      strcpy(confFile,argv[i+1]);
      printf("Given configuration File : %s\n", confFile);
    }
    else if(strcmp(str,"-silhouette")==0){
      printf("silhouette option ON.\n");
      silhouette=1;
    }
    else if(strcmp(str,"-complete")==0){
      printf("complete option ON.\n");
      complete=1;
    }
    else if(strcmp(str,"-o")==0 && (argc > i+1)){
      outputflag++;
      strcpy(outputFile,argv[i+1]);
      printf("Given output File : %s\n", outputFile);
    }
    else if(strcmp(str,"-update")==0 && (argc > i+2)){
      strcpy(update,argv[i+1]); // default
      strcat(update," ");
      strcat(update,argv[i+2]);
      printf("Given Update Method : %s\n", update);
    }
    else if(strcmp(str,"-assignment")==0 && (argc > i+1)){
      strcpy(assignment,argv[i+1]); // default
      printf("Given Assignment Method : %s\n", assignment);
    }
    else if(strcmp(str,"-delta")==0 && (argc > i+1)){
      delta = atof(argv[i+1]);
      printf("Delta Given : %f\n", delta);
    }

  }

  if(!inputflag){
    printf(">Input file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    sscanf(str,"%s\n",inputFile);
    printf("Given input File : %s\n", inputFile);
  }
  if(!confflag){
    printf(">Conf file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    sscanf(str,"%s\n",confFile);
    printf("Given configuration File : %s\n", confFile);
  }
  if(!outputflag){
    printf(">Output file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    sscanf(str,"%s\n",outputFile);
    printf("Given output File : %s\n", outputFile);
  }

  FILE* fptr;
  List list;
  int repeat=1;
  while(1){
    if(repeat){
      clock_t begin = clock();
      int dim = findDim(inputFile);
      printf("DIMENSION = %d\n",dim);
      int numOfClusters=-100,l=3,mHyper=10,probes=2;
      new_dimension=3;
      k_LSH=4;
      w=W_VALUE;
      readConfFile(confFile,&numOfClusters,&l,&mHyper,&probes);
      while(numOfClusters<2){ // clusters number should be >=2
        if(numOfClusters!=-100)
          printf("\n(Clusters number should be >=2)\n");
        printf("\n>Give number of clusters: ");
        fflush(stdin); // clear stdin buffer
        if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
          perror("Error reading string with fgets\n");
          exit(1);
        }
        numOfClusters=atoi(str);
        printf("numOfClusters : %d\n\n", numOfClusters);
      }
      list = initializeList();
      readFile(inputFile,&list,&numOfVecs,dim);
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      printf("Parsed input file in : %f seconds\n",time_spent);
      printf("Number of vectors in input file: %d\n",numOfVecs);
      fflush(stdout);
      fptr = fopen(outputFile, "w");
      if(fptr == NULL){
        /* File not created hence exit */
        printf("Unable to create file.\n");
        exit(EXIT_FAILURE);
      }

      clustering(list,fptr,assignment,update,numOfClusters,l,mHyper,probes,dim,delta);
    }
    repeat=0;


    printOptions(); // just printing the commands options for the user


    char command[200];

      printf("\n");
      printf(">Enter a command: ");
      fflush(stdin); // clear stdin buffer
      if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
        perror("Error reading string with fgets\n");
        exit(1);
      }
      else if(strstr(str, "/repeat") != NULL) {
        repeat=1;
        if(strcmp(assignment,"Classic")==0)
          listDelete(list,1);
        else
          listDelete(list,0);
        sscanf(str,"%s %s %s %s %s %s\n",command,inputFile,confFile,update,assignment,outputFile);
        printf("FILE: %s\n",inputFile);
        printf("Given configuration File : %s\n", confFile);
        printf("Given update method : %s\n", update);
        printf("Given assignment method : %s\n", assignment);
        printf("Given output File : %s\n", outputFile);
        fclose(fptr);
        continue;
      }
      else if(strcmp(str,"/exit\n")==0){
        break;
      }
      else{
        printf("\n\n  --- Wrong command ! Please, try again. ---  \n\n");
        printOptions(); // just printing the commands options for the user
        continue;
      }

  }

  if(strcmp(assignment,"Classic")==0)
    listDelete(list,1);
  else
    listDelete(list,0);

  fclose(fptr);
  return 0;
}
