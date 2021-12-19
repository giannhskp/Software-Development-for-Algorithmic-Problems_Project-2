#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "Vector/vector.h"
#include "./hashTable/hashTable.h"
#include "Hypercube/hypercube.h"
#include "./parsing/parsingCube.h"
#include "./hashTable/hashTableList/hashTableList.h"


int new_dimension;
int m;
int probes;
int w;


void vectorTimeSeriesHypecube(char* arg_inputFile,char* arg_queryFile,int arg_new_dimension,int arg_M,int arg_probes,char* arg_outputFile,int distanceTrueOff) {
  srand(time(NULL));
  char str[200];
  char inputFile[100];
  strcpy(inputFile,arg_inputFile);
  char queryFile[100];
  strcpy(queryFile,arg_queryFile);
  char outputFile[100];
  strcpy(outputFile,arg_outputFile);

  new_dimension=arg_new_dimension;
  int m=arg_M;
  int probes=arg_probes;


  if(new_dimension>25){
    printf("You chose a big value for k (d') and the system may not be capable of allocating all the memory needed by the hypecube\n");
    printf("Do you wish to continue with k=%d ? (Press \"y\" to continue or \"n\" to change the value of k)\n",new_dimension);
    while(1){
      fflush(stdin); // clear stdin buffer
      if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
        perror("Error reading string with fgets\n");
        exit(1);
      }
      char ans[100];
      sscanf(str,"%s\n",ans);
      if(strcmp(ans,"y")==0){
        break;
      }else if(strcmp(ans,"n")==0){
        printf("Give new value for k: ");
        fflush(stdin); // clear stdin buffer
        if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
          perror("Error reading string with fgets\n");
          exit(1);
        }
        sscanf(str,"%s\n",ans);
        new_dimension = atoi(ans);
        break;
      }else{
        printf("Wrong input! Press \"y\" to continue or \"n\" to change the value of k\n");
      }
    }
  }


  List list;
  clock_t begin = clock();

  int dim = findDimCube(inputFile);
  printf("DIMENSION = %d\n",dim);

  list = initializeList();
  int numberOfVectorsInFile = 0;
  readFileCube(inputFile,&list,&numberOfVectorsInFile,dim);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Parsed input file in : %f seconds\n",time_spent);
  printf("Number of vectors in input file: %d\n",numberOfVectorsInFile);

  printf("Finding optimal value of w based on the input file\n");

  w=600;

  printf("Found value of w = %d\n",w );

  HyperCube hc;
  begin = clock();
  hc = initializeHyperCube(dim);
  insertFromListToHyperCube(list,hc);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Created HyperCube in : %f seconds\n",time_spent);

  readQueryFileCube(queryFile,outputFile,hc,list,probes,m,dim,distanceTrueOff);

  deleteHyperCube(hc);
  listDelete(list,0);
}
