#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "Vector/vector.h"
#include "./hashTable/hashTable.h"
#include "LSH/lsh.h"
#include "./parsing/parsingLSH.h"
#include "./hashTable/hashTableList/hashTableList.h"

#define W_DIVIDER 80

int d;
int w;
int k_LSH;
int hashTableSize;
Vector timeVector;

static int wValueCalculation(List list,int numberOfVectorsInFile){
  long double sumDist = 0.0;
  int count=0;
  double persentageToCheck;
  if(numberOfVectorsInFile<=1000){
    persentageToCheck = 0.1;
  }else if(numberOfVectorsInFile<=10000){
    persentageToCheck = 0.001;
  }else if (numberOfVectorsInFile<=100000){
    persentageToCheck = 0.0001;
  }else{
    persentageToCheck = 0.000001;
  }
  persentageToCheck = 0.0001;
  int stopBound = persentageToCheck*numberOfVectorsInFile*numberOfVectorsInFile;
  while(list!=NULL){
    List nested = list;
    while(nested!=NULL){
      if(count>stopBound){
        return floor(sumDist/count);
      }
      sumDist +=distance_metric(getVector(list),getVector(nested),d);
      count++;
      nested = getNext(nested);
    }
    list=getNext(list);
  }
  return floor(sumDist/count);
}


void vectorTimeSeriesLSH(char* arg_inputFile,char* arg_queryFile,int arg_k_LSH,int arg_L,char* arg_outputFile)  {

  char str[200];
  char inputFile[100];
  strcpy(inputFile,arg_inputFile);
  char queryFile[100];
  strcpy(queryFile,arg_queryFile);
  char outputFile[100];
  strcpy(outputFile,arg_outputFile);
  int l=arg_L;
  k_LSH = arg_k_LSH;
  hashTableSize = 1000;



  srand(time(NULL));

  LSH lsh;
  List list;
  clock_t begin = clock();
  d = findDimLSH(inputFile);
  printf("DIMENSION = %d\n",d);
  list = initializeList();
  int numberOfVectorsInFile = 0;
  readFileLSH(inputFile,&list,&numberOfVectorsInFile);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Parsed input file in : %f seconds\n",time_spent);
  printf("Number of vectors in input file: %d\n",numberOfVectorsInFile);
  hashTableSize=numberOfVectorsInFile/16;

  printf("Finding optimal value of w based on the input file\n");
  begin = clock();
  w = wValueCalculation(list,numberOfVectorsInFile);
  w /= W_DIVIDER;
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Found value of w in %f seconds, w = %d\n",time_spent,w );

  begin = clock();
  lsh = initializeLSH(l);
  insertFromListToLSH(list,lsh);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Created LSH in : %f seconds\n",time_spent);

  readQueryFileLSH(queryFile,outputFile,lsh,list);

  destroyLSH(lsh);
  listDelete(list,0);
}


void vectorTimeSeriesLSHFrechetDiscrete(char* arg_inputFile,char* arg_queryFile,int arg_k_LSH,int arg_L,char* arg_outputFile,double arg_delta){
  char str[200];
  char inputFile[100];
  strcpy(inputFile,arg_inputFile);
  char queryFile[100];
  strcpy(queryFile,arg_queryFile);
  char outputFile[100];
  strcpy(outputFile,arg_outputFile);
  int l=arg_L;
  double delta=arg_delta;
  k_LSH = arg_k_LSH;
  hashTableSize = 1000;



  srand(time(NULL));

  double sum=0.0;
  double time[d];
  for(int i=0;i<d;i++){
    time[i]=sum;
    sum+=1.0;
  }
  timeVector=initVector(time,"time");


  LSH lsh;
  List list;
  clock_t begin = clock();
  d = findDimLSH(inputFile);
  printf("DIMENSION = %d\n",d);
  list = initializeList();
  int numberOfVectorsInFile = 0;
  readFileLSH(inputFile,&list,&numberOfVectorsInFile);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Parsed input file in : %f seconds\n",time_spent);
  printf("Number of vectors in input file: %d\n",numberOfVectorsInFile);
  hashTableSize=numberOfVectorsInFile/16;

  printf("Finding optimal value of w based on the input file\n");
  begin = clock();
  w = wValueCalculation(list,numberOfVectorsInFile);
  // w /= W_DIVIDER;
  w = w/10;
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Found value of w in %f seconds, w = %d\n",time_spent,w );

  begin = clock();
  lsh = initializeLSH(l);
  Grids grids = initializeGrids(delta,l);
  insertTimeSeriesFromListToLSH(list,lsh,grids,timeVector,delta);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Created LSH in : %f seconds\n",time_spent);

  readQueryFileLSH_DiscreteFrechet(queryFile,outputFile,lsh,list,grids,timeVector,delta);

  destroyLSH(lsh);
  listDelete(list,0);
}

// void vectorTimeSeriesLSHFrechetContinuous(char* arg_inputFile,char* arg_queryFile,int arg_k_LSH,int arg_L,char* arg_outputFile,double arg_delta){
//
//   // TO DO
// }
