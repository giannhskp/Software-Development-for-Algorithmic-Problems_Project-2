#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "mainLSH.h"
#include "mainCube.h"

#define FILTERING_E 0.1


void printOptions(){
  printf("_________________Options____________________\n\n");
	printf("1. /repeat <new_query_file> <output file> <algorithm> <metric>\n\n");
	printf("2. /exit\n\n");
	printf("_____________________________________\n\n");
}

char *distanceMetric;


int main(int argc, char *argv[])  {
  char str[200];
  char inputFile[100];
  int inputflag=0;
  char queryFile[100];
  int queryflag=0;
  char outputFile[100];
  char algorithm[100];
  char metric[100];
  int m=10;
  double delta=1;
  int outputflag=0;
  int l=5;
  int k_LSH = 4;
  int new_dimension = 14;
  int probes=2;
  int distanceTrueOff = 0;

  for(int i = 1 ; i < argc ; i++){
    strcpy(str,argv[i]);
    if(strcmp(str,"-i")==0 && (argc > i+1)){
      inputflag++;
      strcpy(inputFile,argv[i+1]);
      printf("Given input File : %s\n", inputFile);
    }
    else if(strcmp(str,"-q")==0 && (argc > i+1)){
      queryflag++;
      strcpy(queryFile,argv[i+1]);
      printf("Given query File : %s\n", queryFile);
    }
    else if(strcmp(str,"-k")==0 && (argc > i+1)){
      k_LSH=atoi(argv[i+1]);
      new_dimension = k_LSH;
      printf("k : %d\n", k_LSH);
    }
    else if(strcmp(str,"-L")==0 && (argc > i+1)){
      l=atoi(argv[i+1]);
      printf("L : %d\n", l);
    }
    else if(strcmp(str,"-o")==0 && (argc > i+1)){
      outputflag++;
      strcpy(outputFile,argv[i+1]);
      printf("Given output File : %s\n", outputFile);
    }
    else if(strcmp(str,"-algorithm")==0 && (argc > i+1)){
      strcpy(algorithm,argv[i+1]);
      printf("Given algorithm : %s\n", algorithm);
    }
    else if(strcmp(str,"-metric")==0 && (argc > i+1)){
      strcpy(metric,argv[i+1]);
      printf("Given metric: %s\n", metric);
    }
    else if(strcmp(str,"-delta")==0 && (argc > i+1)){
      delta=atof(argv[i+1]);
      printf("Delta: %f\n", delta);
    }
    else if(strcmp(str,"-probes")==0 && (argc > i+1)){
      probes=atoi(argv[i+1]);
      printf("probes : %d\n", probes);
    }
    else if(strcmp(str,"-M")==0 && (argc > i+1)){
      m=atoi(argv[i+1]);
      printf("M : %d\n", m);
    }
    else if(strcmp(str,"-distanceTrueOff")==0){
      distanceTrueOff=1;
      printf("distanceTrueOff : Enabled\n");
    }
  }

  if(!inputflag){
    printf(">Input file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    strcpy(inputFile,str);
    printf("Given input File : %s\n", inputFile);
  }
  if(!queryflag){
    printf(">Query file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    strcpy(queryFile,str);
    printf("Given query File : %s\n", queryFile);
  }
  if(!outputflag){
    printf(">Output file name: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    strcpy(outputFile,str);
    printf("Given output File : %s\n", outputFile);
  }

  srand(time(NULL));
  char command[100];
  int repeat=1;
  while(1){
    if(repeat){
      if(strcmp(algorithm,"LSH")==0){ // if LSH was given as algorithm
        distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));      // set distance metric to l2
        strcpy(distanceMetric,"l2");
        fprintf(stdout,"Algorithm: LSH\n");
        vectorTimeSeriesLSH(inputFile,queryFile,k_LSH,l,outputFile,distanceTrueOff);  // call the correspoding function of LSH/lsh.c file
      }
      else if(strcmp(algorithm,"Hypercube")==0){ // if Hypercube was given as algorithm
        distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));     // set distance metric to l2
        strcpy(distanceMetric,"l2");
        fprintf(stdout,"Algorithm: LSH\n");
        vectorTimeSeriesHypecube(inputFile,queryFile,new_dimension,m,probes,outputFile,distanceTrueOff);    // call the correspoding function of LSH/lsh.c file
      }
      else if(strcmp(algorithm,"Frechet")==0){ // if Frechet was given as algorithm
        if(strcmp(metric,"discrete")==0){ // if Frechet was given as algorithm and Discrete as metric
          distanceMetric=malloc(sizeof(char)*(strlen("discreteFrechet")+1));    // set distance metric to discrete frechet
          strcpy(distanceMetric,"discreteFrechet");
          fprintf(stdout,"Algorithm: Frechet Discrece\n");
          vectorTimeSeriesLSHFrechetDiscrete(inputFile,queryFile,k_LSH,l,outputFile,delta,distanceTrueOff);   // call the correspoding function of LSH/lsh.c file
        }
        else if(strcmp(metric,"continuous")==0){   // if Frechet was given as algorithm and Continuous as metric
          distanceMetric=malloc(sizeof(char)*(strlen("continuousFrechet")+1));  // set distance metric to continuous frechet
          strcpy(distanceMetric,"continuousFrechet");
          fprintf(stdout,"Algorithm: Frechet Continuous\n");
          vectorTimeSeriesLSHFrechetContinuous(inputFile,queryFile,k_LSH,outputFile,delta,FILTERING_E,distanceTrueOff);   // call the correspoding function of LSH/lsh.c file
        }
        else{
          printf("WRONG METRIC, %s\n",metric);
        }
      }
      else{
        printf("%s\n",algorithm);
        printf("INVALID Algorithm NAME!\n");
        exit(EXIT_FAILURE);
      }
      free(distanceMetric);
    }
    repeat=0;

    printOptions(); // just printing the commands options for the user



    printf("\n");
    printf(">Enter a command: ");
    fflush(stdin); // clear stdin buffer
    if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
      perror("Error reading string with fgets\n");
      exit(1);
    }
    else if(strstr(str, "/repeat") != NULL) {
      repeat=1;
      sscanf(str,"%s %s %s %s %s\n",command,queryFile,outputFile,algorithm,metric);
      printf("query File: %s\n",queryFile);
      printf("output File: %s\n",outputFile);
      printf("algorithm: %s\n",algorithm);
      printf("metric: %s\n",metric);
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


  return 0;
}
