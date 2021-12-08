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
	printf("1. /repeat <new_query_file> <output file>\n\n");
	printf("2. /exit\n\n");
	printf("_____________________________________\n\n");
}

char *distanceMetric;


int main(int argc, char *argv[])  {

  int option;
  char str[200];
  char inputFile[100];
  int inputflag=0;
  char queryFile[100];
  int queryflag=0;
  char outputFile[100];
  char algorithm[100];
  char metric[100];
  int m=10;
  double delta;
  int outputflag=0;
  int l=5;
  int k_LSH = 4;
  int new_dimension = 14;
  int probes=2;

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

  }

  // while((option = getopt(argc, argv, "i:q:k:L:o:M:m:a:d:p:")) != -1){
  //    switch(option){
  //       case 'i':
  //       inputflag++;
  //       strcpy(inputFile,optarg);
  //       printf("Given input File : %s\n", inputFile);
  //       break;
  //
  //       case 'q':
  //       queryflag++;
  //       strcpy(queryFile,optarg);
  //       printf("Given query File : %s\n", queryFile);
  //       break;
  //
  //       case 'k':
  //       k_LSH=atoi(optarg);
  //       new_dimension = k_LSH;
  //       printf("k : %d\n", k_LSH);
  //       break;
  //
  //       case 'L':
  //       l=atoi(optarg);
  //       printf("L : %d\n", l);
  //       break;
  //
  //       case 'o':
  //       outputflag++;
  //       strcpy(outputFile,optarg);
  //       printf("Given output File : %s\n", outputFile);
  //       break;
  //
  //       case 'a':
  //       if(strcmp(argv[optind-1],"-algorithm")==0){
  //         strcpy(algorithm,argv[optind]);
  //         printf("Given algorithm : %s\n", algorithm);
  //       }
  //       break;
  //
  //       case 'm':
  //       if(strcmp(argv[optind-1],"-metric")==0){
  //         strcpy(metric,argv[optind]);
  //         printf("Given metric: %s\n", metric);
  //       }
  //       break;
  //
  //       case 'd':
  //       if(strcmp(argv[optind-1],"-delta")==0){
  //         delta=atof(argv[optind]);
  //         printf("Delta: %f\n", delta);
  //       }
  //       break;
  //
  //       case 'p':
  //       if(strcmp(argv[optind-1],"-probes")==0){
  //         probes=atoi(argv[optind]);
  //         printf("probes : %d\n", probes);
  //       }
  //       break;
  //
  //       case 'M':
  //       m=atoi(optarg);
  //       printf("M : %d\n", m);
  //       break;
  //
  //       case ':':
  //        printf("option needs a value\n");
  //        break;
  //
  //       default:
  //         fprintf(stderr, "Usage: %s –i <input file> –q <query file> –k <int> -L <int> -M <int> -probes <int> -ο <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discreteor continuous | only for –algorithm Frechet> -delta <double>\n",argv[0]);
  //         exit(EXIT_FAILURE);
  //    }
  // }

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
  char command[200];


  if(strcmp(algorithm,"LSH")==0){
    distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));
    strcpy(distanceMetric,"l2");
    fprintf(stdout,"Algorithm: LSH\n");
    vectorTimeSeriesLSH(inputFile,queryFile,k_LSH,l,outputFile);
  }
  else if(strcmp(algorithm,"Hypercube")==0){
    distanceMetric=malloc(sizeof(char)*(strlen("l2")+1));
    strcpy(distanceMetric,"l2");
    fprintf(stdout,"Algorithm: LSH\n");
    vectorTimeSeriesHypecube(inputFile,queryFile,new_dimension,m,probes,outputFile);
  }
  else if(strcmp(algorithm,"Frechet")==0){
    if(strcmp(metric,"discrete")==0){
      distanceMetric=malloc(sizeof(char)*(strlen("discreteFrechet")+1));
      strcpy(distanceMetric,"discreteFrechet");
      fprintf(stdout,"Algorithm: Frechet Discrece\n");
      vectorTimeSeriesLSHFrechetDiscrete(inputFile,queryFile,k_LSH,l,outputFile,delta);
    }
    else if(strcmp(metric,"continuous")==0){
      distanceMetric=malloc(sizeof(char)*(strlen("continuousFrechet")+1));
      strcpy(distanceMetric,"continuousFrechet");
      fprintf(stdout,"Algorithm: Frechet Continuous\n");
      vectorTimeSeriesLSHFrechetContinuous(inputFile,queryFile,k_LSH,outputFile,delta,FILTERING_E);
    }
    else{
      // TODO: GIVE METRIC
      printf("WRONG METRIC, %s\n",metric);
    }
  }
  else{
    printf("%s\n",algorithm);
    printf("INVALID Algorithm NAME!\n");
    exit(EXIT_FAILURE);
  }


  // int repeat=1;
  // while(1){
  //   if(repeat){
  //     // readQueryFile(queryFile,outputFile,lsh,list,n,radius);
  //
  //
  //   }
  //   repeat=0;
  //
  //   printOptions(); // just printing the commands options for the user
  //
  //
  //
  //   printf("\n");
  //   printf(">Enter a command: ");
  //   fflush(stdin); // clear stdin buffer
  //   if (fgets(str, sizeof(char)*200, stdin) == NULL) { // read a command
  //     perror("Error reading string with fgets\n");
  //     exit(1);
  //   }
  //   else if(strstr(str, "/repeat") != NULL) {
  //     repeat=1;
  //     sscanf(str,"%s %s %s\n",command,queryFile,outputFile);
  //     printf("query File: %s\n",queryFile);
  //     printf("output File: %s\n",outputFile);
  //     continue;
  //   }
  //   else if(strcmp(str,"/exit\n")==0){
  //     break;
  //   }
  //   else{
  //     printf("\n\n  --- Wrong command ! Please, try again. ---  \n\n");
  //     printOptions(); // just printing the commands options for the user
  //     continue;
  //   }
  //
  // }

  free(distanceMetric);
  return 0;
}
