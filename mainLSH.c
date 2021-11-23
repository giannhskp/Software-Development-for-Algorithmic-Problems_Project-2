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

int wValueCalculation(List list,int numberOfVectorsInFile){
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
  int stopBound = persentageToCheck*numberOfVectorsInFile*numberOfVectorsInFile;
  while(list!=NULL){
    List nested = list;
    while(nested!=NULL){
      if(count>stopBound){
        return floor(sumDist/count);
      }
      sumDist += distance_metric(getVector(list),getVector(nested),d);
      count++;
      nested = getNext(nested);
    }
    list=getNext(list);
  }
  return floor(sumDist/count);
}

void printOptions(){
  printf("_________________Options____________________\n\n");
	printf("1. /repeat <new_query_file> <output file>\n\n");
	printf("2. /exit\n\n");
	printf("_____________________________________\n\n");
}

int main(int argc, char *argv[])  {

  int option;
  char str[200];
  char inputFile[100];
  int inputflag=0;
  char queryFile[100];
  int queryflag=0;
  char outputFile[100];
  int outputflag=0;
  int l=5;
  int n=1;
  double radius=10000;
  hashTableSize = 1000;
  k_LSH = 4;

  while((option = getopt(argc, argv, "i:q:k:L:o:N:R:")) != -1){
     switch(option){
        case 'i':
        inputflag++;
        strcpy(inputFile,optarg);
        printf("Given input File : %s\n", inputFile);
        break;

        case 'q':
        queryflag++;
        strcpy(queryFile,optarg);
        printf("Given query File : %s\n", queryFile);
        break;

        case 'k':
        k_LSH=atoi(optarg);
        printf("k : %d\n", k_LSH);
        break;

        case 'L':
        l=atoi(optarg);
        printf("L : %d\n", l);
        break;

        case 'o':
        outputflag++;
        strcpy(outputFile,optarg);
        printf("Given output File : %s\n", outputFile);
        break;

        case 'N':
        n=atoi(optarg);
        printf("number of nearest neighbors: %d\n", n);
        break;

        case 'R':
         radius=atof(optarg);
         printf("Radius : %f\n", radius);
         break;

        case ':':
         printf("option needs a value\n");
         break;

        default: /* '?' */
          fprintf(stderr, "Usage: %s –i <input file> –q <query file> –k <int> -L <int> -ο <output file> -Ν <number of nearest> -R <radius>\n",argv[0]);
          exit(EXIT_FAILURE);
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
  char command[200];

  LSH lsh;
  List list;
  int repeat=1;
  clock_t begin = clock();
  d = findDim(inputFile);
  printf("DIMENSION = %d\n",d);
  list = initializeList();
  int numberOfVectorsInFile = 0;
  readFile(inputFile,&list,&numberOfVectorsInFile);
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

  while(1){
    if(repeat){
      readQueryFile(queryFile,outputFile,lsh,list,n,radius);
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
      sscanf(str,"%s %s %s\n",command,queryFile,outputFile);
      printf("query File: %s\n",queryFile);
      printf("output File: %s\n",outputFile);
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


  destroyLSH(lsh);
  listDelete(list,0);

  return 0;
}
