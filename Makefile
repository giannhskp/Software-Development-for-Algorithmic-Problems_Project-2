HASHTABLELIST=./hashTable/hashTableList
HASHTABLE=./hashTable
PARSING=./parsing
VECTOR=./Vector
LSH=./LSH
HASHMAP = ./Hypercube/HashMap
HYPERCUBE = ./Hypercube
CLUSTER = ./Clustering

CC=gcc
CFLAGS= -g -Wall -I$(HASHTABLELIST) -I$(HASHTABLE) -I$(PARSING) -I$(VECTOR) -I$(LSH) -I$(CLUSTER)

OBJ1= mainLSH.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingLSH.o $(VECTOR)/vector.o  $(LSH)/lsh.o $(LSH)/helperFunctions.o
OBJ2= mainCube.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingCube.o $(VECTOR)/vector.o $(HYPERCUBE)/hypercube.o $(HASHMAP)/hashmap.o $(LSH)/helperFunctions.o
OBJ3= mainCluster.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingCluster.o $(VECTOR)/vector.o $(HASHMAP)/hashmap.o $(LSH)/lsh.o $(LSH)/helperFunctions.o $(CLUSTER)/clustering.o $(CLUSTER)/clusterHelpingFuns.o $(CLUSTER)/kmeansPlusPlus.o $(HYPERCUBE)/hypercube.o

EXEC1 = demolsh
EXEC2 = cube
EXEC3 = cluster

all: $(EXEC1) $(EXEC2) $(EXEC3)

$(EXEC1): $(OBJ1)
	$(CC) $(CFLAGS) $(OBJ1) -o $(EXEC1) -lm

$(EXEC2): $(OBJ2)
	$(CC) $(CFLAGS) $(OBJ2) -o $(EXEC2) -lm

$(EXEC3): $(OBJ3)
	$(CC) $(CFLAGS) $(OBJ3) -o $(EXEC3) -lm


.PHONY: clean

clean:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(EXEC1) $(EXEC2) $(EXEC3)
