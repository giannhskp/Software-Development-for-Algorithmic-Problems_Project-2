HASHTABLELIST=./hashTable/hashTableList
HASHTABLE=./hashTable
PARSING=./parsing
VECTOR=./Vector
LSH=./LSH
HASHMAP = ./Hypercube/HashMap
HYPERCUBE = ./Hypercube
CLUSTER = ./Clustering
FRECHET = ./FrechetDistance
BINARYTREE = ./BinaryTree
FRECHERCONT=./Fred-master/src

CC=gcc
CCPLUSPLUS=g++
CFLAGS= -g -Wall -I$(HASHTABLELIST) -I$(HASHTABLE) -I$(PARSING) -I$(VECTOR) -I$(LSH) -I$(CLUSTER) -I$(FRECHET) -I$(BINARYTREE)

OBJ1= mainPart1.o mainLSH.o mainCube.o  $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingLSH.o $(PARSING)/parsingCube.o $(VECTOR)/vector.o $(LSH)/lsh.o $(HYPERCUBE)/hypercube.o $(HASHMAP)/hashmap.o $(LSH)/helperFunctions.o $(FRECHET)/discreteFrechet.o $(BINARYTREE)/binaryTree.o
OBJ2= mainCluster.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingCluster.o $(VECTOR)/vector.o $(HASHMAP)/hashmap.o $(LSH)/lsh.o $(LSH)/helperFunctions.o $(CLUSTER)/clustering.o $(CLUSTER)/clusterHelpingFuns.o $(CLUSTER)/kmeansPlusPlus.o $(HYPERCUBE)/hypercube.o $(FRECHET)/discreteFrechet.o $(BINARYTREE)/binaryTree.o
# OBJ3= $(FRECHERCONT)/simplification.o $(FRECHERCONT)/point.o $(FRECHERCONT)/interval.o $(FRECHERCONT)/frechet.o $(FRECHERCONT)/curve.o $(FRECHERCONT)/config.o

EXEC1 = search
EXEC2 = cluster

all: $(EXEC1) $(EXEC2)
#
# search: $(EXEC1)
#
# cluster: $(EXEC2)

$(EXEC1): $(OBJ1)
	$(CCPLUSPLUS) -c $(FRECHERCONT)/simplification.cpp $(FRECHERCONT)/point.cpp $(FRECHERCONT)/interval.cpp $(FRECHERCONT)/frechet.cpp $(FRECHERCONT)/curve.cpp $(FRECHERCONT)/config.cpp
	# $(CCPLUSPLUS) -o mainPart1.c
	# gcc -c -o lib1.o lib1/lib1.c
	# gcc -c -o lib2.o lib2/lib1.c
	# gcc -c -o lib_main.o lib_main.c
	# g++ -c -o main.o main.cpp
	# g++ -o main lib1.o lib2.o lib_main.o main.o
	$(CC) $(CFLAGS) $(OBJ1) -o $(EXEC1) -lm

$(EXEC2): $(OBJ2)
	$(CC) $(CFLAGS) $(OBJ2) -o $(EXEC2) -lm


.PHONY: clean

clean:
	rm -f $(OBJ1) $(OBJ2) $(EXEC1) $(EXEC2)
	# rm -f $(OBJ1) $(EXEC1)
