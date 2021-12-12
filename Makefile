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
CCPLUSPLUS=g++ -std=c++14
CFLAGS= -g -Wall -I$(HASHTABLELIST) -I$(HASHTABLE) -I$(PARSING) -I$(VECTOR) -I$(LSH) -I$(CLUSTER) -I$(FRECHET) -I$(BINARYTREE) -I$(FRECHERCONT)

OBJ1= mainPart1.o mainLSH.o mainCube.o  $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingLSH.o $(PARSING)/parsingCube.o $(VECTOR)/vector.o $(LSH)/lsh.o $(HYPERCUBE)/hypercube.o $(HASHMAP)/hashmap.o $(LSH)/helperFunctions.o $(FRECHET)/discreteFrechet.o $(BINARYTREE)/binaryTree.o
OBJ2= mainCluster.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(PARSING)/parsingCluster.o $(VECTOR)/vector.o $(HASHMAP)/hashmap.o $(LSH)/lsh.o $(LSH)/helperFunctions.o $(CLUSTER)/clustering.o $(CLUSTER)/clusterHelpingFuns.o $(CLUSTER)/kmeansPlusPlus.o $(HYPERCUBE)/hypercube.o $(FRECHET)/discreteFrechet.o $(BINARYTREE)/binaryTree.o
OBJ3= my_interface.o simplification.o point.o interval.o frechet.o curve.o config.o
OBJ4= unitTesting.o  mainLSH.o mainCube.o $(LSH)/lsh.o $(HYPERCUBE)/hypercube.o $(PARSING)/parsingLSH.o $(PARSING)/parsingCube.o $(HASHTABLE)/hashTable.o $(HASHTABLELIST)/hashTableList.o $(HASHMAP)/hashmap.o $(VECTOR)/vector.o  $(LSH)/helperFunctions.o $(FRECHET)/discreteFrechet.o $(BINARYTREE)/binaryTree.o


EXEC1 = search
EXEC2 = cluster
EXEC3 = unitTesting

all: $(EXEC1) $(EXEC2)
#
# search: $(EXEC1)
#
testing: $(EXEC3)

$(EXEC1): $(OBJ1)
	$(CCPLUSPLUS) -c $(FRECHERCONT)/my_interface.cpp $(FRECHERCONT)/simplification.cpp $(FRECHERCONT)/point.cpp $(FRECHERCONT)/interval.cpp $(FRECHERCONT)/frechet.cpp $(FRECHERCONT)/curve.cpp $(FRECHERCONT)/config.cpp
	$(CC) $(CFLAGS) $(OBJ1) $(OBJ3) -o $(EXEC1) -lm -lstdc++

$(EXEC2): $(OBJ2)
	$(CC) $(CFLAGS) $(OBJ2) $(OBJ3) -o $(EXEC2) -lm -lstdc++

$(EXEC3): $(OBJ4)
	$(CCPLUSPLUS) -c $(FRECHERCONT)/my_interface.cpp $(FRECHERCONT)/simplification.cpp $(FRECHERCONT)/point.cpp $(FRECHERCONT)/interval.cpp $(FRECHERCONT)/frechet.cpp $(FRECHERCONT)/curve.cpp $(FRECHERCONT)/config.cpp
	$(CC) $(CFLAGS) $(OBJ4) $(OBJ3) -o $(EXEC3) -lm -lstdc++ -lcunit


.PHONY: clean

clean:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(EXEC1) $(EXEC2)
	# rm -f $(OBJ1) $(EXEC1)
