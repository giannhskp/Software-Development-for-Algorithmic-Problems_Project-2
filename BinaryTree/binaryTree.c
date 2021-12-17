#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTable.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../FrechetDistance/discreteFrechet.h"



typedef struct tree_n{
  Vector v; // curve that is stored in the tree node
  struct tree_n *left;  // pointer to left child
  struct tree_n *right;  // pointer to right child
}tree_node;
typedef tree_node *TreeNode;


typedef struct tree_head_node{
  int height; // height of the tree
  int count;  // how many curves where inserted at te tree creation
  TreeNode root;  // pointer to the root tree node
}treeHeadNode;
typedef treeHeadNode *Tree;

TreeNode getRoot(Tree tree){
  return tree->root;
}

TreeNode getLeft(TreeNode tn){
  return tn->left;
}

TreeNode getRight(TreeNode tn){
  return tn->right;
}

Vector getTnVector(TreeNode tn){
  return tn->v;
}

TreeNode createTreeNode(Vector v){  // create a tree node
  TreeNode newNode = malloc(sizeof(tree_node)); // allocate space for the node
  newNode->v = copyVector(v); // copy the given curve and store it at the node
  newNode->left=NULL; // assign left child to NULL
  newNode->right=NULL; // assign left child to NULL
  return newNode;
}

void deleteTreeNode(TreeNode tn,int freeVector){
  if(freeVector)
    deleteVector(getTnVector(tn));
  free(tn);
}

Tree initializeTree(){  // initialize an empty tree
  Tree tree = malloc(sizeof(treeHeadNode));
  tree->height=-1;
  tree->count=0;
  tree->root=NULL;
  return tree;
}

TreeNode recursiveTreeCreateFromList(List *list,int height,int *leafCount){
  // recursively creates a tree that contains all the curves in the given list at it's leafs
  if(height==0){  // if we reached at the leaf level
    if((*list)==NULL) return NULL;  // if all the curves of the list were assigned to a leaf, dont create a node
    TreeNode tn = createTreeNode(getVector(*list)); // create a new node, and assing the curve to the node
    (*list) = getNext((*list)); // got the next curve that will be assigned to a leaf
    (*leafCount)++;
    return tn;
  }
  // if we are at an inner node
  // inner nodes are empty in the beggining, as the mean curves will be stored afterwards
  TreeNode tn = createTreeNode(NULL); // create an empty node
  tn->left = recursiveTreeCreateFromList(list,height-1,leafCount);  // create left child at the next level, height is reduced by 1
  tn->right = recursiveTreeCreateFromList(list,height-1,leafCount);  // create right child at the next level, height is reduced by 1
  return tn;
}


Tree createTreeFromList(List list,int count){
  // create a tree that contains all the curves of the list at it's leaves
  if(count==0) return NULL; // if list has no items
  List tempList = list;
  int height = ceil(log2(count)); // compute the height of the tree that is needed in order to store all the curves at the leaves
  Tree tree = initializeTree(); // initialize an empty tree
  tree->height = height;  // store the height of the tree
  int leafCount = 0;
  tree->root = recursiveTreeCreateFromList(&tempList,height,&leafCount);  // create all the tree nodes
  tree->count=leafCount;  // store the count of the leaves that contain a curve
  return tree;
}

TreeNode recursiveTreeCreateFromHt(List **listArray,int *arrayIndex,int height,int *leafCount){
  // create a tree that contains all the curves of the array of lists at it's leaves
  if(height==0){    // if we reached at the leaf level
    if((*arrayIndex)<0 && ((*listArray)[0]==NULL))  // if we are at the last index of the array and at the end of the list
      return NULL;
    while((*listArray)[(*arrayIndex)]==NULL){ // while the index of the array contains an empty list, go to the next array index
      (*arrayIndex)--;
      if((*arrayIndex)<0 )  // if we reached at the end of the array, stop
        return NULL;
    }
     // create a node that contains the vector of the current list node on the current array index
    TreeNode tn = createTreeNode(getVector((*listArray)[(*arrayIndex)]));
    // go to the next node of the list
    (*listArray)[(*arrayIndex)] = getNext((*listArray)[(*arrayIndex)]);
    (*leafCount)++;
    return tn;
  }
  // if we are at an inner node
  // inner nodes are empty in the beggining, as the mean curves will be stored afterwards
  TreeNode tn = createTreeNode(NULL);
  tn->left = recursiveTreeCreateFromHt(listArray,arrayIndex,height-1,leafCount);   // create left child at the next level, height is reduced by 1
  tn->right = recursiveTreeCreateFromHt(listArray,arrayIndex,height-1,leafCount);  // create right child at the next level, height is reduced by 1
  return tn;
}

Tree createTreeFromHt(HashTable ht,int count){
  // create a tree that contains all the curves of the given hash table at it's leaves
  // every bucket of the hash table contains a list
  // so we keep all the lists to an array, an create a tree by getting the curves from this array
  if(count==0) return NULL;
  int numOfBuckets = getNumberOfBuckets(ht);
  List *listArray=malloc(numOfBuckets*sizeof(List));  // create an array with size same as the number of buckets
  for(int i=0;i<numOfBuckets;i++){
    listArray[i] = getListOfBucket(ht,i); // every array index contains a list
  }
  int height = ceil(log2(count)); // compute the height of the tree that is needed in order to store all the curves at the leaves
  Tree tree = initializeTree(); // initialize an empty tree
  tree->height = height;    // store the height of the tree
  int leafCount = 0;
  int starting_index = numOfBuckets-1;  // start from the last index of the array
  tree->root = recursiveTreeCreateFromHt(&listArray,&starting_index,height,&leafCount);   // create all the tree nodes
  tree->count=leafCount;    // store the count of the leaves that contain a curve
  free(listArray);  // free the array
  return tree;
}

void printTreeRecursive(TreeNode tn,int current_height){
  // prints the tree | for debbuging reasons
  if(tn==NULL) return;
  if(current_height==0){
    if(getTnVector(tn)==NULL) return;
    printVectorId(getTnVector(tn));
    return;
  }
  printTreeRecursive(getLeft(tn),current_height-1);
  printTreeRecursive(getRight(tn),current_height-1);
  return;
}

void printTreeDFS(Tree tree){
    // prints the tree | for debbuging reasons
  int height = tree->height;
  printf("------------------------- PRINT TREE -------------------------\n");
  printf("ACTIVE LEAVES = %d\n",tree->count);
  for(int i=0;i<=height;i++){
    printf("* LEVEL %d:\n",i);
    printTreeRecursive(tree->root,i);
  }
}

void destroyTreeRecursive(TreeNode tn){
  // destroys the tree and free's all the allocated memory
  if(tn==NULL) return;
  destroyTreeRecursive(getLeft(tn));  // delete left child
  destroyTreeRecursive(getRight(tn));  // delete right child
  if((getTnVector(tn)!=NULL)) // if it contains a curve
    deleteVector(getTnVector(tn));  // delete the curve
  free(tn); // delete the node
}

void destroyTree(Tree tree){
  // destroys the tree and free's all the allocated memory
  destroyTreeRecursive(getRoot(tree));
  free(tree);
}

void computeMeanCurvesRecursive(TreeNode tn){
  if(tn==NULL) return;
  // DFS on the tree
  // after the DFS the mean final curve will be stored to the root node
  computeMeanCurvesRecursive(getLeft(tn));  // compute mean curve of the left child
  computeMeanCurvesRecursive(getRight(tn));  // compute mean curve of the right child

  if(getLeft(tn)==NULL && getRight(tn)==NULL){  // if the current node has no childs, is a leaf node
    return; // return as the leaf node either contains a input curve or nothing (last remaining leaves of last layer)
  }else if(getLeft(tn)==NULL){  // if the current node only has right child
    tn->v = copyVector(getTnVector(getRight(tn)));  // assing to the current node the curve of the right child
  }else if(getRight(tn)==NULL){  // if the current node only has left child
    tn->v = copyVector(getTnVector(getLeft(tn)));  // assing to the current node the curve of the left child
  }else{  // if the current node hash both child nodes
    if(getTnVector(getRight(tn))==NULL && getTnVector(getLeft(tn))==NULL){  // check if both child nodes contain a curve
      return;
    }else if(getTnVector(getLeft(tn))==NULL){ // if the left node does not contain a curve
      tn->v = copyVector(getTnVector(getRight(tn)));  // assing to the current node the curve of the right child
    }else if (getTnVector(getRight(tn))==NULL){ // if the right node does not contain a curve
      tn->v = copyVector(getTnVector(getLeft(tn)));  // assing to the current node the curve of the left child
    }else{  // if both childs contain curves
      // compute the mean curve between the two curves of the childs, and store the mean curve to the current node
      tn->v = meanCurveBetween2Curves(getTnVector(getLeft(tn)),getTnVector(getRight(tn)));
    }
  }
  return;
}

Vector treeFindMeanCurve(Tree tree){
  // called by the cluster function in order to compute the mean curve of the curves that are inserted at the leaf nodes of the given tree
  if(tree==NULL) { return NULL; }
  computeMeanCurvesRecursive(getRoot(tree));
  return copyVector(getTnVector(getRoot(tree)));
}
