#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../Vector/vector.h"
#include "../hashTable/hashTableList/hashTableList.h"
#include "../FrechetDistance/discreteFrechet.h"



typedef struct tree_n{
  Vector v;
  struct tree_n *left;
  struct tree_n *right;
}tree_node;  //list node
typedef tree_node *TreeNode;  //double pointer to the list node (array of lists)


typedef struct tree_head_node{
  int height;   // the size of the hashTable
  int count;    //how many keys/values are in the HashTable
  TreeNode root;   //the HashTable
}treeHeadNode;    // dummy node for hashTable
typedef treeHeadNode *Tree;    //pointer to the dummy node

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

TreeNode createTreeNode(Vector v){
  TreeNode newNode = malloc(sizeof(tree_node));
  newNode->v = copyVector(v);
  newNode->left=NULL;
  newNode->right=NULL;
  return newNode;
}

void deleteTreeNode(TreeNode tn,int freeVector){
  if(freeVector)
    deleteVector(getTnVector(tn));
  free(tn);
}

Tree initializeTree(){
  Tree tree = malloc(sizeof(treeHeadNode));
  tree->height=-1;
  tree->count=0;
  tree->root=NULL;
  return tree;
}

TreeNode recursiveTreeCreateFromList(List *list,int height,int *leafCount){
  if(height==0){
    if((*list)==NULL) return NULL;
    TreeNode tn = createTreeNode(getVector(*list));
    (*list) = getNext((*list));
    (*leafCount)++;
    return tn;
  }
  TreeNode tn = createTreeNode(NULL);
  tn->left = recursiveTreeCreateFromList(list,height-1,leafCount);
  tn->right = recursiveTreeCreateFromList(list,height-1,leafCount);
  return tn;
}


Tree createTreeFromList(List list,int count){
  if(count==0) return NULL;
  List tempList = list;
  int height = ceil(log2(count));
  Tree tree = initializeTree();
  tree->height = height;
  int leafCount = 0;
  tree->root = recursiveTreeCreateFromList(&tempList,height,&leafCount);
  tree->count=leafCount;
  return tree;
}

void printTreeRecursive(TreeNode tn,int current_height){
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
  int height = tree->height;
  printf("------------------------- PRINT TREE -------------------------\n");
  printf("ACTIVE LEAVES = %d\n",tree->count);
  for(int i=0;i<=height;i++){
    printf("* LEVEL %d:\n",i);
    printTreeRecursive(tree->root,i);
  }
}

void destroyTreeRecursive(TreeNode tn){
  if(tn==NULL) return;
  destroyTreeRecursive(getLeft(tn));
  destroyTreeRecursive(getRight(tn));
  if((getTnVector(tn)!=NULL))
    deleteVector(getTnVector(tn));
  free(tn);
}

void destroyTree(Tree tree){
  destroyTreeRecursive(getRoot(tree));
  free(tree);
}

void computeMeanCurvesRecursive(TreeNode tn){
  if(tn==NULL) return;
  computeMeanCurvesRecursive(getLeft(tn));
  computeMeanCurvesRecursive(getRight(tn));

  if(getLeft(tn)==NULL && getRight(tn)==NULL){
    return;
  }else if(getLeft(tn)==NULL){
    tn->v = copyVector(getTnVector(getRight(tn)));
  }else if(getRight(tn)==NULL){
    tn->v = copyVector(getTnVector(getLeft(tn)));
  }else{
    if(getTnVector(getRight(tn))==NULL && getTnVector(getLeft(tn))==NULL){
      return;
    }else if(getTnVector(getLeft(tn))==NULL){
      tn->v = copyVector(getTnVector(getRight(tn)));
    }else if (getTnVector(getRight(tn))==NULL){
      tn->v = copyVector(getTnVector(getLeft(tn)));
    }else{
      tn->v = meanCurveBetween2Curves(getTnVector(getLeft(tn)),getTnVector(getRight(tn)));
    }
  }
  return;
}

Vector treeFindMeanCurve(Tree tree){
  if(tree==NULL) {printf("TREE EMPTY!\n"); return NULL; }
  computeMeanCurvesRecursive(getRoot(tree));
  printf("computeMeanCurvesRecursive OK!!\n");
  return copyVector(getTnVector(getRoot(tree)));
}
