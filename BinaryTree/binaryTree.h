#ifndef BINARYTREE_H
#define BINARYTREE_H

typedef struct tree_head_node *Tree;

Tree createTreeFromList(List ,int );
Tree createTreeFromHt(HashTable ,int );
void printTreeDFS(Tree );
void destroyTree(Tree );
Vector treeFindMeanCurve(Tree );
#endif
