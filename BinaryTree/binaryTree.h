#ifndef BINARYTREE_H
#define BINARYTREE_H

typedef struct tree_head_node *Tree;

Tree createTreeFromList(List ,int );
void printTreeDFS(Tree );
void destroyTree(Tree ,int );

#endif
