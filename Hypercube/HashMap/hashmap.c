#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hashmap.h"
#include "../../LSH/helperFunctions.h"




typedef struct node{
  Key key;  // the key of the map (the type  of the key is defined in hashmap.h file)
  Value value;  // the value that correspons to the key (the type  of the value is defined in hashmap.h file)
  struct node *next;
}Node;  //list node
typedef Node *Record;    //pointer to the list node
typedef Record *HashMapArray;  //double pointer to the list node (array of lists)

typedef struct hash{
  int size;   // the size of the hashTable
  int count;    //how many keys/values are in the HashTable
  HashMapArray table;   //the HashTable
}HashHeadNode;    // dummy node for hashTable
typedef HashHeadNode *HashMap;    //pointer to the dummy node


Key getKey(Record r){
  return r->key;
}
Value getValue(Record r){
  return r->value;
}

HashMap hmCreate(int size){ //  create/initialize the hash table
  HashMap head=malloc(sizeof(HashHeadNode));  //Create hash table dummy node
  head->size=size;    //store the given size of the hash table
  head->table=malloc(head->size*(sizeof(Record)));   //create a hash table of the given size
  for(int i=0;i<head->size;i++){
    head->table[i]=NULL;    //initialize hash table with null in every space (lists)
  }
  head->count=0;
  return head;
}

Record newNode(Key key,Value value){   // creates a node of a key/value set
  Record node=malloc(sizeof(Node));
  node->key=key;
  node->value=value;
  node->next=NULL;
  return node;
}

void deleteNode(Record node){  // delete a key/value node
  free(node);
}

unsigned int hash(HashMap ht,char *id){      //find the hash table index for the given id
  unsigned long int hashcode = 0;
  int a=33;
  char *tempkey=id;
  for (; *tempkey!='\0'; tempkey++)   // for every char/number of the id
    hashcode=(a*hashcode + *tempkey) % ht->size;
  return hashcode;    //hash code is beiing %hashTableSize so it is 0<= hashcode <= ht->size
}


void hmResizeInsert(HashMap ,Record );

void hmResize(HashMap ht){  // the hash table is resizable, when it is 90% full this function is called
  // in order to double the size of the hash table, so as to retain the search/insert in O(1) time
  HashMapArray oldTable=ht->table;  // store the old hashTable
  int oldSize=ht->size; // store the size of the old hashTable
  int newSize=2*(ht->size); // the new size is double the preivous one
  HashMapArray newTable=malloc(newSize*(sizeof(Record)));  //create the new hashTable
  for(int i=0;i<newSize;i++){ // initialize it
    newTable[i]=NULL;
  }
  ht->table=newTable; // refresh the dummy node with the new hash table
  ht->size=newSize;   // and the new size
  ht->count=0;
  // insert all the records of the old hash table to the new hash table
  for (int i = 0; i < oldSize; ++i){   //for every bucket of the old hash table
    Record slot=oldTable[i];    //the "overflow list" of this bucket
    if(slot==NULL){
        continue;
    }
    while(slot!=NULL){   //for every node (every key/value set) of the overflow list
        Record next=slot->next;
        hmResizeInsert(ht,slot);  // insert it to the new hash table
        // a new insert function is created for the resize so as not to delete every node of the old hash table
        // and create a new node for the new hash table (see below: htResizeInsert)
        slot=next;
    }
  }
  free(oldTable);
}

void hmResizeInsert(HashMap ht,Record rec){
  // this function is created in order to avoid extra memmory allocation and frees
  // with this functions the same nodes from the old hash table are used for the new hash table
  // in this way no extra allocation are done during the resize of the hash table

  // int index=hash(ht,rec->citizenId);
  int index = mod_Int_Int(rec->key,ht->size);
  Record slot=ht->table[index];
  if(slot==NULL){   //no collitions
    ht->table[index]=rec; //insert it to the hash table, as the first (and only) node of the list in this "Bucket"
    rec->next=NULL;
    (ht->count)++;
    return;
  }else{    //if there are already key/value sets in this "bucket"
    Record prev;
    while(slot!=NULL){
      if(getKey(slot)==rec->key){    //check if the key already exists
        return;
      }
      prev=slot;
      slot=slot->next;
    }
    //when we reach the last node of the list (the key does not exist)
    prev->next=rec;    //connect it to the last node of the list
    rec->next=NULL;
    (ht->count)++;
    return;
  }
}

Record hmSearchOrInsert(HashMap ht,Key key,Value value){
  if((double)(ht->count)>(0.9*((double)((ht->size))))){ // if the hash table is 90% full, resize it in order to keep the search/insert O(1)
      hmResize(ht);
  }
  int index = mod_Int_Int(key,ht->size);
  Record slot=ht->table[index];  // get the overflow list of this bucket
  if(slot==NULL){   //no collitions so just insert the node with the citizen
    Record newnode=newNode(key,value);     //create the node of the hashTable list that contains the given key/value set
    ht->table[index]=newnode; //insert it to the hash table, as the first (and only) node of the list in this "Bucket"
    (ht->count)++;
    return newnode;
  }else{    //if there are already key/value sets in this "bucket"
    Record prev;
    while(slot!=NULL){
      if(getKey(slot)==key){    //check if the key already exists
        return slot;
      }
      prev=slot;
      slot=slot->next;
    }
    //when we reach the last node of the list (the key does not exist)
    Record newnode=newNode(key,value); //create the node of the hashTable list that contains the given key/value set
    prev->next=newnode;    // and connect it to the last node of the list
    (ht->count)++;
    return newnode;
  }
}

Record hmSearch(HashMap ht,Key key){
  // int index=hash(ht,citizenId);
  int index = mod_Int_Int(key,ht->size);  //find in wich bucket the key/value set should be (using hash function and key)
  Record slot=ht->table[index];    //go to this "bucket"
  if(slot==NULL){   //no citizen exists in this slot, so it is sure that citizenId does not exist in the hash table
    return NULL;
  }else{
    while(slot!=NULL){    //check every key in this slot (all the collitions)
      if(getKey(slot)==key){   //if the key is found, return a pointer to the node
          return slot;   //we found the key so return a pointer to the citizen
      }
      slot=slot->next;
    }
    //if we traverse the whole list in this bucket and we dont find the key, the key does not exist in the hash table
    return NULL;
  }
}

void hmDestroy(HashMap ht){   //destroys and frees the hash table (and the key/value nodes)
  if(ht==NULL)
    return;
  for (int i=0;i<ht->size;++i){   //for every bucket in the hash table
    Record slot=ht->table[i];    //take the key/value set list
    if(slot==NULL){
        continue;
    }
    for(;;){   //every node in the list
        if(slot->next==NULL){   //for the last node
            deleteNode(slot);     // deleteNode()  deletes the key/value info and the list node itself
            break;
        }
        Record next=slot->next;
        deleteNode(slot);
        slot=next;
    }
  }
  free(ht->table);    //free the board
  free(ht);   //free the dummy node of the hash table
}
