#ifndef HASHMAP_H
#define HASHMAP_H

typedef int Key;
typedef int Value;

typedef struct node *Record;
typedef struct hash *HashMap;


Key getKey(Record );
Value getValue(Record );

HashMap hmCreate(int );

Record hmSearchOrInsert(HashMap ,Key ,Value );

Record hmSearch(HashMap ,Key );

void hmDestroy(HashMap );



#endif
