#ifndef PARSING_H
#define PARSING_H


void readFileLSH(char*,List * ,int *);
void readQueryFileLSH(char*,char*,LSH,List);
void readQueryFileLSH_DiscreteFrechet(char* ,char* ,LSH ,List ,Grids ,Vector ,double );
int findDimLSH(char* );

#endif
