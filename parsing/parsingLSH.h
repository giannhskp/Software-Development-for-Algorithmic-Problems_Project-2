#ifndef PARSING_H
#define PARSING_H


void readFileLSH(char*,List * ,int *);
void readQueryFileLSH(char*,char*,LSH,List);
void readQueryFileLSH_DiscreteFrechet(char* ,char* ,LSH ,List ,Grids ,Vector ,double );
void readQueryFileLSH_ContinuousFrechet(char* ,char* ,LSH ,List ,Vector ,double ,double );
int findDimLSH(char* );

#endif
