#ifndef PARSING_H
#define PARSING_H


void readFileLSH(char*,List * ,int *, int ,double *);
void readQueryFileLSH(char*,char*,LSH,List);
void readQueryFileLSH_DiscreteFrechet(char* ,char* ,LSH ,List ,Grids ,double ,double *);
void readQueryFileLSH_ContinuousFrechet(char* ,char* ,LSH ,List ,double ,double ,double *);
int findDimLSH(char* );

#endif
