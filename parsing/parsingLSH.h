#ifndef PARSING_H
#define PARSING_H


void readFileLSH(char*,List * ,int *, int ,double *,int);
void readQueryFileLSH(char*,char*,LSH,List,int);
void readQueryFileLSH_DiscreteFrechet(char* ,char* ,LSH ,List ,Grids ,double ,double *,int);
void readQueryFileLSH_ContinuousFrechet(char* ,char* ,LSH ,List ,double ,double ,double *,int ,Grids );
int findDimLSH(char* );

#endif
