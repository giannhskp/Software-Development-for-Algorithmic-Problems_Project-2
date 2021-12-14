#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "./Vector/vector.h"
#include "./hashTable/hashTableList/hashTableList.h"
#include "./LSH/lsh.h"
#include "./FrechetDistance/discreteFrechet.h"
//////////////////////////////////////////////////
#include "CUnit/Basic.h"

typedef struct gfunc g_function;

int init_suite1(void){
   return 0;
}

/* The suite cleanup function.
* Closes the temporary file used by the tests.
* Returns zero on success, non-zero otherwise.
*/
int clean_suite1(void){
  return 0;
}


int k_LSH=5;
int hashTableSize=100;
int w=6;


void test_gFunction(void){
  int dim=5;
  double coords1[5]={3.53,7.56,13.34,54.33,10.02};
  double coords2[5]={3.4,7.7,13.5,54.5,10};
  Vector vec1=initVector(coords1,"1",dim);
  Vector vec2=initVector(coords2,"2",dim);
  unsigned int temp=-1;
  LSH lsh = initializeLSH(5,dim);
  getValueOfFirstGFun(lsh,vec1,&temp);
  CU_ASSERT(getValueOfFirstGFun(lsh,vec1,&temp)==getValueOfFirstGFun(lsh,vec2,&temp));
  deleteVector(vec1);
  deleteVector(vec2);
  destroyLSH(lsh);

}

void test_discreteFrechet(void){
  int dim=5;
  double coords1[5]={5.53,75.6,3.34,54.3,10.02};
  double coords2[5]={2.13,67.06,1.44,94.3,31.3};
  double sum=0.0;
  double time[dim];
  for(int i=0;i<dim;i++){
    time[i]=sum;
    sum+=1.0;
  }
  Vector timeSeries1=initTimeSeries(coords1,time,"1",dim);
  Vector timeSeries2=initTimeSeries(coords2,time,"2",dim);
  CU_ASSERT(40.00 == discreteFrechet(timeSeries1,timeSeries2));
  deleteVector(timeSeries1);
  deleteVector(timeSeries2);
}

void test_meanCurveBetween2Curves(void){
  int dim=5;
  double coords1[5]={5.53,75.6,3.34,54.3,10.02};
  double coords2[5]={2.13,67.06,1.44,94.3,31.3};
  double meanCoords[6]={3.83,71.33,2.390,74.3,42.8,20.66};
  double meanTime[6]={0.0,1.0,2.0,3.0,3.5,4.0};
  double sum=0.0;
  double time[dim];
  for(int i=0;i<dim;i++){
    time[i]=sum;
    sum+=1.0;
  }
  Vector timeSeries1=initTimeSeries(coords1,time,"1",dim);
  Vector timeSeries2=initTimeSeries(coords2,time,"2",dim);
  Vector realMean=initTimeSeries(meanCoords,meanTime,"real",6);
  Vector test = meanCurveBetween2Curves(timeSeries1,timeSeries2);
  CU_ASSERT(compareTimeSeries(test,realMean) == 1);


  deleteVector(timeSeries1);
  deleteVector(timeSeries2);
  deleteVector(realMean);
  deleteVector(test);
}



char *distanceMetric;

int main(int argc, char *argv[]) {
  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS != CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", init_suite1, clean_suite1);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  /* add the tests to the suite */
  if ((NULL == CU_add_test(pSuite, "test of Discrete Frechet Distance", test_discreteFrechet)) ||
     (NULL == CU_add_test(pSuite, "test of Mean Curve Between 2 Curves", test_meanCurveBetween2Curves) ||
      NULL == CU_add_test(pSuite, "test of G function of LSH", test_gFunction) )){
    CU_cleanup_registry();
    return CU_get_error();
  }



  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}