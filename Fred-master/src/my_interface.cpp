#include <stdlib.h>
#include "my_interface.hpp"
#include "point.hpp"
#include "curve.hpp"
#include "frechet.hpp"

double compute_continuous_distance(double *y1,double *x1, double *y2, double *x2,int d1,int d2){
  // used to shape the correspodings objects of class Curve that represent the two timeseries
  // and compute the Continuous Frechet Distance between them using the library function,
  // this distance will be returned by this function.
  // this function is essentially the interface between the C++ library and the C code.
  // the first two tables given as arguments correspond to the first time serie
  // the other two to the second one


  // initialize an object point with 2 dimensions for the first time serie
  Points p1 =  Points(2);
  for(int i=0;i<d1;i++){
    if(y1[i]<0){
      continue;
    }
    Point tempPoint = Point(2);
    tempPoint.set(0,x1[i]);
    tempPoint.set(1,y1[i]);
    p1.add(tempPoint);
  }
  // form the first timeseries
  Curve c1 = Curve(p1);

  // initialize an object point with 2 dimensions for the second time serie
  Points p2 =  Points(2);
  for(int i=0;i<d2;i++){
    if(y2[i]<0){
      continue;
    }
    Point tempPoint = Point(2);
    tempPoint.set(0,x2[i]);
    tempPoint.set(1,y2[i]);
    p2.add(tempPoint);
  }
  // form the second timeseries
  Curve c2 = Curve(p2);

  Frechet::Continuous::Distance d = Frechet::Continuous::distance(c1,c2); // call library's distance() function

  return d.value;
}
