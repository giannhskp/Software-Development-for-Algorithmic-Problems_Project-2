#include <stdlib.h>
#include "my_interface.hpp"
#include "point.hpp"
#include "curve.hpp"
#include "frechet.hpp"

double compute_continuous_distance(double *y1,double *x1, double *y2, double *x2,int d1,int d2){

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
  Curve c1 = Curve(p1);

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
  Curve c2 = Curve(p2);

  Frechet::Continuous::Distance d = Frechet::Continuous::distance(c1,c2); // call library's distance() function

  return d.value;
}
