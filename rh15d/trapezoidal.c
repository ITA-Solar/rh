#include<stdio.h>
#include<math.h>

double trapezoidal(double *x, double *y, int n){
  double h,sum=0,integral;
  int i;

//  h = x[1]-x[0];
  
//  for(i=1;i<n;i++){ 
 //   sum=sum+y[i];
  //}

  for (i = 0; i<n-1;i++){
       sum = sum + (x[i+1]-x[i]) * (y[i+1]+y[i])/2.0;
      }
  // integral=(h/2)*(y[0]+y[n-1]+2*sum);
      integral = sum;
//      printf("\nintegral = %.10f",integral);
  return integral;
}
