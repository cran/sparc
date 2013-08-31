#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"

void SGLM(double *A, double *X, double *y, int *NN, int *nn, int *dd, double *lambda){
    
    int N,n,d;
    
    N = *NN;
    n = *nn;
    d = *dd;
    
    int i,j,k,t,blocki,blockj,blockk;
    
    for(t=0;t<d;t++){
        
        k = 0;
        while(k<N){
            for(i=0;i<(n-1);i++)
                for(j=i+1;j<n;j++){
                    A[t*N+k] = (X[t*n+i] - X[t*n+j])*(y[j]-y[i]);
                    k++;
                }
            k++;
        }
        
        lambda[t] = 0;
        for(k=0;k<N;k++)
            lambda[t] = lambda[t] + A[t*N+k];
        
        lambda[t] = lambda[t]/2;
    }
}
