#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"

void get_ep(double *ep, double *A, double *w, int *nn, int *dd, int *idx0, int *ssize0){
    
    int i,j,n,d,block,size0;
    double *p;
    n = *nn;
    d = *dd;
    size0 = *ssize0;
    
    for(i=0;i<n;i++)
        ep[i] = 0;
    
    for(j=0;j<size0;j++){
        block = idx0[j]*n;
        for(i=0;i<n;i++){
            ep[i] = ep[i] + A[block+i]*w[idx0[j]];
        }
    }
    
    for(i=0;i<n;i++)
        ep[i] = exp(ep[i]);
}

void get_gradient(double *g, double *ep, double *A, int *nn, int *dd){
    
    int i,j,n,d,block;
    n = *nn;
    d = *dd;
    
    for(j=0;j<d;j++){
        block = j*n;
        g[j] = 0;
        for(i=0;i<n;i++)
            g[j] = g[j] + ep[i]/(1+ep[i])*A[block+i];
    }
}

void get_obj(double *obj, double *ep, int *nn){
    
    int i,n;
    n = *nn;
    
    *obj = 0;
    for(i=0;i<n;i++){
        *obj = *obj + log(ep[i]+1);
    }
}

void get_sol(double *w1, double *w0, double *g, double *LL, double *iilambda, int *dd, int *idx0, int *ssize0){
    int i,j,k,d,size0;
    double ilambda,L;
    d = *dd;
    L = *LL;
    *ssize0 = 0;
    ilambda = *iilambda/L;
    
    for(j=0;j<d;j++){
        w1[j] = w0[j] - g[j]/L;
        if(w1[j]>=ilambda){
            w1[j] = w1[j] - ilambda;
            idx0[*ssize0] = j;
            *ssize0 = *ssize0 + 1;
        }
        else if(w1[j] <= -ilambda){
            w1[j] = w1[j] + ilambda;
            idx0[*ssize0] = j;
            *ssize0 = *ssize0 + 1;
        }
        else{
            w1[j] = 0;
        }
    }
}


void get_Q(double *Q, double *Q0, double *g, double *w0, double *w1, double *LL, int *dd){
    int j,d;
    double L,tmp;
    d = *dd;
    L = *LL;

    *Q = *Q0;
    
    for(j=0;j<d;j++){
        tmp = w1[j] - w0[j];
        *Q = *Q + g[j]*tmp + L*pow(tmp,2)/2;
    }
}

void get_gap(double *gap, double *g, double *iilambda, int *dd){
    int j,d;
    double tmp,ilambda;
    d = *dd;
    ilambda = *iilambda;
    
    *gap = 0;
    for(j=0;j<d;j++){
        tmp = fabs(g[j]);
        if(tmp>(*gap))
            *gap = tmp;
    }
    
    //printf("%f\n",*gap);
    
    *gap = *gap - ilambda;
}

void SLRM(double *A, double *lambda, int *nnlambda, double *LL0, int *nn, int *dd, double *ww, int *mmax_ite, double *tthol, double *aalpha){
    
    int m,n,d,max_ite,nlambda,block;
    double thol,ilambda,L0,alpha;
    double L;
    double Q,H,Q0;
    n = *nn;
    d = *dd;
    max_ite = *mmax_ite;
    thol = *tthol;
    L0 = *LL0;
    alpha = *aalpha;
    nlambda = *nnlambda;
    
    
    int ite;
    int i,j;
    
    double gap;
    int tracking;
    
    int *idx0;
    int size0;
    idx0 = (int *) malloc(d*sizeof(int));
    
    size0 = 0;
    //for(j=0;j<d;j++){
    //    if(ww[j]!=0){
    //        idx0[size0] = j;
    //        size0 = size0 + 1;
    //    }
    //}
    
    //for(j=0;j<size0;j++)
    //    printf("%d\n",idx0[j]);
    
    double *ep,*g,*w0,*w1;
    ep = (double *) malloc(n*sizeof(double));
    g = (double *) malloc(d*sizeof(double));
    w0 = (double *) malloc(d*sizeof(double));
    w1 = (double *) malloc(d*sizeof(double));
    
    for(j=0;j<d;j++){
        w0[j] = 0;
    }
    
    L = L0;
    
    for(m=0;m<nlambda;m++){
   
        ilambda = lambda[m];
        
        get_ep(ep,A,w0,&n,&d,idx0,&size0);
    
        //for(i=0;i<n;i++)
        //    printf("%f\n",ep[i]);
    
        get_gradient(g,ep,A,&n,&d);
    
        //for(j=0;j<d;j++)
        //    printf("%f\n",g[j]);
    
        get_obj(&Q0,ep,&n);
    
        //printf("%f\n",Q0);
    
        //printf("%f\n",L0);
    
    
        tracking = 1;
    
        while(tracking == 1){
        
            get_sol(w1,w0,g,&L,&ilambda,&d,idx0,&size0);
    
            //for(j=0;j<d;j++)
            //    printf("%f\n",w1[j]);

            get_Q(&Q,&Q0,g,w0,w1,&L,&d);
        
            //printf("%f\n",Q);
        
            get_ep(ep,A,w1,&n,&d,idx0,&size0);
        
            //for(i=0;i<n;i++)
            //    printf("%f\n",ep[i]);
        
            get_obj(&H,ep,&n);
        
            //printf("%f\n",H);
        
            if(Q > H)
                L = L*alpha;
            else{
                L = L/alpha;
                tracking = 0;
            }
            //tracking = 0;
        }
    
        //printf("%f\n",Q);
        //printf("%f\n",H);
        //printf("%f\n",L);
    
    
        get_sol(w1,w0,g,&L,&ilambda,&d,idx0,&size0);
    
    
        //for(j=0;j<d;j++)
        //    printf("%f\n",w1[j]);
    
        for(j=0;j<d;j++)
            w0[j] = w1[j];
    
        gap = 99999999;
        ite = 1;

        while((gap>(thol*lambda[nlambda-1]))&&(ite<max_ite)){
        
            get_ep(ep,A,w0,&n,&d,idx0,&size0);
        
            get_gradient(g,ep,A,&n,&d);
    
            //for(j=0;j<d;j++)
                //printf("%f\n",g[j]);
    
            get_gap(&gap,g,&ilambda,&d);
    
            //printf("%f\n",gap);
    
            tracking = 1;
            
            get_obj(&Q0,ep,&n);
            
            while(tracking == 1){
                
                get_sol(w1,w0,g,&L,&ilambda,&d,idx0,&size0);
                
                get_Q(&Q,&Q0,g,w0,w1,&L,&d);
                
                get_ep(ep,A,w1,&n,&d,idx0,&size0);
                
                get_obj(&H,ep,&n);
                
                if(Q>H)
                    tracking = 0;
                else
                    L = L/alpha;
            }
        
            ite = ite + 1;
        
            for(j=0;j<d;j++)
                w0[j] = w1[j];
        }
        
        block = m*d;
        for(j=0;j<size0;j++)
            ww[block + idx0[j]] = w1[idx0[j]];
    
        //printf("%d\n",block);
        
    }
    
    free(idx0);
    free(ep);
    free(g);
    free(w0);
    free(w1);
}
