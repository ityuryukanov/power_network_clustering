/**********
 *  
 *  mex interface to compute the gradient of the eigenvectors 
 *  alignment quality
 *
 *  To mexify:   
 *               mex  evrot.cpp;
 *   
 *  [clusters,Quality,Vrot] = evrot(V,method);
 *   
 *  Input:
 *    V = eigenvecors, each column is a vector
 *    method = 1   gradient descent
 *             2   approximate gradient descent
 *
 *  Output:
 *    clusts - Resulting cluster assignment
 *    Quality = The final quality
 *    Vr = The rotated eigenvectors
 *
 * 
 *  Lihi Zelnik (Caltech) March 2005 
 *  
 * 
 ************/

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEBUG 0    /* set to 1 and mex to see print outs */
#define  EPS 2.2204e-16 


/**** DEBUG utility *****************************************/
void print_array(mxArray *A)
{
    int i,j,ind;
    
    double *p_A = mxGetPr(A);
    const size_t *idims =  mxGetDimensions(A);

    mexPrintf("\n");
    ind = 0;
    for( j=0; j < idims[1]; j++ ){
        mexPrintf("Column %d:     ",j);
        for( i=0; i < idims[0]; i++ ){
            mexPrintf("%g ",p_A[ind]);
            ind++;
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
}


/************************************************************/
/*** multiply matrices    ***/
mxArray* matrix_mult(mxArray *A, mxArray *B)
{
    mxArray  *lhs[1], *rhs[2];
                                                                                                                              
    rhs[0] = A;
    rhs[1] = B;
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=A*B */
                                                                                                                              
    return lhs[0];
}
        
/*** compute the matrix A given X,U1,V,U2    ***/
mxArray* buildA(mxArray *X, mxArray *U1, mxArray *Vk, mxArray *U2)
{
    mxArray  *lhs[1], *rhs[2];
                                                                                                                    
    rhs[0] = Vk;
    rhs[1] = U2;
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=Vk*U2 */
                                                                                                                              
    rhs[0] = U1;
    rhs[1] = lhs[0];
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=U1*Vk*U2 */
    mxDestroyArray(rhs[1]);
                
    rhs[0] = X;
    rhs[1] = lhs[0];
    mexCallMATLAB(1,lhs,2, rhs, "mtimes"); /* lhs[0]=X*U1*Vk*U2 */
    mxDestroyArray(rhs[1]);

    return lhs[0];
}
        
/******* compute V ***********/
/** Gradient of a single Givens rotation **/
mxArray* gradU(mxArray *theta, int k,
               int* ik, int* jk, int dim)
{
    
    double *p_theta = mxGetPr(theta);
    mxArray *V = mxCreateDoubleMatrix(dim, dim, mxREAL);
    double *p_V = mxGetPr(V);

    p_V[ik[k]+dim*ik[k]] = -sin(p_theta[k]);
    p_V[ik[k]+dim*jk[k]] = cos(p_theta[k]);
    p_V[jk[k]+dim*ik[k]] = -cos(p_theta[k]);
    p_V[jk[k]+dim*jk[k]] = -sin(p_theta[k]);

    return V;
}


/******* add a single Givens rotation to a previous one (calc U1*U2) ***********/
mxArray* U_add_single(mxArray *U1, 
                      mxArray *theta, 
                      int k,
                      int* ik, int* jk, int dim)
{            
    mxArray *Uab = mxDuplicateArray(U1);
    double *p_Uab = mxGetPr(Uab);
    int i,ind_ik,ind_jk;
    
    double *p_theta = mxGetPr(theta);
    double tt,u_ik;
    tt = p_theta[k];
    for( i=0; i<dim; i++ ){
        ind_ik = dim*ik[k] + i;
        ind_jk = dim*jk[k] + i;
        u_ik = p_Uab[ind_ik] * cos(tt) - p_Uab[ind_jk] * sin(tt);
        p_Uab[ind_jk] = p_Uab[ind_ik] * sin(tt) + p_Uab[ind_jk] * cos(tt);
        p_Uab[ind_ik] = u_ik;                               
    }
    return Uab;
}


/******* compute Uab ***********/
/** Givens rotation for angles a to b **/
mxArray* build_Uab(mxArray *theta, int a, int b,
               int* ik, int* jk, int dim)
{        
    mxArray *Uab = mxCreateDoubleMatrix(dim, dim, mxREAL);
    double *p_Uab = mxGetPr(Uab);
    int ind, k,i,j,ind_ik,ind_jk;
    /* set Uab to be an identity matrix */
    for( j=0; j<dim; j++ ){
        ind = dim*j + j;
        p_Uab[ind] = 1.0;
    }        
    
    if( b < a ) {
        return Uab;
    }
    
    double *p_theta = mxGetPr(theta);
    double tt,u_ik;
    for( k=a; k<=b; k++ ){
        tt = p_theta[k];
        for( i=0; i<dim; i++ ){
            ind_ik = dim*ik[k] + i;
            ind_jk = dim*jk[k] + i;
            u_ik = p_Uab[ind_ik] * cos(tt) - p_Uab[ind_jk] * sin(tt);
            p_Uab[ind_jk] = p_Uab[ind_ik] * sin(tt) + p_Uab[ind_jk] * cos(tt);
            p_Uab[ind_ik] = u_ik;
        }                       
    }
    return Uab;
}

/** Rotate vecotrs in X with Givens rotation according to angles in theta **/
mxArray* rotate_givens(mxArray *X, mxArray *theta, 
                       int* ik, int* jk, int angle_num, int dim)
{
    mxArray *G = build_Uab(theta, 0, angle_num-1,ik,jk,dim);
    mxArray *Y = matrix_mult(X,G);        
    mxDestroyArray(G);
    return Y;
}


/****** quality gradient *******************/
double evqualitygrad(mxArray *X, mxArray* theta,
                     int *ik, int *jk,
                     int angle_num, int angle_index,
                     int dim, int ndata)
{   
    /* build V,U,A */
    mxArray *V = gradU(theta,angle_index,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed gradU\n");
    
    mxArray *U1 = build_Uab(theta,0,angle_index-1,ik,jk,dim);
    mxArray *U2 = build_Uab(theta,angle_index+1,angle_num-1,ik,jk,dim);
    if( DEBUG )
        mexPrintf("Computed Uab\n");
    
    mxArray *A = buildA(X, U1, V, U2);
    double *p_A = mxGetPr(A);
    if( DEBUG )
        mexPrintf("Built A\n");
  
    /* get rid of no longer needed arrays */
    mxDestroyArray(V);
    mxDestroyArray(U1);
    mxDestroyArray(U2);
    
    /* rotate vecs according to current angles */   
    mxArray *Y = rotate_givens(X,theta,ik,jk,angle_num,dim);
    double *p_Y = mxGetPr(Y);
    if( DEBUG )
        mexPrintf("Rotated according to Givens successfully\n");
       
    /* find max of each row */
	double *max_values = (double*)mxCalloc(ndata, sizeof(double));
	int *max_index = (int*)mxCalloc(ndata, sizeof(int));
    int i,j, ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
		for( i=0; i<ndata; i++ ){ /* loop over all rows */
            if( max_values[i]*max_values[i] <= p_Y[ind]*p_Y[ind]  ){
                max_values[i] = p_Y[ind];
                max_index[i] = j;
            }
            ind++;
        }
    }
    if( DEBUG )
        mexPrintf("Found max of each row\n");           
    
    /* compute gradient */
    double dJ=0, tmp1, tmp2;
    ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
		for( i=0; i<ndata; i++ ){ /* loop over all rows */
            tmp1 = p_A[ind] * p_Y[ind] / (max_values[i]*max_values[i]);
            tmp2 = p_A[ndata*max_index[i]+i] * (p_Y[ind]*p_Y[ind]) / (max_values[i]*max_values[i]*max_values[i]);
            dJ += tmp1-tmp2;
            ind++;
        }
    }
	dJ = 2*dJ/ndata;
    if( DEBUG )
        mexPrintf("Computed gradient = %g\n",dJ);
    
    mxFree(max_values);
    mxFree(max_index);
    mxDestroyArray(Y);
    mxDestroyArray(A);
    return dJ;
}

  
/******** alignment quality ***********/
double evqual(mxArray *X,
              int *ik, int *jk,
              int dim,int ndata)
{
    double *p_X = mxGetPr(X);
    /* take the square of all entries and find max of each row */
    double *max_values = (double*)mxCalloc(ndata,sizeof(double));
    int *max_index = (int*)mxCalloc(ndata,sizeof(int));
    int i,j,ind = 0; 
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */
            if( max_values[i] <= p_X[ind]*p_X[ind]  ){                
                max_values[i] = p_X[ind]*p_X[ind];
                max_index[i] = j;
            }
            ind++;
        }
    }
    if( DEBUG )
        mexPrintf("Found max of each row\n");
   
    /* compute cost */
    double J=0;
    ind = 0;
    for( j=0; j<dim; j++ ){  /* loop over all columns */
        for( i=0; i<ndata; i++ ){ /* loop over all rows */
            J += p_X[ind]*p_X[ind]/max_values[i];
            ind++;
        }
    }
    J = J/ndata;
    if( DEBUG )
        mexPrintf("Computed quality = %g\n",J);
    
    mxFree(max_values);
    mxFree(max_index);    
    
    return J;
} 

  
/***************************   //////////////// main   */
  
  void
  mexFunction (
  int nlhs, mxArray* plhs[],
  int nrhs, const mxArray* prhs[])
{
    
    if( DEBUG )
        mexPrintf("The entry point!\n");
    
  /* Make sure at most two output arguments are expected */
    if (nlhs < 3) {
        mexErrMsgTxt("Too few output arguments.");
        return;
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
  /* Make sure input number is sufficient */
    if (nrhs < 1) {
        mexErrMsgTxt("Too few input arguments.");
        return;
    }
    if (nrhs > 2) {
        mexErrMsgTxt("Too many input arguments.");
        return;
    }
    
    if( DEBUG )
        mexPrintf("Just starting\n");
        
    
    /* get the number and length of eigenvectors dimensions */
    mxArray *X; X = (mxArray*)prhs[0];
    const size_t *idims =  mxGetDimensions(X);
    const int ndata = idims[0];
    const int dim = idims[1];
    if( DEBUG )
        mexPrintf("Got %d vectors of length %d\n",dim,ndata);
    
    /* get the number of angles */
    int angle_num;
    angle_num = (int)(dim*(dim-1)/2);     
    mxArray *theta = mxCreateDoubleMatrix(1,angle_num,mxREAL);
    if( DEBUG )
        mexPrintf("Angle number is %d\n",angle_num);
    
    /* get iteration scaling */
    int iter_scal;
    if( nrhs >= 2 )
        iter_scal = (int) mxGetScalar(prhs[1]);
    else
        iter_scal = 1;
    if( DEBUG )
        mexPrintf("Got iteration scaling %d\n",iter_scal);     
    
    /* build index mapping */
    int i,j,k,ind;
    int* ik = (int*)mxCalloc(angle_num,sizeof(int));
    int* jk = (int*)mxCalloc(angle_num,sizeof(int));
    k=0;
    for( i=0; i<dim-1; i++ ){
        for( j=i+1; j<=dim-1; j++ ){
            ik[k] = i;
            jk[k] = j;
            k++;
        }
    }
    if( DEBUG )
        mexPrintf("Built index mapping for %d angles\n",k);
    
    /* definitions */
    int max_iter = 100*iter_scal; 
    double dQ,Q,Q_old1;
    double alpha;
    const double dmp = 0.9;
    int iter,d;
    mxArray *theta_est,*Xrot;     
    
    theta_est = mxDuplicateArray(theta);
    double *p_theta = mxGetPr(theta);
    double *p_theta_est = mxGetPr(theta_est);    
    double* grad_curr = (double*)mxCalloc(angle_num,sizeof(double));
    double* grad_prev = (double*)mxCalloc(angle_num,sizeof(double));
    for (d = 0; d < angle_num; d++){
        grad_prev[d] = 0;
    }   
    
    Q = evqual(X,ik,jk,dim,ndata); /* initial quality */
    if( DEBUG )
        mexPrintf("Q = %g\n",Q);
    Q_old1 = Q;
    alpha = 0.5;   
    iter = 0;
    while( iter < max_iter ){ /* iterate to refine quality */                
        iter++;
        for( d = 0; d < angle_num; d++ ){
            p_theta_est[d] = p_theta[d] - dmp*grad_prev[d];    
            dQ = evqualitygrad(X,theta_est,ik,jk,angle_num,d,dim,ndata);  
            grad_curr[d] = dmp*grad_prev[d] + alpha*dQ;
        }
        for( d = 0; d < angle_num; d++ ){
            p_theta[d] = p_theta[d] - grad_curr[d];
            grad_prev[d] = grad_curr[d];
        }
        Xrot = rotate_givens(X,theta,ik,jk,angle_num,dim);
        Q = evqual(Xrot,ik,jk,dim,ndata); 
        mxDestroyArray(Xrot);   
        if (fabs(Q_old1-Q)<1e-3) 
            break;
        else{           
            Q_old1 = Q;
        }              
    }       
    
    if( DEBUG )
        mexPrintf("Done after %d iterations, Quality is %g\n",iter,Q);
    
    Xrot = rotate_givens(X,theta,ik,jk,angle_num,dim);
    
    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(ndata, dim, mxREAL);
    double *p_Xout = mxGetPr(plhs[0]);
    double *p_Xrot = mxGetPr(Xrot);
    memcpy(p_Xout, p_Xrot, ndata*dim*sizeof(double));       
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *p_Q = mxGetPr(plhs[1]);
    *p_Q = Q;
    plhs[2] = mxCreateDoubleMatrix(1, angle_num, mxREAL);
    double *p_ang = mxGetPr(plhs[2]);    
    memcpy(p_ang, p_theta, angle_num*sizeof(double));       
        
    /* free allocated memory */
    mxFree(ik);
    mxFree(jk);   
    mxDestroyArray(theta_est);
    mxFree(grad_curr);
    mxFree(grad_prev);    

    if( DEBUG )
        mexPrintf("Done evrot\n");
        
    return;
}


