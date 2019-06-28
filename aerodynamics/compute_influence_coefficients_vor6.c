/*Authors: Klaus Seywald (klaus.seywald@mytum.de) 
         and Simon Binder (simon.binder@tum.de)

This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
*/
#include "mex.h"
#include <math.h>
/*
 * xtimesy.c - example found in API guide
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 */

/* $Revision: 1.10.6.4 $ */

void compute_influence_coefficients_vor5(double *grid,double *panels, double *colloc,double *te_idx,double *colloc_nvec,double *w,double *w_x,double *w_y,double *w_z,mwSize n_colloc,mwSize n_pan,double * Uinf)
{
    mwSize i,j,count=0;
    
    double nvel;
    double wAB_1,wAB_2,wAB_3;
    double r1xr0_1,r1xr0_2,r1xr0_3;
    double r1ixr0_1,r1ixr0_2,r1ixr0_3;
    double norm_r0;
    double r01,r02,r03,r11,r12,r13,r21,r22,r23,r1xr2_1,r1xr2_2,r1xr2_3,r1dr2_1,r1dr2_2,r1dr2_3,dotp_r0_r1dr2,norm_r1xr2,normr1,normr2;
    double p_i_1,p_i_2,p_i_3,p_o_1,p_o_2,p_o_3;
    double coredist;
    double coredisti;
    double vortex[18];
    int k;
    for (k=0;k<n_pan;k++){
        
        p_i_1=*(grid+((int)*(panels+4*k)-1)*3)*0.75+*(grid+((int)*(panels+4*k+3)-1)*3)*0.25;
        p_i_2=*(grid+((int)*(panels+4*k)-1)*3+1)*0.75+*(grid+((int)*(panels+4*k+3)-1)*3+1)*0.25;
        p_i_3=*(grid+((int)*(panels+4*k)-1)*3+2)*0.75+*(grid+((int)*(panels+4*k+3)-1)*3+2)*0.25;
       
        p_o_1=*(grid+((int)*(panels+4*k+1)-1)*3)*0.75+*(grid+((int)*(panels+4*k+2)-1)*3)*0.25;
        p_o_2=*(grid+((int)*(panels+4*k+1)-1)*3+1)*0.75+*(grid+((int)*(panels+4*k+2)-1)*3+1)*0.25;
        
        p_o_3=*(grid+((int)*(panels+4*k+1)-1)*3+2)*0.75+*(grid+((int)*(panels+4*k+2)-1)*3+2)*0.25;   
        vortex[0]=*(Uinf)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3);
        vortex[1]=*(Uinf+1)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3+1);
        vortex[2]=*(Uinf+2)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3+2);
        vortex[3]=*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3);
        vortex[4]=*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3+1);
        vortex[5]=*(grid+((int)*(te_idx+((int)*(panels+4*k+3)-1)-1)-1)*3+2);
        vortex[6]=p_i_1;
        vortex[7]=p_i_2;
        vortex[8]=p_i_3;
        vortex[9]=p_o_1;
        vortex[10]=p_o_2;
        vortex[11]=p_o_3;
        vortex[12]=*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3);
        vortex[13]=*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3+1);
        vortex[14]=*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3+2);
        vortex[15]=*(Uinf)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3);
        vortex[16]=*(Uinf+1)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3+1);
        vortex[17]=*(Uinf+2)*1E10+*(grid+((int)*(te_idx+((int)*(panels+4*k+2)-1)-1)-1)*3+2);

           
     
      //  mexPrintf("hallooo2   %d, %d \n",((int)*(panels+4*k+3)-1),(((int)*(panels+4*k+2)-1)-1)-1);
        for (j=0;j<n_colloc;j++){
            //Berechnung der 1/4 Punkte 
            //print 1/4 punkt von sowie colloc Punkt Panel 5 
//                   if ((k==0) && (j==0))
//                 {
//                     
//                   mexPrintf("j k   %d %d\n",j,k);
//                   mexPrintf("colloc_1   %f \n",(double)*(colloc+j*3));
//                   mexPrintf("colloc_2   %f \n",(double)*(colloc+1+j*3));
//                   mexPrintf("colloc_3   %f \n",(double)*(colloc+2+j*3));
//                   mexPrintf("p_m_1   %f \n",p_m_1);
//                   mexPrintf("p_m_2   %f \n",p_m_2);
//                   mexPrintf("p_m_3   %f \n",p_m_3);
//                 }   
            
            wAB_1=0;
            wAB_2=0;
            wAB_3=0;
            
            for (i=0;i<5;i++){
                r11=*(colloc+j*3)-*(vortex+i*3);
                r12=*(colloc+1+j*3)-*(vortex+i*3+1);
                r13=*(colloc+2+j*3)-*(vortex+i*3+2);
                
                r21=*(colloc+j*3)-*(vortex+(i+1)*3);
                r22=*(colloc+1+j*3)-*(vortex+(i+1)*3+1);
                r23=*(colloc+2+j*3)-*(vortex+(i+1)*3+2);
                
                
                r01=*(vortex+(i+1)*3)-*(vortex+i*3);
                r02=*(vortex+(i+1)*3+1)-*(vortex+i*3+1);
                r03=*(vortex+(i+1)*3+2)-*(vortex+i*3+2);
                
                r1xr2_1=r12*r23-r13*r22;
                r1xr2_2=r13*r21-r11*r23;
                r1xr2_3=r11*r22-r12*r21;
                
                
                norm_r1xr2=sqrt(r1xr2_1*r1xr2_1+r1xr2_2*r1xr2_2+r1xr2_3*r1xr2_3);
                norm_r0=sqrt(r01*r01+r02*r02+r03*r03);
                
                normr1=sqrt(r11*r11+r12*r12+r13*r13);
                normr2=sqrt(r21*r21+r22*r22+r23*r23);
                
                r1xr0_1=r12*r03-r13*r02;
                r1xr0_2=r13*r01-r11*r03;
                r1xr0_3=r11*r02-r12*r01;
                
                
                coredist=sqrt(r1xr0_1*r1xr0_1+r1xr0_2*r1xr0_2+r1xr0_3*r1xr0_3)/norm_r0;
                 /*if  ( (!((normr1<=1E-5)||(normr2<=1E-5))) && coredist>1E-8)*/ 
                if (i<4 && i>0) 
                {
                    if  ( (!((normr1<=1E-7)||(normr2<=1E-7))) && coredist>1E-7)
                    {
                        r1dr2_1=r11/normr1-r21/normr2;
                        r1dr2_2=r12/normr1-r22/normr2;
                        r1dr2_3=r13/normr1-r23/normr2;

                        dotp_r0_r1dr2=r1dr2_1*r01+r1dr2_2*r02+r1dr2_3*r03;

                        wAB_1=wAB_1-r1xr2_1/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;
                        wAB_2=wAB_2-r1xr2_2/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;
                        wAB_3=wAB_3-r1xr2_3/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;



                    }
                }
                else //Biot_savart law modified for the semi-infinite vortex filaments to remove the singularities during wake vortex interference with downstream surfaces.
                {
                    if  ( (!((normr1<=1E-7)||(normr2<=1E-7))) && coredist>1E-7)
                    {
                        r1dr2_1=r11/normr1-r21/normr2;
                        r1dr2_2=r12/normr1-r22/normr2;
                        r1dr2_3=r13/normr1-r23/normr2;

                        dotp_r0_r1dr2=r1dr2_1*r01+r1dr2_2*r02+r1dr2_3*r03;

                        wAB_1=wAB_1-r1xr2_1/(norm_r1xr2*norm_r1xr2+norm_r0*norm_r0)*dotp_r0_r1dr2;
                        wAB_2=wAB_2-r1xr2_2/(norm_r1xr2*norm_r1xr2+norm_r0*norm_r0)*dotp_r0_r1dr2;
                        wAB_3=wAB_3-r1xr2_3/(norm_r1xr2*norm_r1xr2+norm_r0*norm_r0)*dotp_r0_r1dr2;               
 
                    }
                   
                }
            }
         // *(w+n_colloc*k+j)=*(colloc_nvec+j*3)*wAB_1+*(colloc_nvec+1+j*3)*wAB_2+*(colloc_nvec+2+j*3)*wAB_3;
            *(w_x+n_colloc*k+j)=wAB_1;
            *(w_y+n_colloc*k+j)=wAB_2;
            *(w_z+n_colloc*k+j)=wAB_3;
//   *(wind+n_colloc*k+j)=*(colloc_nvec+j*3)*wAB_1i+*(colloc_nvec+1+j*3)*wAB_2i+*(colloc_nvec+2+j*3)*wAB_3i;
        }
    }
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *grid,*colloc,*colloc_nvec,*w,*w_x,*w_y,*w_z,*Uinf;
    double *panels,*te_idx;
    size_t n_colloc;
    size_t n_pan;
    
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
    if(nrhs!=6)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
                "Two inputs required.");
    if(nlhs!=3)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
                "One output required.");
    
    /* check to make sure the first input argument is a scalar */
    //if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    //    mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    //  mexErrMsgIdAndTxt( "MATLAB:xtimesy:xNotScalar",
    //          "Input x must be a scalar.");
    //}
    
    /*  get the scalar input vortex */
    grid =(double *)mxGetPr(prhs[0]);
    panels =(double *)mxGetPr(prhs[1]);
    te_idx =(double *)mxGetPr(prhs[2]);
    /*  create a pointer to the input matrix colloc */
    
    colloc =(double *)mxGetPr(prhs[3]);
    n_colloc= mxGetN(prhs[3]);
    n_pan= mxGetN(prhs[1]);
    
    
    colloc_nvec =(double *)mxGetPr(prhs[4]);
    Uinf =(double *)mxGetPr(prhs[5]);
    /*  get the dimensions of the matrix input y */
    /* mrows = mxGetM(prhs[0]);
 /* ncols = mxGetN(prhs[0]);
  
  /*  set the output pointer to the output matrix */
    plhs[0]=mxCreateDoubleMatrix(n_colloc,n_pan,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_colloc,n_pan,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_colloc,n_pan,mxREAL);
    /*  create a C pointer to a copy of the output matrix */
    
    w_x = mxGetPr(plhs[0]);
    w_y = mxGetPr(plhs[1]);
    w_z = mxGetPr(plhs[2]);
    
//     mexPrintf("bla %i \n ",n_pan);
    // mexPrintf("Vortex %f,%f,%f \n %f,%f,%f \n %f,%f,%f \n ", *(vortex),*(vortex+1),*(vortex+2),*(vortex+3),*(vortex+4),*(vortex+5),*(vortex+6),*(vortex+7),*(vortex+8));
    /*  call the C subroutine */
    compute_influence_coefficients_vor5(grid,panels,colloc,te_idx,colloc_nvec,w,w_x,w_y,w_z,n_colloc,n_pan,Uinf);
    
}
