
#include "bbfmm3d.hpp"

extern void _main();
#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)





void read_data(const mxArray* source_IN, const mxArray* field_IN, int& Ns, int& Nf, vector3* source, vector3* field){
    
    double* source_pr = mxGetPr(source_IN);
    double* field_pr = mxGetPr(field_IN);
    for (int i = 0; i < Ns; i++) {
        source[i].x = (double)source_pr[i+Ns*0];
        source[i].y = (double)source_pr[i+Ns*1];
        source[i].z = (double)source_pr[i+Ns*2];
    }
    for (int i = 0; i < Nf; i++) {
        field[i].x = (double)field_pr[i+Nf*0];
        field[i].y = (double)field_pr[i+Nf*1];
        field[i].z = (double)field_pr[i+Nf*2];
    }
   
}


void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])  {
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    #define QH_OUT             plhs[0]
    #define QHexact_OUT        plhs[1]
    #define source_IN          prhs[0]
    #define field_IN           prhs[1]
    #define charge_IN          prhs[2]
    #define nCheb_IN           prhs[3]
    #define level_IN           prhs[4]
    #define L_IN               prhs[5]
    #define use_cheby_IN       prhs[6]
    #define singular_IN        prhs[7]  
   
          

    
    // Check number of argument
    if(nrhs != 7) {
        mexErrMsgTxt("Wrong number of input arguments");
    }else if(nlhs > 2){
        mexErrMsgTxt("Too many output arguments");
    }
    
    if( !IS_REAL_2D_FULL_DOUBLE(source_IN)) {
        mexErrMsgTxt("Third input argument is not a real 2D full double array.");
    }
    if( !IS_REAL_2D_FULL_DOUBLE(field_IN)) {
        mexErrMsgTxt("Third input argument is not a real 2D full double array.");
    }
    if( !IS_REAL_2D_FULL_DOUBLE(charge_IN)) {
        mexErrMsgTxt("Third input argument is not a real 2D full double array.");
    }
    if( !IS_REAL_SCALAR(nCheb_IN)){
        mexErrMsgTxt("nChebnotes must be a real double scalar");
    }
    if( !IS_REAL_SCALAR(level_IN)){
        mexErrMsgTxt("level must be a real double scalar");
    }
    if( !IS_REAL_SCALAR(L_IN)){
        mexErrMsgTxt("L must be a real double scalar");
    }
    if( !IS_REAL_SCALAR(use_cheby_IN)){
        mexErrMsgTxt("use_cheby must be a real double scalar");
    }



    double L = *mxGetPr(L_IN);         // Length of simulation cell (assumed to be a cube)
    int n = *mxGetPr(nCheb_IN);        // Number of Chebyshev nodes per dimension
    doft dof;
    dof.f   = 1;
    dof.s   = 1;
    int Ns = mxGetM(source_IN);        // Number of sources in simulation cell
    int Nf = mxGetM(field_IN);         // Number of field points in simulation cell
    int m  = mxGetN(charge_IN);        // Number of r.h.s.
    int level = *mxGetPr(level_IN);
    double eps = 1e-9;
    int use_chebyshev = *mxGetPr(use_cheby_IN);
        
    vector3* source = new vector3[Ns];    // Position array for the source points
    vector3* field = new vector3[Nf];     // Position array for the field points

    read_data(source_IN, field_IN, Ns, Nf, source, field);
    
    double * charges;
    charges = mxGetPr(charge_IN);

    double err;
    QH_OUT = mxCreateDoubleMatrix(Nf*dof.f, m, mxREAL);
    double* stress = mxGetPr(QH_OUT);

    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    mexPrintf("\nStarting FMM computation...\n");

    /*****      Pre Computation     ******/
    clock_t  t0 = clock();
    myKernel Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
    clock_t t1 = clock();
    double tPre = t1 - t0;
    
    /*****      FMM Computation     *******/
    t0 = clock();
    H2_3D_Compute<myKernel> compute(&Atree, field, source, Ns, Nf, charges,m, stress);
    t1 = clock();
    double tFMM = t1 - t0;
    
    

    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/
    
    mexPrintf("\nPre-computation time: %.4f\n",double(tPre) / double(CLOCKS_PER_SEC));
    mexPrintf("FMM computing time: %.4f\n",double(tFMM) / double(CLOCKS_PER_SEC));
    mexPrintf("FMM total time: %.4f\n",double(tFMM+tPre) / double(CLOCKS_PER_SEC));
    if (nlhs == 2) {
        mexPrintf("\nStarting direct computation...\n");
        t0 = clock();
        QHexact_OUT = mxCreateDoubleMatrix(Nf*dof.f, m, mxREAL);
        double* stress_dir = mxGetPr(QHexact_OUT);

        DirectCalc3D(&Atree, field, Nf, source, charges, m, Ns, 0 , L, stress_dir);
        t1 = clock();
        double tExact = t1 - t0;
        mexPrintf("Exact computing time: %.4f\n",double(tExact) / double(CLOCKS_PER_SEC));
        // Compute the 2-norm error
        err = ComputeError(stress,stress_dir,Nf,&dof,m);
        mexPrintf("Relative Error: %e\n",err);
    }

    /*******            Clean Up        *******/
    delete []source;
    delete []field;
    return;
}
