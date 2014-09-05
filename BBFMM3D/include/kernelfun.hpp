class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int level, int n,  double epsilon, int
       use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType,doft*dof) {
       homogen = 0;
       symmetry = 1;
       dof->f = 1;
       dof->s = 1;
       kernelType = "myKernel";}
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                               double *K, doft *dof) {
        double r =  sqrt((sourcepos.x - fieldpos.x)*(sourcepos.x - fieldpos.x) + (sourcepos.y - fieldpos.y)*(sourcepos.y - fieldpos.y) + (sourcepos.z - fieldpos.z)*(sourcepos.z - fieldpos.z));
        double t0;         //implement your own kernel on the next line
         t0 = exp((r*r)*-1.0E4);
       *K =  t0;
    }
};
