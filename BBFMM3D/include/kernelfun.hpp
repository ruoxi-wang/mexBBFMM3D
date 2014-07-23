class myKernel: public H2_3D_Tree {
public:
    myKernel(doft* dof, double L, int level, int n,  double epsilon, int
       use_chebyshev):H2_3D_Tree(dof,L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType) {
       homogen = -2;
       symmetry = 1;
       kernelType = "myKernel";}
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                               double *K, doft *dof) {
        double r =  sqrt((sourcepos.x - fieldpos.x)*(sourcepos.x - fieldpos.x) + (sourcepos.y - fieldpos.y)*(sourcepos.y - fieldpos.y) + (sourcepos.z - fieldpos.z)*(sourcepos.z - fieldpos.z));
        double t0;         //implement your own kernel on the next line
         t0 = 1.0/(r*r);
       *K =  t0;
    }
};
