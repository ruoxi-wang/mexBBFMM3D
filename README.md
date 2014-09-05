#mexBBFMM3D  



###1. INTRODUCTION

mexBBFMM3D is MATLAB interface for an open source package (BBFMM3D) of the <a href="http://www.sciencedirect.com/science/ article/pii/S0021999109004665">Black-box Fast Multipole Method</a> in 3 dimensions.   
The Black-box Fast Multipole Method is an O(N) fast multipole method, which is a technique to calculate sums of the form  

![](http://latex.codecogs.com/gif.latex?f%28x_i%29%20%3D%20%5Cdisplaystyle%20%5Csum_%7Bj%3D1%7D%5EN%20K%28x_i%2Cy_j%29%20%5Csigma_j%2C%20%5C%2C%5C%2C%5C%2C%20%5Cforall%20i%20%5Cin%5C%7B1%2C2%2C%5Cldots%2CN%5C%7D)

where ![](http://latex.codecogs.com/gif.latex?K%28x_i%2Cx_j%29) is kernel function, ![](http://latex.codecogs.com/gif.latex?x_i) are observation points, ![](http://latex.codecogs.com/gif.latex?y_j) are locations of sources, and ![](http://latex.codecogs.com/gif.latex?%5Csigma_i) are charges at corresponding locations.
BBFMM3D provides an O(N) solution to matrix-vector products of the type Ax. In that case the relation between A and K is:
![](http://latex.codecogs.com/gif.latex?A_%7Bij%7D%20%3D%20K%28x_i%2Cy_j%29)



This implementation of the FMM differs from other methods by the fact that it is applicable to all smooth kernels K. [Give examples of RBF kernels, 1/r, log r, Stokes, etc.].

The approximation scheme used in the FMM relies on Chebyshev interplation to construct low-rank approximations for well-separated clusters. In addition the use of Singular Value Decomposition ensures that the computational cost is minimal. In particular the rank is optimally chosen for a given error. 

Please cite the following paper if you use this code:

Fong, William, and Eric Darve. "The black-box fast multipole methodshod." Journal of Computational Physics 228, no. 23 (2009): 8712-8725. You can see details <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">here</a>.

###2. DIRECTORIES AND FILES


	./example.m	:	:	Example of how to use mexBBFMM3D 
	./make.m		:	Makefile 
	./include/		:	Relevant header files  
	./mexFMM3D.cpp	:	mex function  
	./eigen/		:	Eigen Library  
	./BBFMM3D/		: 	BBFMM3D library
	./README.md		:	This file  
	
###3. TUTORIAL
####3.1 To Get Started  
To check whether things are set up correctly, you can go to the directory where this README.m is and run example.m.  Going through example.m should be self-explanatory, but more details are provided below. 
####3.2 Basic usage

#####3.2.1 Complie

To use mexFMM3D, you need to compile first.  

	syms r;                          
	kernel = 1 ./ r.^2;                
	outputfile = 'mexFMM3D';
	homogen = -2;                    
	symmetry = 1;                     
                                 	 
	make(r,kernel,homogen,symmetry,outputfile);
where 

*  r: 		
	the distance of two points (radius)
*  kernel:	
	 the kernel function (radius basis function)
*  outputfile: 	
	 the name of the routine you need to run
*  homogen:		
	The homogeneous property of kernel.(The cost and memory requirements of the pre-computation step can be reduced for the case of homogeneous kernels.)
	* homogen = 0: if the kernel funtion is not homogeneous.
	* homogen = m: if the kernel funtion is homogeneous of degree m, i.e. <img src="http://latex.codecogs.com/gif.latex? $K(\alpha x, \alpha y) = \alpha^m K(x,y)$." border="0"/>
* symmetry:   
	The symmetry property of the kernel.
	* symmetry = 0:  no symmetric property; 
	* symmetry = 1: symmetric kernel; 		[K(x,y) = K(y,x)]
	* symmetry = -1: anti-symmetric kernel.[K(x,y) = - K(y,x)]

#####3.2.2 Set up problems and parameters for BBFMM3D
This is a sample problem, make sure your points stay in the cell of length L.
	
	% Info on dimensions
	Ns  = 10000;    % Number of sources in simulation cell
	Nf  = 10000;    % Number of fields in simulation cell
	m   = 2;        % number of columns of H

	% 3-D locations, stored column-wise i.e. (x | y | z)
	source = rand(Ns,3) - 0.5;
	% field = rand(Nf,3) - 0.5; 
	field = source;

	H = rand(Ns,m); 
	
	L = 1.0;            % Length of simulation cell (assumed to be a cube)
	nCheb = 4;          % Number of Chebyshev nodes per dimension
	level = 3;          % Level of FMM tree
	use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation
	
	
	

* L:   
	Length of simulation cell (assumed to be a cube).
* level:  
	The number of levels in the hierarchy tree
* nCheb:  
	Number of Chebyshev nodes per dimension ( the larger nCheb is, the better accuracy it will achive)
* use_chebyshev(int):  
	Label to indicate whether using chebyshev interpolation formula or uniform interpolation formula.  
	use_chebyshev = 1: chebyshev interplation formula (recommended if there is no memory concern)
	use_chebyshev = 0: uniform interpolation formula.



#####3.2.3 Compute matrix-vector(s) multiplication

	QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
	
If you want to compare it with exact computation:  
	
	[QH, QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);


####3.3 Pre-computation

The power of this package is in the pre-computing part, which is much more computationally expensive than computing part. This package takes advantage of the fact that for a given kernel and number of chebyshev nodes, the precomputing part is the same, so for a fixed kernel and number of chebyshev nodes, it generates 3 files storing information of FMM tree in the folder *BBFMM3D/output/*. Everytime when we use the same kernel type and number of chebyshev nodes, we can directly read from the files, which would save a lot of time.

Note: it is true that sometimes with the pre-computation step, the code will be slower than direct calculation. But if the file already exists, then when doing more computations it will be faster than direct calculation. If you change the kernel, make sure to delete the existed file in *BBFMM3D/output/* before you run the code.

###3.4

The C++ version of BBFMM3D has more options, and is faster than using this MATLAB interface.






