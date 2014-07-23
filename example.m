%% Compile mex code for new kernel
clear all
delete('./*.o');
delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
syms r;                           % distance of two points (radius)
kernel = 1 ./ r.^2;               % radius basis function 
% kernel = exp(-r^2);
outputfile = 'mexFMM3D';
homogen = -2;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
                                 % anti-symmetric: 0
make(r,kernel,homogen,symmetry,outputfile);
 

%% Example: Q*H

% Info on dimensions
Ns  = 10000;    % Number of sources in simulation cell
Nf  = 10000;    % Number of fields in simulation cell
m   = 2;        % number of columns of H
L = 1.0;        % Length of simulation cell (assumed to be a cube)
% 3-D locations, stored column-wise i.e. (x | y | z)

source = (rand(Ns,3) - 0.5) .* L;
% field = (rand(Nf,3) - 0.5) .* L; 
field = source;

H = rand(Ns,m); 

nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 3;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation


%%
% Compute matrix-vectors product QH
% QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
[QH,QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);


 
 


