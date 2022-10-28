README file for CompactLocalSmith 

This directory contains Matlab files to compute the
 
Compact Local Smith form of a mxn polynomial matrix of known normal rank r

    P(s) = P_0 + P_1.s + ... + P_d.s^d

where P(s) should be given in a 3D array of dimensions m x n x (d+1).

The computed compact factorization is P*N=Pc*L at the zero s = 0 
 
where N(s) is a nxr submatrix of a unimodular matrix,
Pc(s) is mxr and Pc(0) has linearly independent columns,
L(s) is a square rxr diagonal matrix with the structural monomials,
and is the compact local Smith form of P(s)

The factorization uses a tolerance tol for all rank calculations, 
which should be of the order of eps*norm(P(:))

This function uses the Matlab functions CGS, PxN and Trim

All files are Matlab function files except for the Test Driver 
which is a Matlab script that produces test results for the codes.

Each Matlab function file is commented and has a built-in help.