function [Pc,N,L]=compactSmith(P,r,tol)
%
% The function [Pc,N,L]=compactSmith(P,r,tol) computes
% the local compact Smith form  P*N=Pc*L at the zero s = 0  
% of a mxn polynomial matrix P(s) = sum P_i s^i of normal rank r.
%
% Here N(s) is a nxr submatrix of a unimodular matrix,
%      Pc(s) is mxr and Pc(0) has linearly independent columns, 
%      L(s) is a rxr diagonal matrix with the structural monomials,
%      and is the local Smith form of P(s)
%
% The factorization uses a tolerance tol for all rank calculations, 
% which should be of the order of eps*norm(P(:))
% The mxn polynomial matrix P of degree d should be given in a 
% 3D array of dimensions m x n x (d+1) 
%
% This function uses the Matlab functions CGS, PxN and Trim
%
m=size(P,1);n=size(P,2);d=size(P,3)-1; % dimensions of P
% Initialization Unimodular factor N with Q from svd of P0
rho=0; % this is in fact rho(-1);
[U,S,Q]=svd(P(:,:,1));
N(:,:,1)=Q;
rhoplus=rank(S,tol);
i=1; e(i)=rhoplus;
L=eye(n,n);
% Multiply the polynomial matrix with Q
P=PxN(P,Q);
% Now perform the subsequent reductions and update factor N
while rhoplus < r,
rho=rhoplus;
% Update L(s) = Lambda(s)
L(rho+1:n,rho+1:n,i+1)=L(rho+1:n,rho+1:n,i);
L(rho+1:n,rho+1:n,i)=zeros(n-rho,n-rho);
% Reduction step
[P,rhoplus,Z,Q]=CGS(P,rho,tol);
% Updating N 
if rho>0,  % construct unimodular Zs 
   Zs=eye(n,n);
   r1=0;
   for j=1:i, 
      r0=r1+1;r1=r1+e(j);
      Zs(r0:r1,rho+1:n,i-j+2)=Z(r0:r1,:);
   end
   N=Trim(PxN(N,Zs),tol);
end
% apply Q to N 
N(:,rho+1:n,:)=PxN(N(:,rho+1:n,:),Q);
% Update also index set
i=i+1;
e(i)=rhoplus-rho;
end
% Now truncate the outputs to its leading r rows/columns
% and trim Pc.
Pc=Trim(P(:,1:r,:),tol);
N=N(:,1:r,:);
L=L(1:r,1:r,:);