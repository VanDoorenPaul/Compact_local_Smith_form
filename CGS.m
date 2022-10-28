function [P,rhoplus,Z,V]=CGS(P,rho,tol)
% Given a polynomial matrix P(s) stored in a 3D array
% 
% P0 = P(:,:,1), ... , Pd= having P(:,:,d+1)
% 
% and with constant coefficient P0=[X, Y]
% where X has rho linearly independent columns
% This routine performs one Toeplitz rank search step :
% 1. orthogonalize Y to X (Gramm-Schmidt ortho) using
%    Y:=Y-XZ where Z=X\Y and apply [I -Z; 0 I] to P(s)
% 2. then compress the columns of the resulting matrix Y
% 3. and return the rank increment rhoplus (using tol)
%
n=size(P,2);m=size(P,1);dp1=size(P,3);d=dp1-1;
% These are the dimensions of P(s)
% Start by shifting down the n-rho last columns
for i=1:d, 
        P(:,rho+1:n,i)=P(:,rho+1:n,i+1);
end
P(:,rho+1:n,dp1)=zeros(m,n-rho);
% Now proceed with the Gram Schmidt elimination
if rho > 0, 
    X=P(:,1:rho,1);
    Y=P(:,rho+1:n,1);
    % orthogonalize with [I Z; 0 I]
    Z=-X\Y;
    for i=1:dp1, 
        P(:,rho+1:n,i)=P(:,rho+1:n,i)+P(:,1:rho,i)*Z;
    end
else
    Z=zeros(0,n);
end
% now compress the columns of Y
[U,S,V]=svd(P(:,rho+1:n,1));
for i=1:dp1, 
    P(:,rho+1:n,i)=P(:,rho+1:n,i)*V;
end
% This is the new rho for the next step :
rhoplus=rank(S,tol)+rho;

