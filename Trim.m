function [P]=Trim(P,tol)
%
% Given a mxn polynomial matrix P of degree d, the function
%
%             [P]=Trim(P,tol)
%
% eliminates the leading coefficients P(:,:,i)
% that are smaller than tol in norm for i>0 and
% thus reduces the degree of P.
%
dp1=size(P,3);
while dp1 > 0,
    if norm(P(:,:,dp1),'fro') <= tol, 
        dp1=dp1-1;
    else
        break,
    end
end
P=P(:,:,1:dp1);