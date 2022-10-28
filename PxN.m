function [PN]=PxN(P,N)
%
% Given two compatible polynomial matrices P(s) and N(s)
% compute their polynomial product P*N via convolution
%
m=size(P,1);r=size(P,2);n=size(N,2);dp=size(P,3)-1;dn=size(N,3)-1;
imax=dp+dn+1;PN=zeros(m,n,imax);
for i=1:imax,
    j1=max([1,i-dn]); % first index of P involved in product
    j2=min([dp+1,i]); % last index of P involved in product
    for isum=j1:j2,
        PN(:,:,i)=PN(:,:,i)+P(:,:,isum)*N(:,:,i-isum+1);
        end
    end
end