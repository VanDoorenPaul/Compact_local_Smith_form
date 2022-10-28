% This generates a mxn polynomial matrix with Smith form
% and size mxn with r structural indices at 0 (r=normal rank)
d=3;m=4;n=5;s=[0,1,3];r=length(s);
La=zeros(m,n,d+1);
for i=1:r,
    La(i,i,s(i)+1)=1;
end
% Now construct random polynomial matrices M and N of
% degree k and form P= M.La.N with local Smith form La.
k=2;
for iruns=1:10,
M=randn(m,m,k+1).^iruns;
N=randn(n,n,k+1).^iruns;
P=PxN(PxN(M,La),N);
% Its degree will be 2k+d
% We also set tol 
normP(iruns)=norm(P(:));
tol=10*eps*normP(iruns);
% now run the compact Smith algorithm
[Pc,N,L]=compactSmith(P,r,tol);
% residual errors
Res=PxN(P,N)-PxN(Pc,L);Res=Res(:,1:r,:);resP(iruns)=norm(Res(:));
normN(iruns)=norm(N(:));
end
[order,Ind]=sort(normP);
normP=normP(Ind);normN=normN(Ind);resP=resP(Ind);
format short e
Table=[normP(:) resP(:)./normP(:) normN(:)]
