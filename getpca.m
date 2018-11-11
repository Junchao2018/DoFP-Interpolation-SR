function [P, mx]=getpca(X)

%X: MxN matrix (M dimensions, N trials)
%Y: Y=P*X
%P: the transform matrix
%V: the variance vector
[M,N]=size(X);
mx   =  mean(X,2);
mx2  =  repmat(mx,1,N);
SubX=X-mx2;
CovX=SubX*SubX'/(N-1);
ind = CovX<0;
CovX(ind) =0.0001;

[P,V]=eig(CovX);

V=diag(V);
[t,ind]=sort(-V);
% V=V(ind);
P=P(:,ind);
P=P';
return;

