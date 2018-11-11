function [x, niters] = cgsolve2(A, b)
% cgsolve : Solve Ax=b by conjugate gradients
%
% [x, niters] = cgsolve(A, b);
%
% Given symmetric positive definite sparse matrix A and vector b, 
% this runs conjugate gradient to solve for x in A*x=b.
% It iterates until the residual norm is reduced by 10^-6,
% or for at most max(100,sqrt(n)) iterations.
%
% The code below follows the slides on the CS140 web site.
%
% John R. Gilbert    12 Feb 2006

n = length(b);
maxiters = max(100,sqrt(n));
normb = norm(b);
x = zeros(n,1);
r = b;
rtr = r'*r;
d = r;
niters = 0;
while sqrt(rtr)/normb > 2e-6  &&  niters < maxiters
    niters = niters+1;
    Ad = A(d);
    alpha = rtr / (d'*Ad);
    x = x + alpha * d;
    r = r - alpha * Ad;
    rtrold = rtr;
    rtr = r'*r;
    beta = rtr / rtrold;
    d = r + beta * d;
end;