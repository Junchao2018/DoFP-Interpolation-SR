% ========================================================================
% Sparse representation for DoFP image interpolation, Version 1.0
% Copyright(c) 2018  Junchao Zhang, Haibo Luo, Rongguang Liang, 
% Ashfaq Ahmed,Xiangyue Zhang, Bin Hui and Zheng Chang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%
%----------------------------------------------------------------------
% This is an demo of "Sparse representation-based demosaicing
% method for microgrid polarimeter imagery"
% 
% Please cite the following paper when you use it:
%
% Junchao Zhang, Haibo Luo, Rongguang Liang, Ashfaq Ahmed, Xiangyue Zhang, 
% Bin Hui, and Zheng Chang, "Sparse representation-based demosaicing method 
% for microgrid polarimeter imagery," Optics Letters 43(14), 3265-3268 (2018).
%----------------------------------------------------------------------
function  [im1, im_out, im_wei] = SP_Interp( im, Dict, par, Tau, blk_arr, wei_arr)
opts = par.s1;
[h, w, ch] = size(im);

b = opts.win;
b2 = b*b*ch;
k = 0;
s = opts.step;

N =  h-b+1;
M =  w-b+1;
L =  N*M;
r =  [1:s:N];
r =  [r r(end)+1:N];
c =  [1:s:M];
c =  [c c(end)+1:M];
X =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
         for tt =1:ch
            blk  =  im(i:h-b+i,j:w-b+j,tt);
            blk  =  blk(:);
            X(k+(tt-1)*b*b,:) =  blk';
        end
    end
end
X_m = zeros(length(r)*length(c),b2,'single');
X0 = X';
for i = 1:opts.nblk
    v = wei_arr(:,i);
    X_m = X_m(:,:) + X0(blk_arr(:,i),:) .*v(:, ones(1,b2));
end
X_m = X_m';
ind = zeros(N,M);
ind(r,c) = 1;
X = X(:, ind~=0);
N = length(r);
M = length(c);
L = N*M;
Y = zeros( b2, L );

for i  = 1:L   
    P = reshape(Dict(:, i), b*b*ch, b*b*ch);
    tau = Tau(:, i);      
    XX = X(:, i);
    XXm = X_m(:,i);
    tmp = (P*XX + 2*tau.*(P*XXm))./(1+2*tau);
    Y(:, i) = P'*tmp;
end

im_out   =  zeros(h,w,ch);
im_wei   =  zeros(h,w,ch);
k        =  0;
for tt = 1:ch
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            im_out(r-1+i,c-1+j,tt)  =  im_out(r-1+i,c-1+j,tt) + reshape( Y(k,:)', [N M]);
            im_wei(r-1+i,c-1+j,tt)  =  im_wei(r-1+i,c-1+j,tt) + 1;
        end
    end
end
im1  =  im_out./(im_wei+eps);
return;
