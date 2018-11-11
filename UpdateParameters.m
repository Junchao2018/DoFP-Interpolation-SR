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
function  [Tau]   =  UpdateParameters( im, Dict, par, blk_arr, wei_arr)
b        =  par.s1.win;
s        =  par.s1.step;
b2       =  b*b;
[h  w ch]   =  size(im);

N       =  h-b+1;
M       =  w-b+1;
r       =  [1:s:N];
r       =  [r r(end)+1:N];
c       =  [1:s:M];
c       =  [c c(end)+1:M];
X0      =  zeros(b*b*ch, N*M);
X_m     =  zeros(b*b*ch,length(r)*length(c),'single');

N       =  length(r);
M       =  length(c);
L       =  N*M;

k    =  0;
for i  = 1:b
    for j  = 1:b
        k        =  k+1;        
        for tt =1:ch
            blk      =  im(i:end-b+i,j:end-b+j,tt);
            X0(k+(tt-1)*b*b,:)  =  blk(:)';
        end
    end
end
set        =   1:size(X_m,2);
X_m        =   zeros(length(r)*length(c),b*b*ch,'single');
X = X0';
for i = 1:par.nblk
   v            =  wei_arr(set,i); 
   X_m(set,:)   =  X_m(set,:) + X(blk_arr(set,i),:) .*v(:, ones(1,b2*ch));
end
X_m=X_m';

Cu0         =   zeros(b2*ch, L, 'single' );
b0          =   zeros(b2*ch, L, 'single');
ind         =   zeros(h-b+1,w-b+1);
ind(r,c)    =   1;
X1           =   X0(:, ind~=0);

for i=1:L
    P    =   reshape(Dict(:, i), b*b*ch, b*b*ch);
    coe         =   P*(X0(:, blk_arr(i, :)) - repmat(X_m(:, i), 1, par.s1.nblk) );
    Cu0(:,i)    =   mean(coe.^2, 2);
    b0(:,i)     =   P*(X1(:,i)-X_m(:,i));
end
Cu0      =   max(0, Cu0-par.nSig^2);
b0       =   (Cu0.*(b0.^2)).^(1/4);
Tau     =   par.c./(b0.^2 + par.eps);
return;

