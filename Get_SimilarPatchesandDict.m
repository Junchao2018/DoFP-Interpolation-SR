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
function [blk_arr, wei_arr, Dict]  =  Get_SimilarPatchesandDict(fs, par)
S  = 21;
f = par.win;
f2 =  f^2;
nv =  par.nblk;
s =  par.step;
hp =  par.hp;
ch =  size(fs,3);

N =  size(fs,1)-f+1;
M =  size(fs,2)-f+1;
r =  [1:s:N];
r =  [r r(end)+1:N];
c =  [1:s:M];
c =  [c c(end)+1:M];
L =  N*M;
X =  zeros(f2*ch, L, 'single');

k    =  0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;        
        for i1 = 1:ch
            blk  =  fs(i:end-f+i,j:end-f+j, i1);
            X(k + (i1-1)*f2,:) =  blk(:)';
        end       
    end
end

I = (1:L);
I = reshape(I, N, M);
N1 = length(r);
M1 = length(c);

blk_arr   =  ones(nv, N1*M1 );
wei_arr   =  zeros(nv, N1*M1 ); 
Dict      =  zeros((f2*ch)^2, N1*M1);
X         =  X';
INDEX = zeros(1, N1*M1 );


for  i  =  1 : N1
    for  j  =  1 : M1
        nv        =  par.nblk;
        row     =   r(i);
        col     =   c(j);
        off     =  (col-1)*N + row;
        off1    =  (j-1)*N1 + i;
                
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
         
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        B       =   X(idx, :);
        v       =   X(off, :);
        
        dis     =   (B(:,1) - v(1)).^2;
        for k = 2:f2*ch
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        dis   =  dis./(f2*ch);
        [val,ind]   =  sort(dis);
        
        maxv = max(val);
        NN = sum(val<maxv*0.4);
        
        if NN<f2*ch
            P = dctmtx(f2*ch);
            Dict(:,off1) = P(:);
            indc        =  idx( ind(2:nv+1) );
            blk_arr(1:nv,off1)  =  indc;
            wei_arr(1:nv,off1)  =  0;
            INDEX(off1) = off;
        else
            if NN >3*f2*ch
                Block = B(ind(1:3*f2*ch),:);
            else
                Block = B(ind(1:NN),:);
            end
            [P, mx] = getpca(Block');
            Dict(:,off1) = P(:);
            wei         =  exp( -dis(ind(2:nv+1))./hp );  
            wei         =  wei./(sum(wei)+eps);
            indc        =  idx( ind(2:nv+1) );
            blk_arr(1:nv,off1)  =  indc;
            wei_arr(1:nv,off1)  =  wei;
            INDEX(off1) = off;
        end    
    end
end
blk_arr  = blk_arr';
wei_arr  = wei_arr';
