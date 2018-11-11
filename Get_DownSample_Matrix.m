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
function D = Get_DownSample_Matrix(par)
[lh,lw,ch] = size(par.LR);
D = cell(ch,1);
s = par.scale;
hh = lh*s;
hw = lw*s;
M = lh*lw;
N = hh*hw;
%%
R = zeros(M,1);
C = zeros(M,1);
V = zeros(M,1);
cnt = 1;
pos = (1:N);
pos = reshape(pos,[hh hw]);
for lrow = 1:lh
    for lcol = 1:lw
        row = (lrow-1)*s + 1;
        col = (lcol-1)*s + 1;
        row_idx = (lcol-1)*lh + lrow;
        col_ind = pos(row, col);
        nn =  1;
        R(cnt:cnt+nn-1) = row_idx;
        C(cnt:cnt+nn-1) = col_ind;
        V(cnt:cnt+nn-1) = 1;
        cnt = cnt+nn;
    end
end
R = R(1:cnt-1);
C = C(1:cnt-1);
V = V(1:cnt-1);
D{1} = sparse(R, C, V, M, N);
%%
R = zeros(M,1);
C = zeros(M,1);
V = zeros(M,1);
cnt = 1;
pos = (1:N);
pos = reshape(pos,[hh hw]);
for lrow = 1:lh
    for lcol = 1:lw
        row = (lrow-1)*s + 1;
        col = (lcol-1)*s + 2;
        row_idx = (lcol-1)*lh + lrow;
        col_ind = pos(row, col);
        nn =  1;
        R(cnt:cnt+nn-1) = row_idx;
        C(cnt:cnt+nn-1) = col_ind;
        V(cnt:cnt+nn-1) = 1;
        cnt = cnt+nn;
    end
end
R = R(1:cnt-1);
C = C(1:cnt-1);
V = V(1:cnt-1);
D{2} = sparse(R, C, V, M, N);
%%
R = zeros(M,1);
C = zeros(M,1);
V = zeros(M,1);
cnt = 1;
pos = (1:N);
pos = reshape(pos,[hh hw]);
for lrow = 1:lh
    for lcol = 1:lw
        row = (lrow-1)*s + 2;
        col = (lcol-1)*s + 2;
        row_idx = (lcol-1)*lh + lrow;
        col_ind = pos(row, col);
        nn =  1;
        R(cnt:cnt+nn-1) = row_idx;
        C(cnt:cnt+nn-1) = col_ind;
        V(cnt:cnt+nn-1) = 1;
        cnt = cnt+nn;
    end
end
R = R(1:cnt-1);
C = C(1:cnt-1);
V = V(1:cnt-1);
D{3} = sparse(R, C, V, M, N);
%%
R = zeros(M,1);
C = zeros(M,1);
V = zeros(M,1);
cnt = 1;
pos = (1:N);
pos = reshape(pos,[hh hw]);
for lrow = 1:lh
    for lcol = 1:lw
        row = (lrow-1)*s + 2;
        col = (lcol-1)*s + 1;
        row_idx = (lcol-1)*lh + lrow;
        col_ind = pos(row, col);
        nn =  1;
        R(cnt:cnt+nn-1) = row_idx;
        C(cnt:cnt+nn-1) = col_ind;
        V(cnt:cnt+nn-1) = 1;
        cnt = cnt+nn;
    end
end
R = R(1:cnt-1);
C = C(1:cnt-1);
V = V(1:cnt-1);
D{4} = sparse(R, C, V, M, N);
end

