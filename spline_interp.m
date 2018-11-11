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
function [imgout] = spline_interp(img)
%   0   45
%   135 90
%%
[m, n, ch]=size(img);
M = 2*m;
N = 2*n;
imgout = zeros(M,N,ch);
img1=zeros(M,N);
img2=zeros(M,N);
img3=zeros(M,N);
img4=zeros(M,N);
img1(1:2:M,1:2:N)=img(:,:,1);
img2(1:2:M,2:2:N)=img(:,:,2);
img3(2:2:M,2:2:N)=img(:,:,3);
img4(2:2:M,1:2:N)=img(:,:,4);
%%
m = M;
n = N;
%0бу
x=1:2:n;
y=img1(1:2:m,1:2:n);
pp=csapi(x,y);
img1(1:2:m,1:n-1)=fnval(pp,1:n-1);
x1=1:2:m;
y1=img1(1:2:m,1:n-1);
y1=y1';
pp=csapi(x1,y1);
img1(1:m-1,1:n-1)=(fnval(pp,1:m-1))';
img1(m,:)=img1(m-1,:);
img1(:,n)=img1(:,n-1);
%45бу
x=2:2:n;
y=img2(1:2:m,2:2:n);
pp=csapi(x,y);
img2(1:2:m,2:n)=fnval(pp,2:n);
x1=1:2:m;
y1=img2(1:2:m,2:n);
y1=y1';
pp=csapi(x1,y1);
img2(1:m-1,2:n)=(fnval(pp,1:m-1))';
img2(m,:)=img2(m-1,:);
img2(:,1)=img2(:,2);
%90бу
x=2:2:n;
y=img3(2:2:m,2:2:n);
pp=csapi(x,y);
img3(2:2:m,2:n)=fnval(pp,2:n);
x1=2:2:m;
y1=img3(2:2:m,2:n);
y1=y1';
pp=csapi(x1,y1);
img3(2:m,2:n)=(fnval(pp,2:m))';
img3(1,:)=img3(2,:);
img3(:,1)=img3(:,2);
%135бу
x=1:2:n;
y=img4(2:2:m,1:2:n);
pp=csapi(x,y);
img4(2:2:m,1:n-1)=fnval(pp,1:n-1);
x1=2:2:m;
y1=img4(2:2:m,1:n);
y1=y1';
pp=csapi(x1,y1);
img4(2:m,1:n)=(fnval(pp,2:m))';
img4(1,:)=img4(2,:);
img4(:,n)=img4(:,n-1);
%%
imgout(:,:,1) = img1;
imgout(:,:,2) = img2;
imgout(:,:,3) = img3;
imgout(:,:,4) = img4;
end

