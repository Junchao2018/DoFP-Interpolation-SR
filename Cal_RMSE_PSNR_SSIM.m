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
function [result] = Cal_RMSE_PSNR_SSIM(im, ori_im, row, col, scale)
[h,w,ch]=size(im);
rmse = zeros(ch,1);
psnr = zeros(ch,1);
Ssim = zeros(ch,1);
for i =1:ch
    e = im(:,:,i) - ori_im(:,:,i);
    e = e(row+1:h-row,col+1:w-col);
    me = mean(mean(e.^2));
    rmse(i) = sqrt(me);
    psnr(i) = 10*log10(scale^2/me);
    Ssim(i) = ssim(im(row+1:end-row,col+1:end-col,i),ori_im(row+1:end-row,col+1:end-col,i));
end
result = [rmse psnr Ssim];

