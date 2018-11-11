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

clc;
clear;
close all;
warning off all

%Micro-polarizator layout(2-by-2 super pixel)
% 0    45
% 135  90
%% load Images
img0 = double( imread('./Test images/Building1/0.bmp') );
img45 = double( imread('./Test images/Building1/45.bmp') );
img90 = double( imread('./Test images/Building1/90.bmp') );
img135 = double( imread('./Test images/Building1/135.bmp') );

%% for building 1, crop image
img0=img0(233:end,89:436);
img45=img45(233:end,89:436);
img90=img90(233:end,89:436);
img135=img135(233:end,89:436);

%% for building 2, crop image
% img0=img0(211:end,231:end);
% img45=img45(211:end,231:end);
% img90=img90(211:end,231:end);
% img135=img135(211:end,231:end);

%%
scale = 2;
par = Parameters_setting();
I = zeros(size(img0,1),size(img0,2),4);
I(:,:,1) = img0;
I(:,:,2) = img45;
I(:,:,3) = img90;
I(:,:,4) = img135;
%% Calculate Stokes Vector
Ori_S0 = (I(:,:,1)+I(:,:,2)+I(:,:,3)+I(:,:,4))*0.5;
Ori_S1 = I(:,:,1)-I(:,:,3);
Ori_S2 = I(:,:,2)-I(:,:,4);
Ori_DoLP = sqrt(Ori_S1.^2+Ori_S2.^2)./Ori_S0;
figure(1);imshow(graystretch(Ori_S0),[]);
figure(2);imshow(graystretch(Ori_DoLP),[]);
%% Parameters setting
par.I = I;
lr = par.I(1:scale:end,1:scale:end,1);
par.LR = zeros(size(lr,1),size(lr,2),4);
par.LR(:,:,1) = lr;
par.LR(:,:,2) = par.I(1:scale:end,2:scale:end,2);
par.LR(:,:,3) = par.I(2:scale:end,2:scale:end,3);
par.LR(:,:,4) = par.I(2:scale:end,1:scale:end,4);
%%
par.D = Get_DownSample_Matrix(par);
par.S0 = Ori_S0;
par.DoLP = Ori_DoLP;
%% Main Code
im = SR_Interpolation_DOFP(par);

imOut0 = im(:,:,1);
imOut45 = im(:,:,2);
imOut90 = im(:,:,3);
imOut135 = im(:,:,4);

S0 = (imOut0+imOut45+imOut90+imOut135)*0.5;
S1 = imOut0-imOut90;
S2 = imOut45-imOut135;
DoLP = sqrt(S1.^2+S2.^2)./S0;
figure;imshow(graystretch(S0),[]);
figure;imshow(graystretch(DoLP),[]);
%%
[res] = Cal_RMSE_PSNR_SSIM(im, I, 5, 5, 255);
fprintf('-----------------Final Result---------------------\n');
fprintf( 'RMSE PSNR SSIM Of <I0 I45 I90 I135>: \n');
disp(res);
[S0, DoLP] = Cal_Stokes(im);
[res] = Cal_RMSE_PSNR_SSIM(S0, par.S0, 5, 5, 255);
fprintf( 'RMSE PSNR SSIM Of S0:\n');
disp(res);
[res] = Cal_RMSE_PSNR_SSIM(DoLP, par.DoLP, 5, 5, 1);
fprintf( 'RMSE PSNR SSIM Of DoLP:\n');
disp(res);
fprintf('============================================================================\n');