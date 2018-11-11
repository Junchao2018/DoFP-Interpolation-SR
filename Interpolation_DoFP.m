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
function [hr_im] = Interpolation_DoFP( lr_im, par, ori_im, hr_im0, L)
hr_im = hr_im0;
[hh, hw, ch] = size(hr_im);
%% Calculate RMSE PSNR SSIM
fprintf('==============================SR Interpolation==============================\n');
[res] = Cal_RMSE_PSNR_SSIM(hr_im, ori_im, 5, 5, 255);
fprintf( '\nInitial RMSE PSNR SSIM Of <I0 I45 I90 I135>: \n');
disp(res);
[S0, DoLP] = Cal_Stokes(hr_im);
[res] = Cal_RMSE_PSNR_SSIM(S0, par.S0, 5, 5, 255);
fprintf( 'Initial RMSE PSNR SSIM Of S0:\n');
disp(res);
[res] = Cal_RMSE_PSNR_SSIM(DoLP, par.DoLP, 5, 5, 1);
fprintf( 'Initial RMSE PSNR SSIM Of DoLP:\n');
disp(res);
fprintf('---------------------------------------------------\n');
%%
gamma = par.gamma;
y = lr_im;
D = par.D;
cnt = 0;
im0 = hr_im0;
Y = zeros( numel(y(:,:,1)), ch);
for i = 1:ch
    tmp = y(:,:,i);
    Y(:,i) = tmp(:);
end
mu = 1.7;
for k   =  1:L
    [blk_arr, wei_arr, Dict] = Get_SimilarPatchesandDict( im0, par.s1);
    Tau = UpdateParameters(im0, Dict, par, blk_arr, wei_arr);
    S = @(x) SP_Interp( x, Dict, par, Tau./gamma, blk_arr, wei_arr);
    f = hr_im;
    for  iter = 1 : par.iters
        if (mod(cnt, 10) == 0 && cnt~=0)
            if isfield(par,'I')
                [res] = Cal_RMSE_PSNR_SSIM(f, ori_im, 5, 5, 255);
                fprintf('\nIter %d :\n',cnt);
                fprintf( 'RMSE PSNR SSIM Of <I0 I45 I90 I135>: \n');
                disp(res);
                [S0, DoLP] = Cal_Stokes(f);
                [res] = Cal_RMSE_PSNR_SSIM(S0, par.S0, 5, 5, 255);
                fprintf( 'RMSE PSNR SSIM Of S0:\n');
                disp(res);
                [res] = Cal_RMSE_PSNR_SSIM(DoLP, par.DoLP, 5, 5, 1);
                fprintf( 'RMSE PSNR SSIM Of DoLP:\n');
                disp(res);
                fprintf('---------------------------------------------------\n');
            end
        end    
        [f, T1, P1] = S( f );
        for i = 1:ch
            tmp = T1(:,:,i);
            TT(:,i) = tmp(:);
        end
        for i =1:ch
            DTD = D{i}'*D{i};
            b = gamma*TT(:,i) + mu*D{i}'*(Y(:,i));
            f1(:,:,i) = reshape( cgsolve2( @(x) fun(x, gamma, P1(:,:,i), size(f(:,:,i)), mu, DTD), b), hh, hw);
        end
        mu = mu*1.1;
        f = f1;
        cnt = cnt  +  1;
    end
    hr_im =  f;
    im0 =  f;
    gamma = gamma*1.25;
end
function  y  =  fun(x, gamma, P1, sz, mu, DTD)
y      =   gamma*P1.*reshape(x, sz);
y      =   y(:) + mu*DTD*x;
return;


