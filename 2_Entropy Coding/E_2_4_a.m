clc;clear;close all;

% Read Image
imageLena = double(imread('lena.tif'));
% create the predictor and obtain the residual image
resImage  = get_residual_image(imageLena);
% get the PMF of the residual image
pmfRes    = stats_marg(resImage,-255:255);
% calculate the entropy of the residual image
H_res     = calc_entropy(pmfRes);

fprintf('H_err_OnePixel   = %.2f bit/pixel\n',H_res);

%% Put all sub-functions which are called in your script here.
function res_im= get_residual_image(image)
[m,n,c]=size(image);
res_im=zeros(m,n,c);
for i=1:c
    res_im(:,1,i)=image(:,1,i);
    for j=1:m
        for k=2:n
            res_im(j,k,i)=image(j,k,i)-image(j,k-1,i);
        end
    end
end
      
end

function pmf = stats_marg(image, range)
pmf=hist(image(:),range); % a vector
pmf=pmf/sum(pmf);
end

function H = calc_entropy(pmf)
pmf_wo0=pmf(pmf>0);
H=-sum(pmf_wo0.*log2(pmf_wo0));
end
