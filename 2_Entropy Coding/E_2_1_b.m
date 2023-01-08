clc;clear;close all;

imageLena     = double(imread('lena.tif'));
imageSail     = double(imread('sail.tif'));
imageSmandril = double(imread('smandril.tif'));

pmfLena       = stats_marg(imageLena, 0:255);
HLena         = calc_entropy(pmfLena);

pmfSail       = stats_marg(imageSail, 0:255);
HSail         = calc_entropy(pmfSail);

pmfSmandril   = stats_marg(imageSmandril, 0:255);
HSmandril     = calc_entropy(pmfSmandril);

fprintf('--------------Using individual code table--------------\n');
fprintf('lena.tif      H = %.2f bit/pixel\n', HLena);
fprintf('sail.tif      H = %.2f bit/pixel\n', HSail);
fprintf('smandril.tif  H = %.2f bit/pixel\n', HSmandril);



%% Put all sub-functions which are called in your script here.
function pmf = stats_marg(image, range)
pmf=hist(image(:),range); % a vector
pmf=pmf/sum(pmf);
end

function H = calc_entropy(pmf)
pmf_wo0=pmf(pmf>0);
H=-sum(pmf_wo0.*log2(pmf_wo0));
end

