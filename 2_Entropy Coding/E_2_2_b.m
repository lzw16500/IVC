clc;clear;close all;

% Read Image
imageLena = double(imread('lena.tif'));
% Calculate Joint PMF
jpmfLena  = stats_joint(imageLena);
% Calculate Joint Entropy
Hjoint    = calc_entropy(jpmfLena);
fprintf('H_joint = %.2f bit/pixel pair\n', Hjoint);

mesh(jpmfLena);

%% Put all sub-functions which are called in your script here.
function pmf = stats_joint(image)
%  Input         : image (Original Image)
%
%  Output        : pmf   (Probability Mass Function)
pmf=zeros(256,256);
s_left=image(:,1:2:end);
s_right=image(:,2:2:end);
s_left=s_left(:);
s_right=s_right(:);
for i=1:numel(s_left)
        pmf(s_left(i)+1,s_right(i)+1)=pmf(s_left(i)+1,s_right(i)+1)+1;
end
pmf=pmf/sum(pmf(:));
end

function H = calc_entropy(pmf)
%  Input         : pmf   (Probability Mass Function)
%
%  Output        : H     (Entropy in bits)
pmf_wo0=pmf(pmf>0);
H=-sum(pmf_wo0.*log2(pmf_wo0));
end