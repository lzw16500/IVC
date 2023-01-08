clc;clear;close all;

% Read Image
imageLena = double(imread('lena.tif'));
% Calculate Joint Entropy
H_cond    = stats_cond(imageLena);
fprintf('H_cond = %.2f bit/pixel\n', H_cond);


function H = stats_cond(image)
%  Input         : image (Original Image)
%
%  Output        : H   (Conditional Entropy)

[m,n,c]=size(image);
joint_pmf=zeros(256,256);
cond_pmf=zeros(256,256);
for i=1:c
    for j=1:m
        for k=1:n-1
            joint_pmf(image(j,k,i)+1,image(j,k+1,i)+1)=joint_pmf(image(j,k,i)+1,image(j,k+1,i)+1)+1;
        end
    end
end

joint_pmf=joint_pmf/sum(joint_pmf(:));

for i=1:256
    if sum(joint_pmf(i,:))~=0
    cond_pmf(i,:)=joint_pmf(i,:)/sum(joint_pmf(i,:));
    end
end

cond_pmf=cond_pmf(cond_pmf>0);
joint_pmf=joint_pmf(joint_pmf>0);
H=-sum(joint_pmf.*log2(cond_pmf));

end