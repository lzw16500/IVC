clc;clear;close all;

% Read Image
imageLena = double(imread('lena.tif'));
YCbCr_Lena=ictRGB2YCbCr(imageLena); % convert into YCbCr
% create the predictor and obtain the residual image
resImage  = get_residual_YCbCrimage(YCbCr_Lena);
% get the PMF of the residual image
pmfRes    = stats_marg(resImage,-383:383);
% calculate the entropy of the residual image
H_res     = calc_entropy(pmfRes);

fprintf('H_err_OnePixel   = %.2f bit/pixel\n',H_res);

%% build Huffman Code
% input:        PMF               - probabilty mass function of the source
%
% return value:  BinaryTree        - cell structure containing the huffman tree
%                HuffCode          - Array of integers containing the huffman tree
%                BinCode           - Matrix containing the binary version of the code
%                Codelengths       - Array with number of bits in each Codeword

[BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmfRes);
figure;
plot(Codelengths); title('length of codeword');

%% Put all sub-functions which are called in your script here.
function yuv = ictRGB2YCbCr(rgb)
% Input         : rgb (Original RGB Image)
% Output        : yuv (YCbCr image after transformation)
% YOUR CODE HERE
r=rgb(:,:,1);
g=rgb(:,:,2);
b=rgb(:,:,3);
y=0.299*r+0.587*g+0.114*b;
cb=-0.169*r-0.331*g+0.5*b;
cr=0.5*r-0.419*g-0.081*b;
yuv(:,:,1)=y;
yuv(:,:,2)=cb;
yuv(:,:,3)=cr;
end

function res_im= get_residual_YCbCrimage(image)
[m,n,c]=size(image);
res_im=zeros(m,n,c);
predict_im=zeros(m,n,c);
res_im(:,1,:)=image(:,1,:);
res_im(1,:,:)=image(1,:,:);
predict_im(:,1,:)=image(:,1,:);
predict_im(1,:,:)=image(1,:,:);
for i=1:c
    if i==1
        for j=2:m
            for k=2:n
                predict_im(j,k,i)=(7/8)*predict_im(j,k-1,i)-(1/2)*predict_im(j-1,k-1,i)+(5/8)*predict_im(j-1,k,i);
                res_im(j,k,i)=round(image(j,k,i)-predict_im(j,k,i));
                predict_im(j,k,i)=predict_im(j,k,i)+res_im(j,k,i);
            end
        end
    else
        for j=2:m
            for k=2:n
                predict_im(j,k,i)=(3/8)*predict_im(j,k-1,i)-(1/4)*predict_im(j-1,k-1,i)+(7/8)*predict_im(j-1,k,i);
                res_im(j,k,i)=round(image(j,k,i)-predict_im(j,k,i));
                predict_im(j,k,i)=predict_im(j,k,i)+res_im(j,k,i);
            end
        end
    end
end

end

function pmf = stats_marg(image, range)
pmf=hist(image(:),range); % a vector
pmf=pmf/sum(pmf);
end

function H = calc_entropy(pmf)
H=-sum(pmf.*log2(pmf),'omitnan');
end