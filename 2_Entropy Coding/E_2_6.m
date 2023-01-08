clc;clear;close all;

% read image
imageLena       = double(imread('lena.tif'));
imageLena_small = double(imread('lena_small.tif'));

%% calculate the residual from lena_small
YCbCr_Lena_small=ictRGB2YCbCr(imageLena_small); % convert into YCbCr
% create the predictor and obtain the residual image
resImage_small  = get_residual_YCbCrimage(YCbCr_Lena_small);
% get the PMF of the residual image
pmfRes_small    = stats_marg(resImage_small,-383:383);

%% calculate the residual from lena
YCbCr_Lena=ictRGB2YCbCr(imageLena); % convert into YCbCr
[Y_lena,down_Cb_lena,down_Cr_lena]=YCbCr_downsample(YCbCr_Lena,1,2,3);

% create the predictor and obtain the residual image
res_Y_Lena=get_residual_Y(Y_lena);
res_down_Cb_lena=get_residual_C(down_Cb_lena);
res_down_Cr_lena=get_residual_C(down_Cr_lena);


%% Coding
% codebook construction
[BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmfRes_small);  % build huffman table of residual of small_lena

% Encoding
bytestream_y = enc_huffman_new( round(res_Y_Lena(:))+383+1, BinCode, Codelengths); % use the above huffman table to encode the residual of the original Lena 
bytestream_cb = enc_huffman_new( round(res_down_Cb_lena(:))+383+1, BinCode, Codelengths);
bytestream_cr = enc_huffman_new( round(res_down_Cr_lena(:))+383+1, BinCode, Codelengths);

% Decoding
decoded_res_lena_y = double(reshape( dec_huffman_new ( bytestream_y, BinaryTree, max(size(res_Y_Lena(:))) ), size(res_Y_Lena))) -384;
decoded_res_lena_cb = double(reshape( dec_huffman_new ( bytestream_cb, BinaryTree, max(size(res_down_Cb_lena(:))) ), size(res_down_Cb_lena))) -384;
decoded_res_lena_cr = double(reshape( dec_huffman_new ( bytestream_cr, BinaryTree, max(size(res_down_Cr_lena(:))) ), size(res_down_Cr_lena))) -384;

%% Reconstruction
[rec_y,rec_down_cb,rec_down_cr]=reconstruction_from_res(decoded_res_lena_y,decoded_res_lena_cb,decoded_res_lena_cr);
rec_image=YCbCr_upsample(rec_y,rec_down_cb,rec_down_cr,2,1,3);
rec_image=ictYCbCr2RGB(rec_image);

%% evaluation and show results
figure
subplot(121)
imshow(uint8(imageLena)), title('Original Image')
subplot(122)

PSNR = calcPSNR(imageLena, rec_image);
imshow(uint8(rec_image)), title(sprintf('Reconstructed Image, PSNR = %.2f dB', PSNR))
bits=(numel(bytestream_y)+numel(bytestream_cr)+numel(bytestream_cr))*8;
BPP = bits / (numel(imageLena)/3);
CompressionRatio = 24/BPP;

fprintf('Bit Rate         = %.2f bit/pixel\n', BPP);
fprintf('CompressionRatio = %.2f\n', CompressionRatio);
fprintf('PSNR             = %.2f dB\n', PSNR);

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

function rgb = ictYCbCr2RGB(yuv)
% Input         : yuv (Original YCbCr image)
% Output        : rgb (RGB Image after transformation)
% YOUR CODE HERE
y=yuv(:,:,1);
cb=yuv(:,:,2);
cr=yuv(:,:,3);
rgb(:,:,1)=y+1.402*cr;
rgb(:,:,2)=y-0.344*cb-0.714*cr;
rgb(:,:,3)=y+1.772*cb;
end

function res_im= get_residual_YCbCrimage(image) % res_im has the same size as image
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

function res_y=get_residual_Y(y)
%input: Y
res_y=zeros(size(y));
pred_y=zeros(size(y));
[m,n]=size(y);
res_y(1,:)=y(1,:);
res_y(:,1)=y(:,1);
pred_y(1,:)=y(1,:);
pred_y(:,1)=y(:,1);
for i=2:m
    for j=2:n
        pred_y(i,j)=(7/8)*pred_y(i,j-1)-(1/2)*pred_y(i-1,j-1)+(5/8)*pred_y(i-1,j);
        res_y(i,j)=round(y(i,j)-pred_y(i,j));
        pred_y(i,j)=pred_y(i,j)+res_y(i,j);
    end
end    
end

function res_c=get_residual_C(c)
% input: Cb or Cr
res_c=zeros(size(c));
pred_c=zeros(size(c));
[m,n]=size(c);
res_c(1,:)=c(1,:);
res_c(:,1)=c(:,1);
pred_c(1,:)=c(1,:);
pred_c(:,1)=c(:,1);
for i=2:m
    for j=2:n
        pred_c(i,j)=(3/8)*pred_c(i,j-1)-(1/4)*pred_c(i-1,j-1)+(7/8)*pred_c(i-1,j);
        res_c(i,j)=round(c(i,j)-pred_c(i,j));
        pred_c(i,j)=pred_c(i,j)+res_c(i,j);
    end
end   
end

function [rec_y,rec_down_cb,rec_down_cr]=reconstruction_from_res(res_y,res_cb,res_cr) % reconstruct the image based on the residual matrix
[m1,n1]=size(res_y);
[m2,n2]=size(res_cb);
rec_y=zeros(m1,n1);
rec_down_cb=zeros(m2,n2);
rec_down_cr=zeros(m2,n2);
rec_y(1,:)=res_y(1,:);
rec_down_cb(1,:)=res_cb(1,:);
rec_down_cr(1,:)=res_cr(1,:);
rec_y(:,1)=res_y(:,1);
rec_down_cb(:,1)=res_cb(:,1);
rec_down_cr(:,1)=res_cr(:,1);

for i=2:m1
    for j=2:n1
        rec_y(i,j)=(7/8)*rec_y(i,j-1)-(1/2)*rec_y(i-1,j-1)+(5/8)*rec_y(i-1,j)+res_y(i,j);
    end
end

for i=2:m2
    for j=2:n2
        rec_down_cb(i,j)=(3/8)*rec_down_cb(i,j-1)-(1/4)*rec_down_cb(i-1,j-1)+(7/8)*rec_down_cb(i-1,j)+res_cb(i,j);
        rec_down_cr(i,j)=(3/8)*rec_down_cr(i,j-1)-(1/4)*rec_down_cr(i-1,j-1)+(7/8)*rec_down_cr(i-1,j)+res_cr(i,j);
    end
end

end

function pmf = stats_marg(image, range)
pmf=hist(image(:),range); % a vector
pmf=pmf/sum(pmf);
end

function MSE = calcMSE(Image, recImage)
% Input         : Image    (Original Image)
%                 recImage (Reconstructed Image)
% Output        : MSE      (Mean Squared Error)
% YOUR CODE HERE

[h,w,c]=size(Image);
Image=double(Image);
recImage=double(recImage);
diff=Image(:)-recImage(:);
MSE=sum((diff).^2)/(h*w*c);
end

function PSNR = calcPSNR(Image, recImage)
% Input         : Image    (Original Image)
%                 recImage (Reconstructed Image)
%
% Output        : PSNR     (Peak Signal to Noise Ratio)
% YOUR CODE HERE
% call calcMSE to calculate MSE

MSE=calcMSE(Image, recImage);
PSNR=10*log10((2^8-1)^2/MSE);
end

function [y,cb,cr]=YCbCr_downsample(yuv_im,p,q,filter_length)
yuv_wrap=padarray(yuv_im,[4 4],'symmetric','both');
y=yuv_im(:,:,1);
for i=2:3
    yuv_wrap_samp1(:,:,i)=resample(yuv_wrap(:,:,i),p,q,filter_length)';
    yuv_wrap_samp2(:,:,i)=resample(yuv_wrap_samp1(:,:,i),p,q,filter_length)';
end
cb=yuv_wrap_samp2(3:end-2,3:end-2,2);
cr=yuv_wrap_samp2(3:end-2,3:end-2,3);
end

function upsampled_yuv=YCbCr_upsample(y,cb,cr,p,q,filter_length)
cb_wrap=padarray(cb,[2 2],'symmetric','both');
cr_wrap=padarray(cr,[2 2],'symmetric','both');

cb_wrap_samp1=resample(cb_wrap,p,q,filter_length)';
cb_wrap_samp2=resample(cb_wrap_samp1,p,q,filter_length)';

cr_wrap_samp1=resample(cr_wrap,p,q,filter_length)';
cr_wrap_samp2=resample(cr_wrap_samp1,p,q,filter_length)';

upsampled_yuv(:,:,1)=y;
upsampled_yuv(:,:,2)=cb_wrap_samp2(5:end-4,5:end-4);
upsampled_yuv(:,:,3)=cr_wrap_samp2(5:end-4,5:end-4);
end


