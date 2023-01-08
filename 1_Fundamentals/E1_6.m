clc;clear;close all;

% read original RGB image 
I_ori = double(imread('sail.tif'));
[m,n,c]=size(I_ori);

% YOUR CODE HERE for chroma subsampling
I_yuv=ictRGB2YCbCr(I_ori);

% keep luminance unchanged
I_u_wrap=padarray(I_yuv(:,:,2),[4 4],'symmetric','both');
I_v_wrap=padarray(I_yuv(:,:,3),[4 4],'symmetric','both');

% Resample(subsample)
% YOUR CODE HERE
I_u_wrap_samp1=resample(I_u_wrap,1,2,3)';
I_u_wrap_samp2=resample(I_u_wrap_samp1,1,2,3)';
I_v_wrap_samp1=resample(I_v_wrap,1,2,3)';
I_v_wrap_samp2=resample(I_v_wrap_samp1,1,2,3)';

% Crop Back
% YOUR CODE HERE
I_u_samp=I_u_wrap_samp2(3:end-2,3:end-2);
I_v_samp=I_v_wrap_samp2(3:end-2,3:end-2);

% Wrap Round
% YOUR CODE HERE
I_u_samp_wrap=padarray(I_u_samp,[2 2],'symmetric','both'); % the padded broder should be a half of the preivous one
I_v_sample_wrap=padarray(I_v_samp,[2 2],'symmetric','both');

% Resample (upsample)
% YOUR CODE HERE
    I_u_samp_wrap_samp1=resample(I_u_samp_wrap,2,1,3)';
    I_u_samp_wrap_samp2=resample(I_u_samp_wrap_samp1,2,1,3)';
    I_v_samp_wrap_samp1=resample(I_v_sample_wrap,2,1,3)';
    I_v_samp_wrap_samp2=resample(I_v_samp_wrap_samp1,2,1,3)';

% Crop back
% YOUR CODE HERE
I_rec_u=I_u_samp_wrap_samp2(5:end-4,5:end-4);
I_rec_v=I_v_samp_wrap_samp2(5:end-4,5:end-4);

I_rec_yuv(:,:,1)=I_yuv(:,:,1);
I_rec_yuv(:,:,2)=I_rec_u;
I_rec_yuv(:,:,3)=I_rec_v;
I_rec=ictYCbCr2RGB(I_rec_yuv);

% Evaluation
% I_rec is the reconstructed image in RGB color space
PSNR = calcPSNR(I_ori, I_rec);
fprintf('PSNR is %.2f dB\n', PSNR);
% BitRate=(m*n*c*8)/8*c;
% fprintf('BitRate is %.2f\n', BitRate);

figure;
subplot(1,2,1); imshow(uint8(I_ori)); title('original image');
subplot(1,2,2); imshow(uint8(I_rec)); title('reconstructed image by subsampling Cb and Cr');

%% put all the sub-functions called in your script here
function rgb = ictYCbCr2RGB(yuv)
y=yuv(:,:,1);
cb=yuv(:,:,2);
cr=yuv(:,:,3);
rgb(:,:,1)=y+1.402*cr;
rgb(:,:,2)=y-0.344*cb-0.714*cr;
rgb(:,:,3)=y+1.772*cb;
end

function yuv = ictRGB2YCbCr(rgb)
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

function MSE = calcMSE(Image, recImage)
[h,w,c]=size(Image);
Image=double(Image);
recImage=double(recImage);
diff=Image(:)-recImage(:);
MSE=sum((diff).^2)/(h*w*c);
end

function PSNR = calcPSNR(Image, recImage)
MSE=calcMSE(Image, recImage);
PSNR=10*log10((2^8-1)^2/MSE);
end