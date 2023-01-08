clc;clear;close all;

% image read
I_lena = double(imread('lena.tif'));
I_sail = double(imread('sail.tif'));

[~,~,c]=size(I_lena);

%% Reconstruction Pipeline
% Wrap Round
% YOUR CODE HERE
I_lena_wrap=padarray(I_lena,[4 4],'symmetric','both');
I_sail_wrap=padarray(I_sail,[4 4],'symmetric','both');

% Resample(subsample)
% YOUR CODE HERE
for i=1:c
    I_lena_wrap_samp1(:,:,i)=resample(I_lena_wrap(:,:,i),1,2,3)';
    I_lena_wrap_samp2(:,:,i)=resample(I_lena_wrap_samp1(:,:,i),1,2,3)';
    I_sail_wrap_samp1(:,:,i)=resample(I_sail_wrap(:,:,i),1,2,3)';
    I_sail_wrap_samp2(:,:,i)=resample(I_sail_wrap_samp1(:,:,i),1,2,3)';
end

% Crop Back
% YOUR CODE HERE
I_lena_samp=I_lena_wrap_samp2(3:end-2,3:end-2,:);
I_sail_samp=I_sail_wrap_samp2(3:end-2,3:end-2,:);

% Wrap Round
% YOUR CODE HERE
I_lena_samp_wrap=padarray(I_lena_samp,[2 2],'symmetric','both'); % the padded broder should be a half of the preivous one
I_sail_sample_wrap=padarray(I_sail_samp,[2 2],'symmetric','both');

% Resample (upsample)
% YOUR CODE HERE
for i=1:c
    I_lena_samp_wrap_samp1(:,:,i)=resample(I_lena_samp_wrap(:,:,i),2,1,3)';
    I_lena_samp_wrap_samp2(:,:,i)=resample(I_lena_samp_wrap_samp1(:,:,i),2,1,3)';
    I_sail_samp_wrap_samp1(:,:,i)=resample(I_sail_sample_wrap(:,:,i),2,1,3)';
    I_sail_samp_wrap_samp2(:,:,i)=resample(I_sail_samp_wrap_samp1(:,:,i),2,1,3)';
end

% Crop back
% YOUR CODE HERE
I_rec_lena=I_lena_samp_wrap_samp2(5:end-4,5:end-4,:);
I_rec_sail=I_sail_samp_wrap_samp2(5:end-4,5:end-4,:);

% Distortion Analysis
PSNR_lena        = calcPSNR(I_lena, I_rec_lena);
PSNR_sail        = calcPSNR(I_sail, I_rec_sail);
fprintf('PSNR lena subsampling = %.3f dB\n', PSNR_lena)
fprintf('PSNR sail subsampling = %.3f dB\n', PSNR_sail)

%% put all the sub-functions called in your script here
function MSE = calcMSE(Image, recImage)
% YOUR CODE HERE
[h,w,c]=size(Image);
Image=double(Image);
recImage=double(recImage);
diff=Image(:)-recImage(:);
MSE=sum((diff).^2)/(h*w*c);
end

function PSNR = calcPSNR(Image, recImage)
% YOUR CODE HERE
MSE=calcMSE(Image, recImage);
PSNR=10*log10((2^8-1)^2/MSE);
end