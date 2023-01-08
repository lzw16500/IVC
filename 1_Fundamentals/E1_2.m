clc;clear;close all;

% Read Image
I = double(imread('satpic1.bmp'));
[m,n,c]=size(I);

% YOUR CODE HERE
kernel=[1 2 1;2 4 2;1 2 1];
kernel=kernel/(sum(kernel(:)));  % normalization of the filter

%% b. plot the frequency response of the filter
figure;
imagesc(abs(fftshift(fft2(kernel)))); % 
title('Frequency response of the filter');

%% Show the effect of filtering
I_filtered=prefilterlowpass2d(I,kernel);
figure;
subplot(1,2,1); imshow(uint8(I)); title('Original Image')
subplot(1,2,2); imshow(uint8(I_filtered)); title('Filtered image');
% subplot(1,2,1); imshow(I); title('Original Image')
% subplot(1,2,2); imshow(I_filtered); title('Filtered image');


%% compare the PSNR of reconstructed Images with/ without prefiltering
% without prefiltering
I_filtered_down1=zeros(512,256,3);
I_filtered_down2=zeros(256,256,3);
I_filtered_up1=zeros(256,512,3);
I_filtered_up2=zeros(512,512,3);
for i=1:c
    I_filtered_down1(:,:,i)=downsample(I(:,:,i),2)';
    I_filtered_down2(:,:,i)=downsample(I_filtered_down1(:,:,i),2)'; % get the image after downsampling
end
for i=1:c
    I_filtered_up1(:,:,i)=upsample(I_filtered_down2(:,:,i),2)';
    I_filtered_up2(:,:,i)=upsample(I_filtered_up1(:,:,i),2)'; % get the image after upsampling
end
I_rec_notpre=prefilterlowpass2d(I_filtered_up2,kernel);
I_rec_notpre=I_rec_notpre*4; %energy

% Evaluation without prefiltering
% I_rec_notpre is the reconstructed image WITHOUT prefiltering
PSNR_notpre = calcPSNR(I, I_rec_notpre);
fprintf('Reconstructed image, not prefiltered, PSNR = %.2f dB\n', PSNR_notpre)

% with prefiltering
% YOUR CODE HERE
% picture_pad=padarray(I,[1 1],'symmetric','both'); % pad 3-dimensional array
I_filtered_down1=zeros(512,256,3);
I_filtered_down2=zeros(256,256,3);
I_filtered_up1=zeros(256,512,3);
I_filtered_up2=zeros(512,512,3);
for i=1:c
    I_filtered_down1(:,:,i)=downsample(I_filtered(:,:,i),2)';
    I_filtered_down2(:,:,i)=downsample(I_filtered_down1(:,:,i),2)'; % get the image after downsampling
end
% I_rec_pre_down=I_rec_pre_down2(2:end,2:end,:); %crop, s.t. the final downsampled image is 1/4 of original imoage in szie
% picture_pad=padarray(I_rec_pre_down,[1 1],'symmetric','both'); % pad 3-dimensional array
% I_rec_pre_up1=zeros(258,516,c);
% I_rec_pre_up2=zeros(516,516,c);
for i=1:c
    I_filtered_up1(:,:,i)=upsample(I_filtered_down2(:,:,i),2)';
    I_filtered_up2(:,:,i)=upsample(I_filtered_up1(:,:,i),2)'; % get the image after downsampling
end
I_rec_pre=prefilterlowpass2d(I_filtered_up2,kernel);
I_rec_pre=I_rec_pre*4; %energy
% I_rec_pre=I_filtered(2:end-1,2:end-1,:);

% Evaluation with prefiltering
% I_rec_pre is the reconstructed image WITH prefiltering
PSNR_pre = calcPSNR(I, I_rec_pre);
fprintf('Reconstructed image, prefiltered, PSNR = %.2f dB\n', PSNR_pre)

figure;
subplot(1,2,1); imshow(uint8(I_rec_notpre)); title('reconstructed image without prefiltering');
subplot(1,2,2); imshow(uint8(I_rec_pre)); title('reconstructed image with prefiltering');

%% put all the sub-functions called in your script here
function pic_pre = prefilterlowpass2d(picture, kernel)
% YOUR CODE HERE
[~,~,c]=size(picture);
% picture_pad=padarray(picture,[1 1],'symmetric','both'); % pad 3-dimensional array
for i=1:c
    pic_pre(:,:,i)=conv2(picture(:,:,i),kernel,'same');
end
% pic_pre=picture_pad;
end

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