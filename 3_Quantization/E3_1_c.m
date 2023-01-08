clc;clear;
% general script, e.g. image loading, function calls, etc.
imageLena_small = double(imread('lena_small.tif'));
imageLena = double(imread('lena.tif'));
bits_small      = [1 2 3 5 7];
bits = [3 5];
PSNR_small = [];
for bit = bits_small
    qImageLena_small = UniQuant(imageLena_small, bit);
    recImage_small   = InvUniQuant(qImageLena_small, bit);
    PSNR_small = [PSNR_small calcPSNR(imageLena_small, recImage_small)];
end

PSNR = [];
for bit = bits
    qImageLena = UniQuant(imageLena, bit);
    recImage   = InvUniQuant(qImageLena, bit);
    PSNR = [PSNR calcPSNR(imageLena, recImage)];
end


%% define your functions, e.g. calcPSNR, UniQuant, InvUniQuant
function qImage = UniQuant(image, bits)
norm_image=image/256;
qImage=floor(norm_image*(2^bits));
end

function image = InvUniQuant(qImage, bits)
norm_image=(qImage+0.5)/(2^bits);
image=floor(norm_image*256);
end

function PSNR = calcPSNR( image1, image2 )
MSE=calcMSE(image1, image2);
PSNR=10*log10((2^8-1)^2/MSE);
end

function MSE = calcMSE( image1, image2 )
[h,w,c]=size(image1);
image1=double(image1);
image2=double(image2);
diff=image1(:)-image2(:);
MSE=sum((diff).^2)/(h*w*c);
end