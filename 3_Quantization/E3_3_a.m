clc;clear;

% general script, e.g. image loading, function calls, etc.
epsilon          = 0.001;

imageLena_small  = double(imread('lena_small.tif'));
[qImageLena_small, clusters_small] = LloydMax(imageLena_small, 3, epsilon);
recImage_small   = InvLloydMax(qImageLena_small, clusters_small);
PSNR_small       = calcPSNR(imageLena_small, recImage_small);

imageLena  = double(imread('lena.tif'));
[qImageLena, clusters] = LloydMax(imageLena, 3, epsilon);
recImageLena   = InvLloydMax(qImageLena, clusters);
PSNR       = calcPSNR(imageLena, recImageLena);


%% define your functions, e.g. calcPSNR, LloydMax, InvLloydMax
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