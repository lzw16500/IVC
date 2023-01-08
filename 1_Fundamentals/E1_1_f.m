% Read Images
Image_lena = imread('lena.tif');
Image_lena_compressed = imread('lena_compressed.tif');
Image_monarch = imread('monarch.tif');
Image_monarch_compressed = imread('monarch_compressed.tif');

% YOUR CODE HERE
% do NOT change the name of variables (PSNR_lena, PSNR_monarch), the assessment code will check these values with our reference answers, same for all the script assignment.
PSNR_lena = calcPSNR(Image_lena, Image_lena_compressed);
PSNR_monarch = calcPSNR(Image_monarch, Image_monarch_compressed);

fprintf('PSNR of lena.tif is %.3f dB\n', PSNR_lena)
fprintf('PSNR of monarch.tif is %.3f dB\n', PSNR_monarch)

subplot(221), imshow(Image_lena), title('Original Image Lena')
subplot(222), imshow(Image_lena_compressed), title('Compressed Image Lena')
subplot(223), imshow(Image_monarch), title('Original Image Monarch')
subplot(224), imshow(Image_monarch_compressed), title('Compressed Image Monarch')

% put all the sub-functions called in your script here
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
MSE=calcMSE(Image, recImage);
PSNR=10*log10((2^8-1)^2/MSE);
end