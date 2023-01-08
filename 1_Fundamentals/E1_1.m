% read images
Image_lena = imread('lena.tif');
Image_lena_gray = rgb2gray(Image_lena);
Image_smandril = imread('smandril.tif');
Image_smandril_gray = rgb2gray(Image_smandril);

% show images in a subplot figure
subplot(221), imshow(Image_lena), title('Original Image Lena')
subplot(222), imshow(Image_lena_gray), title('Compressed Image Lena')
subplot(223), imshow(Image_smandril), title('Original Image Smandril')
subplot(224), imshow(Image_smandril_gray), title('Compressed Image Smandril')

% check the property of variables
whos I*