clc;clear;close all;

lena_small = double(imread('lena_small.tif'));
Lena       = double(imread('lena.tif'));

scales = 1 : 0.6 : 1; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1
for scaleIdx = 1 : numel(scales)
    qScale   = scales(scaleIdx);
    k_small  = IntraEncode(lena_small, qScale);
    k        = IntraEncode(Lena, qScale);
    %% use pmf of k_small to build and train huffman table
    %your code here
    pmf_ind=hist(k_small(:),min(k):max(k));
    pmf_ind=pmf_ind/sum(pmf_ind);
    [BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind);
    %% use trained table to encode k to get the bytestream
    % your code here
    bytestream = enc_huffman_new(k-min(k)+1, BinCode, Codelengths);
    bitPerPixel(scaleIdx) = (numel(bytestream)*8) / (numel(Lena)/3);
    %% image reconstruction
    k_rec=double(dec_huffman_new(bytestream,BinaryTree,size(k(:),1)))+min(k)-1;
    I_rec = IntraDecode(k_rec, size(Lena),qScale);
    PSNR(scaleIdx) = calcPSNR(Lena, I_rec);
    fprintf('QP: %.1f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, bitPerPixel(scaleIdx), PSNR(scaleIdx))
end
figure;
plot(bitPerPixel,PSNR);
xlabel('bit/pixel');
ylabel('PSNR [dB]');
%% put all used sub-functions here.
function dst = IntraDecode(image, img_size , qScale)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image)
m=img_size(1);
n=img_size(2);
ch=img_size(3);
dst_temp=zeros(m,n,ch);
eob=1000;
zz=ZeroRunDec_EoB(image,eob);
count=1;
for a=1:8:m
    for b=1:8:n
        temp_zz=reshape(zz(count:count+191),64,3);
        dst_temp(a:a+7,b:b+7,:)=DeZigZag8x8(temp_zz);
        count=count+192;
    end
end
dct_coeff=blockproc(dst_temp,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale));
dst=blockproc(dct_coeff,[8,8],@(block_struct) IDCT8x8(block_struct.data));
dst = ictYCbCr2RGB(dst);
end

function dst = IntraEncode(image, qScale)
%  Function Name : IntraEncode.m
%  Input         : image (Original RGB Image)
%                  qScale(quantization scale)
%  Output        : dst   (sequences after zero-run encoding, 1xN)
image=ictRGB2YCbCr(image);
image_dct=blockproc(image,[8,8],@(block_struct) DCT8x8(block_struct.data));
image_dct_quantized=blockproc(image_dct,[8,8],@(block_struct) Quant8x8(block_struct.data,qScale));
zig_zag=blockproc(image_dct_quantized,[8 8],@(block_struct) ZigZag8x8(block_struct.data)); % zig_zag is still a matrix, with each block of size 64*3
eob=1000;
zz=[];
for i=1:64:size(zig_zag,1)
    for j=1:3:size(zig_zag,2)
        zigzag_block=zig_zag(i:i+63,j:j+2);
        zigzag_vec=zigzag_block(:)';
        zz=[zz,zigzag_vec];
    end
end
dst=ZeroRunEnc_EoB(zz,eob);
end

%% and many more functions
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

function coeff = DCT8x8(block)
%  Input         : block    (Original Image block, 8x8x3)
%
%  Output        : coeff    (DCT coefficients after transformation, 8x8x3)

coeff=zeros(size(block));
for i=1:3
    coeff(:,:,i)=dct2(block(:,:,i));
end

end

function block = IDCT8x8(coeff)
%  Function Name : IDCT8x8.m
%  Input         : coeff (DCT Coefficients) 8*8*3
%  Output        : block (original image block) 8*8*3

block=zeros(size(coeff));
for i=1:3
    block(:,:,i)=idct2(coeff(:,:,i));
end

end

function quant = Quant8x8(dct_block, qScale)
%  Input         : dct_block (Original Coefficients, 8x8x3)
%                  qScale (Quantization Parameter, scalar)
%
%  Output        : quant (Quantized Coefficients, 8x8x3)

Y_table = [16 11 10 16 24 40 51 61; 
           12 12 14 19 26 58 60 55; 
           14 13 16 24 40 57 69 56; 
           14 17 22 29 51 87 80 62; 
           18 55 37 56 68 109 103 77; 
           24 35 55 64 81 104 113 92; 
           49 64 78 87 103 121 120 101; 
           72 92 95 98 112 100 103 99];
CbCr_table = [17 18 24 47 99 99 99 99; 
              18 21 26 66 99 99 99 99; 
              24 13 56 99 99 99 99 99; 
              47 66 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99];
          
quant=zeros(size(dct_block));
quant(:,:,1)=round(dct_block(:,:,1)./(Y_table*qScale));
for i=2:3
    quant(:,:,i)=round(dct_block(:,:,i)./(CbCr_table*qScale));
end

end

function dct_block = DeQuant8x8(quant_block, qScale)
%  Function Name : DeQuant8x8.m
%  Input         : quant_block  (Quantized Block, 8x8x3)
%                  qScale       (Quantization Parameter, scalar)
%
%  Output        : dct_block    (Dequantized DCT coefficients, 8x8x3)

Y_table = [16 11 10 16 24 40 51 61; 
           12 12 14 19 26 58 60 55; 
           14 13 16 24 40 57 69 56; 
           14 17 22 29 51 87 80 62; 
           18 55 37 56 68 109 103 77; 
           24 35 55 64 81 104 113 92; 
           49 64 78 87 103 121 120 101; 
           72 92 95 98 112 100 103 99];
CbCr_table = [17 18 24 47 99 99 99 99; 
              18 21 26 66 99 99 99 99; 
              24 13 56 99 99 99 99 99; 
              47 66 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99; 
              99 99 99 99 99 99 99 99];
          
dct_block=zeros(size(quant_block));
dct_block(:,:,1)=quant_block(:,:,1).*(Y_table*qScale);
for i=2:3
    dct_block(:,:,i)=quant_block(:,:,i).*(CbCr_table*qScale);
end

end

function zz = ZigZag8x8(quant)
%  Input         : quant (Quantized Coefficients, 8x8x3)
%
%  Output        : zz (zig-zag scaned Coefficients, 64x3)

Zigzag = [1 2 6 7 15 16 28 29;
    3 5 8 14 17 27 30 43;
    4 9 13 18 26 31 42 44;
    10 12 19 25 32 41 45 54;
    11 20 24 33 40 46 53 55;
    21 23 34 39 47 52 56 61;
    22 35 38 48 51 57 60 62;
    36 37 49 50 58 59 63 64];
zz=zeros(64,3);
for i=1:3
    quant_i=quant(:,:,i);
    zz(Zigzag(:),i)=quant_i(:);
end

end

function coeffs = DeZigZag8x8(zz)
%  Function Name : DeZigZag8x8.m
%  Input         : zz    (Coefficients in zig-zag order)
%
%  Output        : coeffs(DCT coefficients in original order)

Zigzag = [1 2 6 7 15 16 28 29;
    3 5 8 14 17 27 30 43;
    4 9 13 18 26 31 42 44;
    10 12 19 25 32 41 45 54;
    11 20 24 33 40 46 53 55;
    21 23 34 39 47 52 56 61;
    22 35 38 48 51 57 60 62;
    36 37 49 50 58 59 63 64];
coeffs=zeros(8,8,3);
for i=1:3
    coeffs_temp=zz(Zigzag(:),i);
    coeffs(:,:,i)=reshape(coeffs_temp,8,8);
end
    
end

function zze = ZeroRunEnc_EoB(zz, EOB)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)

zze=[];
count=0;
i=1;
zz=[zz,1];
n=length(zz);
while i<=n
    if zz(i)~=0
        zze=[zze,zz(i)];
        i=i+1;
    else
        if mod(i,64)==0
            zze=[zze,EOB];
        else
            temp=[0,count];
            while zz(i+1)==0
                count=count+1;
                temp=[0,count];
                i=i+1;
                if mod(i,64)==0
                    temp=EOB;
                    break;
                end
            end
            zze=[zze,temp];
        end
        count=0;
        i=i+1;
    end
end
zze=zze(1:end-1);
end

function dst_eob = ZeroRunDec_EoB(src, EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)

dst_eob=[];
flag=0;
for i=1:length(src)
    if flag==0
        if src(i)~=0
            dst_eob=[dst_eob,src(i)];
        else
            dst_eob=[dst_eob,src(i)];
            flag=1;
        end
    else
        dst_eob=[dst_eob,zeros(1,src(i))];
        flag=0;
    end
end


index_EoB=find(dst_eob==EoB);
for i=1:length(index_EoB)
    index_EoB=find(dst_eob==EoB);
    num=index_EoB(1)-1;
    num_zeros_in_eob=64-mod(num,64);
    dst_eob=[dst_eob(1:index_EoB(1)-1),zeros(1,num_zeros_in_eob),dst_eob(index_EoB(1)+1:end)];
end
end

function [bytestream] = enc_huffman_new( data, BinCode, Codelengths)

a = BinCode(data(:),:)';
b = a(:);
mat = zeros(ceil(length(b)/8)*8,1);
p  = 1;
for i = 1:length(b)
    if b(i)~=' '
        mat(p,1) = b(i)-48;
        p = p+1;
    end
end
p = p-1;
mat = mat(1:ceil(p/8)*8);
d = reshape(mat,8,ceil(p/8))';
multi = [1 2 4 8 16 32 64 128];
bytestream = sum(d.*repmat(multi,size(d,1),1),2);

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

