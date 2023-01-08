%% submission chapter 5
% Zhiwei Liang, 03734179

% please add foreman sequence and lena_small.tif to path

%% initialization
clc;clear;close all;

tic 
lena_small = double(imread('lena_small.tif'));
lena_small_ycbcr=ictRGB2YCbCr(lena_small);

% load image sequence
sequence=cell(1,21);
for i=1:21
    file_name=['foreman00',num2str(i+19),'.bmp'];
    sequence{i}=double(imread(file_name));
end

scales = [0.07 0.2 0.4 0.8 1.0 1.5 2 3 4 4.5]; % quantization scale factor, the larger, the lower PSNR
% scales = [0.07]; % quantization scale factor, the larger, the lower PSNR

%% video Codec
fprintf('starting video Codec\n');
rate_mat_video=zeros(21,numel(scales));
psnr_mat_video=zeros(21,numel(scales));
for scaleIdx=1:numel(scales)
    rec_frame_ycbcr=cell(1,21);
    %% first frame
    qScale   = scales(scaleIdx);
    k_small  = IntraEncode(lena_small, qScale);
    k        = IntraEncode(sequence{1}, qScale);
    
    % train Huffman table
    pmf_ind=hist(k_small(:),min(k):max(k));
    pmf_ind=pmf_ind/sum(pmf_ind);
    [BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind);
    
    % encode the first image
    bytestream = enc_huffman_new(k-min(k)+1, BinCode, Codelengths);
    bitPerPixel = (numel(bytestream)*8) / (numel(sequence{1})/3);
    
    % decode the first image
    k_rec=double(dec_huffman_new(bytestream,BinaryTree,size(k(:),1)))+min(k)-1;
    I_rec = IntraDecode(k_rec, size(sequence{1}),qScale);  % rgb
    rec_frame_ycbcr{1}=ictRGB2YCbCr(I_rec);
    
    % calcualte bitrate and PSNR
    rate_mat_video(1,scaleIdx)=bitPerPixel;
    psnr_mat_video(1,scaleIdx)=calcPSNR(sequence{1},I_rec);
    
    %% following frames
    for i=2:21
        ref_frame_ycbcr=rec_frame_ycbcr{i-1};
        current_frame_rgb=sequence{i};
        current_frame_ycbcr=ictRGB2YCbCr(current_frame_rgb);
        motion_vector=SSD(ref_frame_ycbcr(:,:,1),current_frame_ycbcr(:,:,1));
        motion_compensation_image_ycbcr=SSD_rec(ref_frame_ycbcr,motion_vector);
        residual_image_ycbcr=current_frame_ycbcr-motion_compensation_image_ycbcr;
        residual_image_rgb=ictYCbCr2RGB(residual_image_ycbcr);
        zeroRun_residual=IntraEncode(residual_image_rgb,qScale);  % residual image
        
        % train Huffman code
        if i==2
            % Huffman code for zero run
%             pmf_ind_zerorun=hist(zeroRun(:),min(zeroRun):max(zeroRun));
            pmf_ind_zerorun=hist(zeroRun_residual(:),-1000:4000);  % why -1000:4000?
            pmf_ind_zerorun=pmf_ind_zerorun/sum(pmf_ind_zerorun);
            [BinaryTree_zerorun,HuffCode_zerorun,BinCode_zerorun,Codelengths_zerorun]=buildHuffman(pmf_ind_zerorun);
            
            % Huffman code for motion vectors
            pmf_ind_mv=hist(motion_vector(:),1:81);
            pmf_ind_mv=pmf_ind_mv/sum(pmf_ind_mv);
            [BinaryTree_mv,HuffCode_mv,BinCode_mv,Codelengths_mv]=buildHuffman(pmf_ind_mv);
        end
        
        % encode Huffman
%         bytestream1=enc_huffman_new(zeroRun_residual-min(zeroRun_residual)+1, BinCode_zerorun, Codelengths_zerorun); % rate for zero run
        bytestream1=enc_huffman_new(zeroRun_residual+1001, BinCode_zerorun, Codelengths_zerorun); % rate for zero run  % Why 1001?
        bytestream2=enc_huffman_new(motion_vector(:)-min(motion_vector(:))+1,BinCode_mv,Codelengths_mv); % rate for motion vector
        bitPerPixel1 = (numel(bytestream1)*8) / (numel(sequence{i})/3);
        bitPerPixel2 = (numel(bytestream2)*8) / (numel(sequence{i})/3);
        rate_mat_video(i,scaleIdx)=bitPerPixel1+bitPerPixel2;
        
        % decode Huffman
        image_size=[size(motion_vector,1)*8,size(motion_vector,2)*8,3];
        decoded_frame_rgb=IntraDecode(zeroRun_residual,image_size,qScale)+ictYCbCr2RGB(motion_compensation_image_ycbcr);
        rec_frame_ycbcr{i}=ictRGB2YCbCr(decoded_frame_rgb);
        psnr_mat_video(i,scaleIdx)=calcPSNR(sequence{i},decoded_frame_rgb);
    end
    fprintf('Quantization scale: %.1f, bitrate: %.2f bits/pixel, PSNR: %.2fdB\n',qScale,mean(rate_mat_video(:,scaleIdx)),mean(psnr_mat_video(:,scaleIdx)));
end

rate_video=mean(rate_mat_video,1);
psnr_video=mean(psnr_mat_video,1);

% plot
figure;
plot(rate_video,psnr_video,'b-*');hold on;
xlabel('bit/pixel');
ylabel('PSNR [dB]');

%% still image Codec
scales_still=[0.15,0.3,0.7,1.0,1.5,3,5,7,10];
fprintf('starting still image Codec\n');
rate_mat_still=zeros(21,numel(scales_still));
psnr_mat_still=zeros(21,numel(scales_still));

for scaleIdx=1:numel(scales_still)
    
    % train Huffman table
    qScale  = scales_still(scaleIdx);
    k_base  = IntraEncode(lena_small, qScale);
    pmf_ind=hist(k_base(:),-1000:4000);
    pmf_ind=pmf_ind/sum(pmf_ind);
    [BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind);
    
    for i=1:21
    still_im_zerorun  = IntraEncode(sequence{i}, qScale);
    
    % encode the image
    bytestream = enc_huffman_new(still_im_zerorun+1001, BinCode, Codelengths);
    bitPerPixel = (numel(bytestream)*8) / (numel(sequence{i})/3);
    
    % decode the image
    k_rec=double(dec_huffman_new(bytestream,BinaryTree,size(still_im_zerorun(:),1)))-1001;
    I_rec = IntraDecode(k_rec, size(sequence{i}),qScale);  % rgb
%     I_rec=ictRGB2YCbCr(I_rec);
    
    % calcualte bitrate and PSNR
    rate_mat_still(i,scaleIdx)=bitPerPixel;
    psnr_mat_still(i,scaleIdx)=calcPSNR(sequence{i},I_rec);
    end
    fprintf('Quantization scale: %.1f, bitrate: %.2f bits/pixel, PSNR: %.2fdB\n',qScale,mean(rate_mat_still(:,scaleIdx)),mean(psnr_mat_still(:,scaleIdx)));
end
rate_still=mean(rate_mat_still,1);
psnr_still=mean(psnr_mat_still,1);
plot(rate_still,psnr_still,'r-*');hold on;
axis([0 4 20 45])
legend('Video Codec','Still image Codec');
toc

%% put all used sub-functions here.
function motion_vectors_indices = SSD(ref_image, image)
%  Input         : ref_image(Reference Image, size: height x width)
%                  image (Current Image, size: height x width)
%
%  Output        : motion_vectors_indices (Motion Vector Indices, size: (height/8) x (width/8) x 1 )

index_mat=zeros(9,9);
temp=1;
for i=1:9
    for j=1:9
        index_mat(i,j)=temp;
        temp=temp+1;
    end
end

[m,n]=size(image);
motion_vectors_indices=zeros(m/8,n/8);


for i=1:8:m
    for j=1:8:n
        start_dx=-4;
        start_dy=-4;
        end_dx=4;
        end_dy=4;
        
        current_block=image(i:i+7,j:j+7);
        ssd=255^255;
        if j==1 start_dx=0; end
        if i==1 start_dy=0; end
        if j==n-7 end_dx=0; end
        if i==m-7 end_dy=0; end
        for dx=start_dx:end_dx
            for dy=start_dy:end_dy
                ref_block=ref_image(i+dy:i+dy+7,j+dx:j+dx+7);
                ssd_new=sum(sum((current_block-ref_block).^2));
                if ssd_new<ssd
                    ssd=ssd_new;
                    motion_vectors_indices((i+7)/8,(j+7)/8)=index_mat(5+dy,5+dx);
                end
            end
        end
    end
end
end

function rec_image = SSD_rec(ref_image, motion_vectors)
%  Input         : ref_image(Reference Image, YCbCr image)
%                  motion_vectors
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

rec_image=zeros(size(ref_image));
[m,n,c]=size(ref_image);

for j=1:8:m  % row
    for k=1:8:n  % column
        dy=ceil(motion_vectors((j+7)/8,(k+7)/8)/9)-5;
        dx=mod(motion_vectors((j+7)/8,(k+7)/8),9)-5;
        if dx==-5
            dx=4;
        end
        first_pixel_row=j+dy;
        first_pixel_column=k+dx;
        for i=1:c
            rec_image(j:j+7,k:k+7,i)=ref_image(first_pixel_row:first_pixel_row+7,first_pixel_column:first_pixel_column+7,i);
        end
    end
end
end

function dst = IntraDecode(image, img_size , qScale)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image),rgb
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

zze=zeros(size(zz));
count=1;
for i=1:64:length(zz)
    count_following_zero=0;
    flag=0; % flag=0 if previous one is non-zero or next one is non-zero
    zz_block=zz(i:i+63);
    for j=1:64
        if zz_block(j)~=0 && flag==0
            zze(count)=zz_block(j);
            count=count+1;
            flag=0;
        elseif zz_block(j)==0 && flag==0
            while(j+count_following_zero+1<=64 && zz_block(j+count_following_zero+1)==0)
                count_following_zero=count_following_zero+1;
            end
            flag=1;
            zze(count:count+1)=[0,count_following_zero];
            count=count+2;
        end
        if j<=63 && zz_block(j)==0 && zz_block(j+1)~=0
            count_following_zero=0;
            flag=0;
        end
    end
    if count==3
        zze(1)=EOB;
        count=2;
    elseif zze(count-3)~=0 && zze(count-2)==0
        zze=zze(1:count-3);
        zze(count-2)=EOB;
        count=count-1;
    end
end
    eob_index=find(zze~=0);
    zze=zze(1:max(eob_index));
end

function dst_eob = ZeroRunDec_EoB(src, EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)

dst_eob=[];
count_dst=1;
count_src=1;
while (count_src<=length(src))
    if src(count_src)~=0 
        dst_eob(count_dst)=src(count_src);
        count_dst=count_dst+1;
        count_src=count_src+1;
    elseif src(count_src)==0
        num_following_zero=src(count_src+1);
        zero_array=zeros(1,num_following_zero+1);
        dst_eob(count_dst:count_dst+num_following_zero)=zero_array;
        count_dst=count_dst+num_following_zero+1;
        count_src=count_src+2;
    end

end

while (ismember(EoB,dst_eob))
    index_eob=find(dst_eob==EoB);
    array_length_before_eob=length(dst_eob(1:index_eob(1)-1));
    num_zeros_in_eob=64-mod(array_length_before_eob,64);
    temp=dst_eob(index_eob(1)+1:end);
    dst_eob(index_eob(1)+num_zeros_in_eob:length(dst_eob)+num_zeros_in_eob-1)=temp;
    dst_eob(index_eob(1):index_eob(1)+num_zeros_in_eob-1)=zeros(1,num_zeros_in_eob);


end

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

function [output] = dec_huffman_new (bytestream, BinaryTree, nr_symbols)
output = zeros(1,nr_symbols);
ctemp = BinaryTree;
dec = zeros(size(bytestream,1),8);
for i = 8:-1:1
    dec(:,i) = rem(bytestream,2);
    bytestream = floor(bytestream/2);
end
dec = dec(:,end:-1:1)';
a = dec(:);
i = 1;
p = 1;
while(i <= nr_symbols)&&p<=max(size(a))
    while(isa(ctemp,'cell'))
        next = a(p)+1;
        p = p+1;
        ctemp = ctemp{next};
    end;
    output(i) = ctemp;
    ctemp = BinaryTree;
    i=i+1;
end;
end
%%
function [ BinaryTree, HuffCode, BinCode, Codelengths] = buildHuffman( p );
global y
p=p(:)/sum(p)+eps;              % normalize histogram
p1=p;                           % working copy
c=cell(length(p1),1);			% generate cell structure 
for i=1:length(p1)				% initialize structure
   c{i}=i;						
end
while size(c)-2					% build Huffman tree
	[p1,i]=sort(p1);			% Sort probabilities
	c=c(i);						% Reorder tree.
	c{2}={c{1},c{2}};           % merge branch 1 to 2
    c(1)=[];	                % omit 1
	p1(2)=p1(1)+p1(2);          % merge Probabilities 1 and 2 
    p1(1)=[];	                % remove 1
end
getcodes(c,[]);                  % recurse to find codes
code=char(y);
[numCodes maxlength] = size(code); % get maximum codeword length
length_b=0;
HuffCode=zeros(1,numCodes);
for symbol=1:numCodes
    for bit=1:maxlength
        length_b=bit;
        if(code(symbol,bit)==char(49)) HuffCode(symbol) = HuffCode(symbol)+2^(bit-1)*(double(code(symbol,bit))-48);
        elseif(code(symbol,bit)==char(48))
        else 
            length_b=bit-1;
            break;
        end;
    end;
    Codelengths(symbol)=length_b;
end;
BinaryTree = c;
BinCode = code;
clear global y;
return
end
%----------------------------------------------------------------
function getcodes(a,dum)       
global y                            % in every level: use the same y
if isa(a,'cell')                    % if there are more branches...go on
         getcodes(a{1},[dum 0]);    % 
         getcodes(a{2},[dum 1]);
else   
   y{a}=char(48+dum);   
end
end
