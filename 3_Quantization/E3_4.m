clc;clear;

%% Main
bits         = 8;
epsilon      = 0.1;
block_size   = 2;
%% lena small for VQ training
image_small  = double(imread('lena_small.tif'));
[clusters, Temp_clusters] = VectorQuantizer(image_small, bits, epsilon, block_size);
qImage_small              = ApplyVectorQuantizer(image_small, clusters, block_size);
%% Huffman table training
pmf_ind=hist(qImage_small(:),1:256);
pmf_ind=pmf_ind/sum(pmf_ind);
[BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind); 
%% lena quantization
image  = double(imread('lena.tif'));
qImage = ApplyVectorQuantizer(image, clusters, block_size);
%% Huffman encoding
bytestream = enc_huffman_new(qImage(:), BinCode, Codelengths);
%% bit rate
bpp  = (numel(bytestream) * 8) / (numel(image)/3);
%% Huffman decoding
qReconst_image=reshape(dec_huffman_new (bytestream, BinaryTree, numel(qImage)), size(qImage));

%% lena reconstruction
reconst_image  = InvVectorQuantizer(qReconst_image, clusters, block_size);
PSNR = calcPSNR(image, reconst_image);

%% sub-functions
function [clusters, Temp_clusters]= VectorQuantizer(image, bits, epsilon, bsize)

[m,n,c]=size(image);
reshape_image=zeros(m*n*c/bsize^2,bsize^2);
t=1;
for i=1:c
    for j=1:bsize:n
        for k=1:bsize:m
            reshape_image(t,:)=[image(k,j,i) image(k+1,j,i) image(k,j+1,i) image(k+1,j+1,i)]; % reshape the image
            t=t+1;
        end
    end
end

Temp_clusters=zeros(2^bits,bsize^2);
for i=1:2^bits
    Temp_clusters(i,:)=floor((i-0.5)*(256/2^bits));
end
clusters=Temp_clusters;

num_im_vectors=size(reshape_image,1);
num_rep_vectors=size(Temp_clusters,1);
mse_old=10000;
sum_cell=zeros(num_rep_vectors,bsize^2);
count=zeros(num_rep_vectors,1);

while 1
    [index,distance]=knnsearch(Temp_clusters,reshape_image,'Distance','euclidean');
    mse_new=sum(distance.^2)/numel(index);
    for i=1:num_im_vectors
        sum_cell(index(i),:)=sum_cell(index(i),:)+reshape_image(i,:);
        count(index(i))=count(index(i))+1;
    end
    
    for i=1:num_rep_vectors
        if count(i)~=0
            clusters(i,:)=round(sum_cell(i,:)/count(i));
        end
    end
    for i=1:num_rep_vectors
        if count(i)==0
            [~,max_num]=max(count);
            clusters(i,:)=clusters(max_num,:)+[0 0 0 1];
            count(i)=floor(count(max_num)/2);
            count(max_num)=count(max_num)-count(i);
        end

    end
    if abs(mse_old-mse_new)/mse_new<epsilon
        break;
    end
    sum_cell=zeros(num_rep_vectors,bsize^2);
    count=zeros(num_rep_vectors,1);
    Temp_clusters=clusters;
    mse_old=mse_new;
end

end

function qImage = ApplyVectorQuantizer(image, clusters, bsize)
[m,n,c]=size(image);
reshape_image=zeros(m*n*c/bsize^2,bsize^2);
t=1;
for i=1:c
    for j=1:2:n
        for k=1:2:m
            reshape_image(t,:)=[image(k,j,i) image(k+1,j,i) image(k,j+1,i) image(k+1,j+1,i)]; % reshape the image
            t=t+1;
        end
    end
end
[index,~]=knnsearch(clusters,reshape_image,'Distance','euclidean');
qImage=reshape(index,[m/bsize,n/bsize,c]);
end

function image = InvVectorQuantizer(qImage, clusters, block_size)
[m,n,c]=size(qImage);
image=zeros(m*block_size,n*block_size,c);
reshape_image=zeros(numel(qImage),block_size^2);
vector_qImage=qImage(:);
for i=1:numel(qImage)
    reshape_image(i,:)=clusters(vector_qImage(i),:);
end

t=1;
for i=1:c
   for j=1:block_size:n*block_size
       for k=1:block_size:m*block_size
           image(k,j,i)=reshape_image(t,1);
           image(k+1,j,i)=reshape_image(t,2);
           image(k,j+1,i)=reshape_image(t,3);
           image(k+1,j+1,i)=reshape_image(t,4);
           t=t+1;
       end
   end
end
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