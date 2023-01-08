clc;clear;

image=imread('lena_small.tif');


reshape_image=[];
[m,n,c]=size(image);
for i=1:c
    for j=1:2:n
        for k=1:2:m
            reshape_image=[reshape_image; image(k,j,i) image(k+1,j,i) image(k,j+1,i) image(k+1,j+1,i)]; % reshape the image
        end
    end
end


[w,h,c]=size(image);
for i=1:c
    Cell(:,:,i)=mat2cell(image(:,:,i),repmat(2,[1 w/2]),repmat(2,[1 h/2]));
end
numcell=numel(Cell);
reshaped_block=zeros(numcell,2^2);
for q=1:numel(Cell)
    A=cell2mat(Cell(q));
    reshaped_block(q,:)=reshape(A,1,[]);
end