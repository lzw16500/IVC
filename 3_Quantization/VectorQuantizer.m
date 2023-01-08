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