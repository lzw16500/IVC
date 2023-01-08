function qImage = ApplyVectorQuantizer(image, clusters, bsize)
%  Function Name : ApplyVectorQuantizer.m
%  Input         : image    (Original Image)
%                  clusters (Quantization Representatives)
%                  bsize    (Block Size)
%  Output        : qImage   (Quantized Image)


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