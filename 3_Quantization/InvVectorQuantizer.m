function image = InvVectorQuantizer(qImage, clusters, block_size)
%  Function Name : VectorQuantizer.m
%  Input         : qImage     (Quantized Image)
%                  clusters   (Quantization clusters)
%                  block_size (Block Size)
%  Output        : image      (Dequantized Images)

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