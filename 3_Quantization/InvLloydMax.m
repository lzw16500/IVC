function image = InvLloydMax(qImage, clusters)
%  Input         : qImage   (Quantized Image)
%                  clusters (Quantization Table)
%  Output        : image    (Recovered Image)

vector_qImage=qImage(:);
image=[];
for i=1:numel(vector_qImage)
    image=[image;clusters(qImage(i),1)];
end
image=reshape(image,size(qImage));

end