function [qImage, clusters] = LloydMax(image, bits, epsilon)
%  Input         : image (Original RGB Image)
%                  bits (bits for quantization)
%                  epsilon (Stop Condition)
%  Output        : qImage (Quantized Image)
%                  clusters (Quantization Table)
num_rep_values=2^bits;
interval_length=256/(2^bits);
clusters=zeros(num_rep_values,3);
for i=1:num_rep_values
    clusters(i,1)=i*interval_length-interval_length/2;
end
vector_image=image(:);
mse_old=10000;  % initialize mse to be very large

while 1
    [distance,index]=pdist2(clusters(:,1),vector_image,'euclidean','Smallest',1); % index starts from 1
    mse_new=sum(distance.^2)/numel(vector_image);
    for i=1:numel(vector_image)  % set up the quantization table
        clusters(index(i),2)=clusters(index(i),2)+vector_image(i);
        clusters(index(i),3)=clusters(index(i),3)+1;
    end
    
    for i=1:num_rep_values
        if clusters(i,3)~=0
            clusters(i,1)=round(clusters(i,2)/clusters(i,3)); % update representative values. Representatives can be non-integer.
        end
    end
    for i=1:num_rep_values
        if clusters(i,3)==0
            [~,max_num]=max(clusters(:,3));
            clusters(i,1)=clusters(max_num,1)+1;
            clusters(i,3)=floor(clusters(max_num,3)/2);
            clusters(max_num,3)=clusters(max_num,3)-clusters(i,3);
        end
    end
    
    clusters(:,2:3)=0; % reset the table
    
    if abs(mse_old-mse_new)/mse_new<epsilon
        break;
    end
    mse_old=mse_new;
end

[~,index]=pdist2(clusters(:,1),vector_image,'euclidean','Smallest',1);
qImage=reshape(index,size(image));

end