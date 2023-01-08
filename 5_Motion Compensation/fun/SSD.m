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