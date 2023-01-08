function motion_vectors_indices = SSD_HalfPel_2(ref_image, image)
%  Input         : ref_image(Reference Image, size: height x width)
%                  image (Current Image, size: height x width)
%
%  Output        : motion_vectors_indices (Motion Vector Indices, size: (height/8) x (width/8) x 1 )

index_mat=zeros(17,17);
temp=1;
for i=1:17
    for j=1:17
        index_mat(i,j)=temp;
        temp=temp+1;
    end
end

[m,n]=size(image);
motion_vectors_indices=zeros(m/8,n/8);

%%%%%%%%%%%%% reference image interpolation %%%%%%%%%%%%%%%%%
ref_image_double=double(ref_image);
[x,y]=size(ref_image_double);
[X,Y]=meshgrid(1:y,1:x);
[Xq,Yq] = meshgrid(1:0.5:y,1:0.5:x);
ref_image_interp=interp2(X,Y,ref_image_double,Xq,Yq); % double,[0,255]!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(image);
motion_vectors_indices=zeros(m/8,n/8);


for i=1:8:m
    for j=1:8:n
        start_dx=-8;
        start_dy=-8;
        end_dx=8;
        end_dy=8;
        
        current_block=image(i:i+7,j:j+7);
        ssd=255^255;
        if j==1 start_dx=0; end
        if i==1 start_dy=0; end
        if j==n-7 end_dx=0; end
        if i==m-7 end_dy=0; end
        for dx=start_dx:end_dx
            for dy=start_dy:end_dy
                ref_block=ref_image_interp(2*i-1+dy:2:2*i-1+dy+15,2*j-1+dx:2:2*j-1+dx+15);
                ssd_new=sum(sum((current_block-ref_block).^2));
                if ssd_new<ssd
                    ssd=ssd_new;
                    motion_vectors_indices((i+7)/8,(j+7)/8)=index_mat(9+dy,9+dx);
                end
            end
        end % having found the motion vector for a single block in full-pel case
    end
end
end