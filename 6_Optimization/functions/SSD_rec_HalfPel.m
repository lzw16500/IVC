function rec_image = SSD_rec_HalfPel(ref_image, motion_vectors)
%  Input         : ref_image(Reference Image, YCbCr image)
%                  motion_vectors (indices matrix)
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

%%%%%%%%%%%%% reference image interpolation %%%%%%%%%%%%%%%%%
[x,y,~]=size(ref_image);
[X,Y]=meshgrid(1:y,1:x);
[Xq,Yq] = meshgrid(1:0.5:y,1:0.5:x);
for i=1:3
    ref_image_double=double(ref_image(:,:,i));
    ref_image_interp(:,:,i)=interp2(X,Y,ref_image_double,Xq,Yq); % double,[0,255]!!
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rec_image=zeros(size(ref_image));
[m,n,c]=size(ref_image);

for j=1:8:m  % row
    for k=1:8:n  % column
        dy=ceil(motion_vectors((j+7)/8,(k+7)/8)/17)-9;
        dx=mod(motion_vectors((j+7)/8,(k+7)/8),17)-9;
        if dx==-9
            dx=8;
        end
        first_pixel_row=2*j-1+dy;
        first_pixel_column=2*k-1+dx;
        for i=1:c
            rec_image(j:j+7,k:k+7,i)=ref_image_interp(first_pixel_row:2:first_pixel_row+15,first_pixel_column:2:first_pixel_column+15,i);
        end
    end
end
end

        
        
        
