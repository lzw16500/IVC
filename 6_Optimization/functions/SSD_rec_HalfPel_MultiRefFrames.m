function rec_image = SSD_rec_HalfPel_MultiRefFrames(ref_image, motion_vectors,ref_frame_indices)
%  Input         : ref_image(Reference Image, YCbCr image), cell type
%                  motion_vectors (indices matrix)
%                  ref_frame_indices
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

%%%%%%%%%%%%% reference image interpolation %%%%%%%%%%%%%%%%%
[x,y,~]=size(ref_image{1});
[X,Y]=meshgrid(1:y,1:x);
[Xq,Yq] = meshgrid(1:0.5:y,1:0.5:x);
ref_image_interp_group=cell(1,numel(ref_image));
for j=1:numel(ref_image)
    ref_image_temp=ref_image{j};
    for i=1:3
        ref_image_double=double(ref_image_temp(:,:,i));
        ref_image_interp(:,:,i)=interp2(X,Y,ref_image_double,Xq,Yq); % double,[0,255]!!
    end
    ref_image_interp_group{j}=ref_image_interp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rec_image=zeros(size(ref_image{1}));
[m,n,c]=size(ref_image{1});

for j=1:8:m  % row
    for k=1:8:n  % column
        block_no_row=(j+7)/8;
        block_no_column=(k+7)/8;
        ref_image_interp=ref_image_interp_group{ref_frame_indices(block_no_row,block_no_column)};
        dy=ceil(motion_vectors(block_no_row,block_no_column)/17)-9;
        dx=mod(motion_vectors(block_no_row,block_no_column),17)-9;
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




