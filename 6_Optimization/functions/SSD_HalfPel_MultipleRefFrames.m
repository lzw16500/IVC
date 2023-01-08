function [motion_vectors_indices,ref_frame_indices] = SSD_HalfPel_MultipleRefFrames(ref_image, image)
%  Input         : ref_image(Reference Images, cell type, each cell with size: height*width*channel), ycbcr
%                  image (Current Image, size: height*width*channel), ycbcr
%
%  Output        : motion_vectors_indices (Motion Vector Indices, size: (height/8) x (width/8) x 1 )
%                 ref_frame_indices: matrix with the same size as "motion_vectors_indice", e.g.1: the farthest previous frame

index_mat=zeros(17,17);
temp=1;
for i=1:17
    for j=1:17
        index_mat(i,j)=temp;
        temp=temp+1;
    end
end

[m,n,~]=size(image);
motion_vectors_indices=zeros(m/8,n/8);
ref_frame_indices=zeros(m/8,n/8);

%%%%%%%%%%%%% reference images interpolation %%%%%%%%%%%%%%%%%
ref_image_interp_group=cell(1,numel(ref_image)); 
for i=1:numel(ref_image)
    temp_ref_image=ref_image{i};
    [x,y,~]=size(temp_ref_image);
    [X,Y]=meshgrid(1:y,1:x);
    [Xq,Yq] = meshgrid(1:0.5:y,1:0.5:x);
    ref_image_double_channel=double(temp_ref_image(:,:,1)); % 2 dimensional, luminance channel only
    ref_image_interp_group{i}=interp2(X,Y,ref_image_double_channel,Xq,Yq); % double,[0,255]!!
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each block: first full-pel full search across all reference images,
% then find the best full-pel search in each reference image, then hal-pel
% neighboring search in each reference image

% [ref_interp_row,ref_interp_column]=size(ref_image_interp);
% ref_ssd_min=zeros(1,numel(ref_image));

for i=1:8:m  % row
    for j=1:8:n  % column
        ssd_min_global=255^255;
        for k=1:numel(ref_image)
            ref_image_interp=ref_image_interp_group{k};
            
            start_dx=-8;
            start_dy=-8;
            end_dx=8;
            end_dy=8;
            
            current_block=image(i:i+7,j:j+7,1);
            ssd_min=255^255;
            optimal_dx_integer=9;
            optimal_dy_integer=9;
            if j==1 start_dx=0; end
            if i==1 start_dy=0; end
            if j==n-7 end_dx=0; end
            if i==m-7 end_dy=0; end
            
            
            for dx=start_dx:2:end_dx
                for dy=start_dy:2:end_dy
                    ref_block=ref_image_interp(2*i-1+dy:2:2*i-1+dy+15,2*j-1+dx:2:2*j-1+dx+15);
                    ssd_new=sum(sum((current_block-ref_block).^2));
                    if ssd_new<ssd_min
                        ssd_min=ssd_new;
%                         if k==1
%                         ref_frame_indices((i+7)/8,(j+7)/8)=k;
%                         motion_vectors_indices((i+7)/8,(j+7)/8)=index_mat(9+dy,9+dx);
%                         end
                        optimal_dx_integer=dx;
                        optimal_dy_integer=dy;
                    end
                end
            end % having found the motion vector for a single block in full-pel case in each reference image
            
            % refine in the neighborhood (half-pel)
            refine_start_dx=-1;
            refine_start_dy=-1;
            refine_end_dx=1;
            refine_end_dy=1;
            if j==1 refine_start_dx=0; end
            if i==1 refine_start_dy=0; end
            if j==n-7 refine_end_dx=0; end
            if i==m-7 refine_end_dy=0; end
            if optimal_dx_integer==-8 refine_start_dx=0; end
            if optimal_dy_integer==-8 refine_start_dy=0; end
            if optimal_dx_integer==8 refine_end_dx=0; end
            if optimal_dy_integer==8 refine_end_dy=0; end
            for dx=optimal_dx_integer+refine_start_dx:optimal_dx_integer+refine_end_dx
                for dy=optimal_dy_integer+refine_start_dy:optimal_dy_integer+refine_end_dy
%                     if dx~=optimal_dx_integer || dy~=optimal_dy_integer
                        ref_block_halfPel=ref_image_interp(2*i-1+dy:2:2*i-1+dy+15,2*j-1+dx:2:2*j-1+dx+15);
                        ssd_new=sum(sum((current_block-ref_block_halfPel).^2));
                        if ssd_new<ssd_min_global
                            ssd_min_global=ssd_new;
                            ref_frame_indices((i+7)/8,(j+7)/8)=k;
                            motion_vectors_indices((i+7)/8,(j+7)/8)=index_mat(9+dy,9+dx);
                        end
%                     end
                end
            end

            
        end
    end
end

end