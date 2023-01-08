function rec_image = SSD_rec(ref_image, motion_vectors)
%  Input         : ref_image(Reference Image, YCbCr image)
%                  motion_vectors
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)

rec_image=zeros(size(ref_image));
[m,n,c]=size(ref_image);

for j=1:8:m  % row
    for k=1:8:n  % column
        dy=ceil(motion_vectors((j+7)/8,(k+7)/8)/9)-5;
        dx=mod(motion_vectors((j+7)/8,(k+7)/8),9)-5;
        if dx==-5
            dx=4;
        end
        first_pixel_row=j+dy;
        first_pixel_column=k+dx;
        for i=1:c
            rec_image(j:j+7,k:k+7,i)=ref_image(first_pixel_row:first_pixel_row+7,first_pixel_column:first_pixel_column+7,i);
        end
    end
end
end
