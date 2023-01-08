function [y,cb,cr]=YCbCr_downsample(yuv_im,p,q,filter_length)
yuv_wrap=padarray(yuv_im,[4 4],'symmetric','both');
y=yuv_im(:,:,1);
for i=2:3
    yuv_wrap_samp1(:,:,i)=resample(yuv_wrap(:,:,i),p,q,filter_length)';
    yuv_wrap_samp2(:,:,i)=resample(yuv_wrap_samp1(:,:,i),p,q,filter_length)';
end
cb=yuv_wrap_samp2(3:end-2,3:end-2,2);
cr=yuv_wrap_samp2(3:end-2,3:end-2,3);
end