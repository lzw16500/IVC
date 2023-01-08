function upsampled_yuv=YCbCr_upsample(y,cb,cr,p,q,filter_length)
cb_wrap=padarray(cb,[2 2],'symmetric','both');
cr_wrap=padarray(cr,[2 2],'symmetric','both');

cb_wrap_samp1=resample(cb_wrap,p,q,filter_length)';
cb_wrap_samp2=resample(cb_wrap_samp1,p,q,filter_length)';

cr_wrap_samp1=resample(cr_wrap,p,q,filter_length)';
cr_wrap_samp2=resample(cr_wrap_samp1,p,q,filter_length)';

upsampled_yuv(:,:,1)=y;
upsampled_yuv(:,:,2)=cb_wrap_samp2(5:end-4,5:end-4);
upsampled_yuv(:,:,3)=cr_wrap_samp2(5:end-4,5:end-4);
end