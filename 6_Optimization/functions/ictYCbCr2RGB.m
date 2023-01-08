function rgb = ictYCbCr2RGB(yuv)
% Input         : yuv (Original YCbCr image)
% Output        : rgb (RGB Image after transformation)
% YOUR CODE HERE
y=yuv(:,:,1);
cb=yuv(:,:,2);
cr=yuv(:,:,3);
rgb(:,:,1)=y+1.402*cr;
rgb(:,:,2)=y-0.344*cb-0.714*cr;
rgb(:,:,3)=y+1.772*cb;
end