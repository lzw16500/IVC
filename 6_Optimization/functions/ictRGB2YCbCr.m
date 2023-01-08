function yuv = ictRGB2YCbCr(rgb)
% Input         : rgb (Original RGB Image)
% Output        : yuv (YCbCr image after transformation)
% YOUR CODE HERE
r=rgb(:,:,1);
g=rgb(:,:,2);
b=rgb(:,:,3);
y=0.299*r+0.587*g+0.114*b;
cb=-0.169*r-0.331*g+0.5*b;
cr=0.5*r-0.419*g-0.081*b;
yuv(:,:,1)=y;
yuv(:,:,2)=cb;
yuv(:,:,3)=cr;
end