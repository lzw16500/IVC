function MSE = calcMSE(Image, recImage)
% Input         : Image    (Original Image)
%                 recImage (Reconstructed Image)
% Output        : MSE      (Mean Squared Error)
% YOUR CODE HERE

[h,w,c]=size(Image);
Image=double(Image);
recImage=double(recImage);
diff=Image(:)-recImage(:);
MSE=sum((diff).^2)/(h*w*c);
end