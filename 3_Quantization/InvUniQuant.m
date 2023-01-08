function image = InvUniQuant(qImage, bits)
%  Input         : qImage (Quantized Image)
%                : bits (bits available for representatives)
%
%  Output        : image (Mid-rise de-quantized Image)

norm_image=(qImage+0.5)/(2^bits);
image=floor(norm_image*256);

end