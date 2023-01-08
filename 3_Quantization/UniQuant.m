function qImage = UniQuant(image, bits)
%  Input         : image (Original Image)
%                : bits (bits available for representatives)
%
%  Output        : qImage (Quantized Image)
norm_image=image/256;
qImage=floor(norm_image*(2^bits));
end