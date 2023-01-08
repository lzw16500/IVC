function dst = IntraEncode(image, qScale)
%  Function Name : IntraEncode.m
%  Input         : image (Original RGB Image)
%                  qScale(quantization scale)
%  Output        : dst   (sequences after zero-run encoding, 1xN)
image=ictRGB2YCbCr(image);
image_dct=blockproc(image,[8,8],@(block_struct) DCT8x8(block_struct.data));
image_dct_quantized=blockproc(image_dct,[8,8],@(block_struct) Quant8x8(block_struct.data,qScale));
zig_zag=blockproc(image_dct_quantized,[8 8],@(block_struct) ZigZag8x8(block_struct.data)); % zig_zag is still a matrix, with each block of size 64*3
eob=1000;
zz=[];
for i=1:64:size(zig_zag,1)
    for j=1:3:size(zig_zag,2)
        zigzag_block=zig_zag(i:i+63,j:j+2);
        zigzag_vec=zigzag_block(:)';
        zz=[zz,zigzag_vec];
    end
end
dst=ZeroRunEnc_EoB(zz,eob);
end