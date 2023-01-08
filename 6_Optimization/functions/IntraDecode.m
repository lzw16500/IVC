function dst = IntraDecode(image, img_size , qScale)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image),rgb
m=img_size(1);
n=img_size(2);
ch=img_size(3);
dst_temp=zeros(m,n,ch);
eob=1000;
zz=ZeroRunDec_EoB(image,eob);
count=1;
for a=1:8:m
    for b=1:8:n
        temp_zz=reshape(zz(count:count+191),64,3);
        dst_temp(a:a+7,b:b+7,:)=DeZigZag8x8(temp_zz);
        count=count+192;
    end
end
dct_coeff=blockproc(dst_temp,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale));
dst=blockproc(dct_coeff,[8,8],@(block_struct) IDCT8x8(block_struct.data));
dst = ictYCbCr2RGB(dst);
end