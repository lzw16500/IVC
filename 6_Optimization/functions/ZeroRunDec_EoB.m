function dst = ZeroRunDec_EoB(src, EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)
dst=[];
for i=1:length(src)
    if src(i)==EoB
        length_dst=length(dst);
        num_zero=64-mod(length_dst,64);
        dst(length_dst+1:length_dst+num_zero)=0;
    else 
        if src(i)~=0 
            if  i==1
                length_dst=length(dst);
                dst(length_dst+1)=src(i);
            elseif (i>=2 && src(i-1)~=0) || (i>=3 && src(i-1)==0 && src(i-2)==0)
                length_dst=length(dst);
                dst(length_dst+1)=src(i);
            end
        else
            if  i==1
                length_dst=length(dst);
                dst(length_dst+1:length_dst+1+src(i+1))=0;
            elseif i>=2 && src(i-1)~=0
                length_dst=length(dst);
                dst(length_dst+1:length_dst+1+src(i+1))=0;
            end    
        end
    end
end                
end