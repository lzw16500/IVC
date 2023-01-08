function dst_eob = ZeroRunDec_EoB(src, EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)

dst_eob=[];
flag=0;
for i=1:length(src)
    if flag==0
        if src(i)~=0
            dst_eob=[dst_eob,src(i)];
        else
            dst_eob=[dst_eob,src(i)];
            flag=1;
        end
    else
        dst_eob=[dst_eob,zeros(1,src(i))];
        flag=0;
    end
end


index_EoB=find(dst_eob==EoB);
for i=1:length(index_EoB)
    index_EoB=find(dst_eob==EoB);
    num=index_EoB(1)-1;
    num_zeros_in_eob=64-mod(num,64);
    dst_eob=[dst_eob(1:index_EoB(1)-1),zeros(1,num_zeros_in_eob),dst_eob(index_EoB(1)+1:end)];
end
end