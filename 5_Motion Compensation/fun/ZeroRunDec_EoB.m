function dst_eob = ZeroRunDec_EoB(src, EoB)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)

dst_eob=[];
count_dst=1;
count_src=1;
while (count_src<=length(src))
    if src(count_src)~=0 
        dst_eob(count_dst)=src(count_src);
        count_dst=count_dst+1;
        count_src=count_src+1;
    elseif src(count_src)==0
        num_following_zero=src(count_src+1);
        zero_array=zeros(1,num_following_zero+1);
        dst_eob(count_dst:count_dst+num_following_zero)=zero_array;
        count_dst=count_dst+num_following_zero+1;
        count_src=count_src+2;
    end

end

while (ismember(EoB,dst_eob))
    index_eob=find(dst_eob==EoB);
    array_length_before_eob=length(dst_eob(1:index_eob(1)-1));
    num_zeros_in_eob=64-mod(array_length_before_eob,64);
    temp=dst_eob(index_eob(1)+1:end);
    dst_eob(index_eob(1)+num_zeros_in_eob:length(dst_eob)+num_zeros_in_eob-1)=temp;
    dst_eob(index_eob(1):index_eob(1)+num_zeros_in_eob-1)=zeros(1,num_zeros_in_eob);


end

end