function zze = ZeroRunEnc_EoB(zz, EOB)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)

zze=zeros(size(zz));
count=1;
for i=1:64:length(zz)
    count_following_zero=0;
    flag=0; % flag=0 if previous one is non-zero or next one is non-zero
    zz_block=zz(i:i+63);
    for j=1:64
        if zz_block(j)~=0 && flag==0
            zze(count)=zz_block(j);
            count=count+1;
            flag=0;
        elseif zz_block(j)==0 && flag==0
            while(j+count_following_zero+1<=64 && zz_block(j+count_following_zero+1)==0)
                count_following_zero=count_following_zero+1;
            end
            flag=1;
            zze(count:count+1)=[0,count_following_zero];
            count=count+2;
        end
        if j<=63 && zz_block(j)==0 && zz_block(j+1)~=0
            count_following_zero=0;
            flag=0;
        end
    end
    if count==3
        zze(1)=EOB;
        count=2;
    elseif zze(count-3)~=0 && zze(count-2)==0
        zze=zze(1:count-3);
        zze(count-2)=EOB;
        count=count-1;
    end
end
    eob_index=find(zze~=0);
    zze=zze(1:max(eob_index));
end