function zze = ZeroRunEnc_EoB(zz, EOB)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)

zze=[];
count=0;
i=1;
zz=[zz,1];
n=length(zz);
while i<=n
    if zz(i)~=0
        zze=[zze,zz(i)];
        i=i+1;
    else
        if mod(i,64)==0
            zze=[zze,EOB];
        else
            temp=[0,count];
            while zz(i+1)==0
                count=count+1;
                temp=[0,count];
                i=i+1;
                if mod(i,64)==0
                    temp=EOB;
                    break;
                end
            end
            zze=[zze,temp];
        end
        count=0;
        i=i+1;
    end
end
zze=zze(1:end-1);
end
