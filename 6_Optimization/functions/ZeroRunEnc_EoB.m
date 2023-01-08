function zze = ZeroRunEnc_EoB(zz, EoB)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)
zze=[];
num_zero=0;
for i=1:length(zz)
    if zz(i)~=0
        zze(length(zze)+1)= zz(i);
        num_zero=0;
    else
        if num_zero==0
            length_zze=length(zze);
            zze(length_zze+1:length_zze+2)=[0,0];
            num_zero=num_zero+1;
        else
            zze(length(zze))= num_zero;
            num_zero=num_zero+1;
        end
    end
    
    if mod(i,64)==0 && zz(i)==0
        length_zze=length(zze);
        zze(length_zze)=[];
        zze(length_zze-1)=EoB;
        num_zero=0;
    end
end
end

