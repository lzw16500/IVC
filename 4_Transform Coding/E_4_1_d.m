load('foreman1_residual_zig_zag');
load('foreman2_residual_zig_zag');
load('foreman10_residual_zig_zag');
load('lena_small_zig_zag');

% test on zig-zag sequences

% convert the zig-zag scan into a 1xN vector
lena_small_zig_zag = lena_small_zig_zag(:)';
foreman1_residual_zig_zag = foreman1_residual_zig_zag(:)';
foreman2_residual_zig_zag = foreman2_residual_zig_zag(:)';
foreman10_residual_zig_zag = foreman10_residual_zig_zag(:)';

% perform ZeroRunEnc and ZeroRunDec on the sequences
zero_run_enc = ZeroRunEnc_EoB(lena_small_zig_zag);
zig_zag_result1 = ZeroRunDec_EoB(zero_run_enc);

zero_run_enc = ZeroRunEnc_EoB(foreman1_residual_zig_zag);
zig_zag_result2 = ZeroRunDec_EoB(zero_run_enc);

zero_run_enc = ZeroRunEnc_EoB(foreman2_residual_zig_zag);
zig_zag_result3 = ZeroRunDec_EoB(zero_run_enc);

zero_run_enc = ZeroRunEnc_EoB(foreman10_residual_zig_zag);
zig_zag_result4 = ZeroRunDec_EoB(zero_run_enc);

function zze = ZeroRunEnc_EoB(zz)
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)
EOB=1000;
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


function dst_eob = ZeroRunDec_EoB(src)
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)
EoB=1000;
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