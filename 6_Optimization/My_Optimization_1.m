%% initialization
clc;clear;close all;
addpath(strcat(pwd,'\functions'));
addpath(strcat(pwd,'\reference curves'));

%%
tic 
lena_small = double(imread('lena_small.tif'));
lena_small_ycbcr=ictRGB2YCbCr(lena_small);

% load image sequence
sequence=cell(1,21);
for i=1:21
    file_name=['foreman00',num2str(i+19),'.bmp'];
    sequence{i}=double(imread(file_name));
end

% scales = [0.07]; % quantization scale factor, the larger, the lower PSNR

%% video Codec
scales = [0.07 0.2 0.4 0.8 1.0 1.5 2 3 4 4.5]; % quantization scale factor, the larger, the lower PSNR
% scales = [4 4.5]; % quantization scale factor, the larger, the lower PSNR


fprintf('starting video Codec\n');
rate_mat_video=zeros(21,numel(scales));
psnr_mat_video=zeros(21,numel(scales));

for scaleIdx=1:numel(scales)
    rec_frame_ycbcr=cell(1,21);
    %% first frame
    qScale   = scales(scaleIdx);
    k_small  = IntraEncode(lena_small, qScale);
    k        = IntraEncode(sequence{1}, qScale);
    
    % train Huffman table
    pmf_ind=hist(k_small(:),min(k):max(k));
    pmf_ind=pmf_ind/sum(pmf_ind);
    [BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind);
    
    % encode the first image
    bytestream = enc_huffman_new(k-min(k)+1, BinCode, Codelengths);
    bitPerPixel = (numel(bytestream)*8) / (numel(sequence{1})/3);
    
    % decode the first image
    k_rec=double(dec_huffman_new(bytestream,BinaryTree,size(k(:),1)))+min(k)-1;
    I_rec = IntraDecode(k_rec, size(sequence{1}),qScale);  % rgb
    rec_frame_ycbcr{1}=ictRGB2YCbCr(I_rec);
    
    % calcualte bitrate and PSNR
    rate_mat_video(1,scaleIdx)=bitPerPixel;
    psnr_mat_video(1,scaleIdx)=calcPSNR(sequence{1},I_rec);
    
    %% following frames
    for i=2:21
        ref_frame_ycbcr=rec_frame_ycbcr{i-1};
        current_frame_rgb=sequence{i};
        current_frame_ycbcr=ictRGB2YCbCr(current_frame_rgb);
%         motion_vector=SSD(ref_frame_ycbcr(:,:,1),current_frame_ycbcr(:,:,1));
        motion_vector=SSD_HalfPel(ref_frame_ycbcr(:,:,1),current_frame_ycbcr(:,:,1));
%         motion_compensation_image_ycbcr=SSD_rec(ref_frame_ycbcr,motion_vector);
        motion_compensation_image_ycbcr=SSD_rec_HalfPel(ref_frame_ycbcr,motion_vector);
        residual_image_ycbcr=current_frame_ycbcr-motion_compensation_image_ycbcr;
        residual_image_rgb=ictYCbCr2RGB(residual_image_ycbcr);
        zeroRun_residual=IntraEncode(residual_image_rgb,qScale);  % residual image
        
        % train Huffman code
        if i==2
            % Huffman code for zero run
%             pmf_ind_zerorun=hist(zeroRun(:),min(zeroRun):max(zeroRun));
            pmf_ind_zerorun=hist(zeroRun_residual(:),-2000:4000);
            pmf_ind_zerorun=pmf_ind_zerorun/sum(pmf_ind_zerorun);
            [BinaryTree_zerorun,HuffCode_zerorun,BinCode_zerorun,Codelengths_zerorun]=buildHuffman(pmf_ind_zerorun);
            
            % Huffman code for motion vectors
%             pmf_ind_mv=hist(motion_vector(:),1:81);
            pmf_ind_mv=hist(motion_vector(:),1:289);
            pmf_ind_mv=pmf_ind_mv/sum(pmf_ind_mv);
            [BinaryTree_mv,HuffCode_mv,BinCode_mv,Codelengths_mv]=buildHuffman(pmf_ind_mv);
        end
        
        % encode Huffman
%         bytestream1=enc_huffman_new(zeroRun_residual-min(zeroRun_residual)+1, BinCode_zerorun, Codelengths_zerorun); % rate for zero run
        bytestream1=enc_huffman_new(zeroRun_residual+2001, BinCode_zerorun, Codelengths_zerorun); % rate for zero run  % Why 1001?
        bytestream2=enc_huffman_new(motion_vector(:)-min(motion_vector(:))+1,BinCode_mv,Codelengths_mv); % rate for motion vector
        bitPerPixel1 = (numel(bytestream1)*8) / (numel(sequence{i})/3);
        bitPerPixel2 = (numel(bytestream2)*8) / (numel(sequence{i})/3);
        rate_mat_video(i,scaleIdx)=bitPerPixel1+bitPerPixel2;
        
        % decode Huffman
        image_size=[size(motion_vector,1)*8,size(motion_vector,2)*8,3];
        decoded_frame_rgb=IntraDecode(zeroRun_residual,image_size,qScale)+ictYCbCr2RGB(motion_compensation_image_ycbcr);
        rec_frame_ycbcr{i}=ictRGB2YCbCr(decoded_frame_rgb);
        psnr_mat_video(i,scaleIdx)=calcPSNR(sequence{i},decoded_frame_rgb);
    end
    fprintf('Quantization scale: %.1f, bitrate: %.2f bits/pixel, PSNR: %.2fdB\n',qScale,mean(rate_mat_video(:,scaleIdx)),mean(psnr_mat_video(:,scaleIdx)));
end
rate_video=mean(rate_mat_video,1);
psnr_video=mean(psnr_mat_video,1);

% plot
figure;
plot(rate_video,psnr_video,'-d','LineWidth',1);hold on;
xlabel('bitrate [bit/pixel]','FontSize',16);
ylabel('PSNR [dB]','FontSize',16);

%% still image Codec
% scales_still=[0.15,0.3,0.7,1.0,1.5,3,5,7,10];
% fprintf('starting still image Codec\n');
% rate_mat_still=zeros(21,numel(scales_still));
% psnr_mat_still=zeros(21,numel(scales_still));
% 
% for scaleIdx=1:numel(scales_still)
%     
%     % train Huffman table
%     qScale  = scales_still(scaleIdx);
%     k_base  = IntraEncode(lena_small, qScale);
%     pmf_ind=hist(k_base(:),-1000:4000);
%     pmf_ind=pmf_ind/sum(pmf_ind);
%     [BinaryTree,HuffCode,BinCode,Codelengths]=buildHuffman(pmf_ind);
%     
%     for i=1:21
%     still_im_zerorun  = IntraEncode(sequence{i}, qScale);
%     
%     % encode the image
%     bytestream = enc_huffman_new(still_im_zerorun+1001, BinCode, Codelengths);
%     bitPerPixel = (numel(bytestream)*8) / (numel(sequence{i})/3);
%     
%     % decode the image
%     k_rec=double(dec_huffman_new(bytestream,BinaryTree,size(still_im_zerorun(:),1)))-1001;
%     I_rec = IntraDecode(k_rec, size(sequence{i}),qScale);  % rgb
% %     I_rec=ictRGB2YCbCr(I_rec);
%     
%     % calcualte bitrate and PSNR
%     rate_mat_still(i,scaleIdx)=bitPerPixel;
%     psnr_mat_still(i,scaleIdx)=calcPSNR(sequence{i},I_rec);
%     end
%     fprintf('Quantization scale: %.1f, bitrate: %.2f bits/pixel, PSNR: %.2fdB\n',qScale,mean(rate_mat_still(:,scaleIdx)),mean(psnr_mat_still(:,scaleIdx)));
% end
% rate_still=mean(rate_mat_still,1);
% psnr_still=mean(psnr_mat_still,1);
% plot(rate_still,psnr_still,'r-*');hold on;
% axis([0 4 20 45])
% legend('Video Codec','Still image Codec');
toc

%% reference curves comparision
hold on;
load('solution_values5.mat');
load('solution_values4.mat');
plot(bpp_solution, psnr_solution, '--+','LineWidth',1); % reference curve Chapter 5
plot(bpp_solution_ch4, psnr_solution_ch4, '--+','LineWidth',1); % reference curve Chapter 4
title('RD performance of Optimization, Foreman Sequence','FontSize',16);
legend('My Optimization1','Baseline2 (chapter5)','Baseline1 (chapter4)','FontSize',12);