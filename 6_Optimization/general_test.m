clc,clear; close all;
%%
% [X,Y] = meshgrid(1:4,1:3);
% V=[1 2 3 5;4 5 6 8;10 0 10 12];
% [Xq,Yq] = meshgrid(1:0.5:4,1:0.5:3);
% Vq = interp2(X,Y,V,Xq,Yq);

%%
% index_mat=zeros(17,17);
% temp=1;
% for i=1:17
%     for j=1:17
%         index_mat(i,j)=temp;
%         temp=temp+1;
%     end
% end

%%
% im=imread('foreman0020.bmp');
% [x,y,~]=size(im);
% [X,Y]=meshgrid(1:y,1:x);
% [Xq,Yq] = meshgrid(1:0.5:y,1:0.5:x);
% for i=1:3
%     ref_image_double=double(im(:,:,i));
%     ref_image_interp(:,:,i)=interp2(X,Y,ref_image_double,Xq,Yq); % double,[0,255]!!
% end
% imshow(ref_image_interp/max(ref_image_interp(:))); % convert the values back into [0,1]

%%
hold on;
load('solution_values5.mat');
load('solution_values4.mat');
plot(bpp_solution, psnr_solution, '--*','LineWidth',1); % reference curve Chapter 5
plot(bpp_solution_ch4, psnr_solution_ch4, '--+','LineWidth',1); % reference curve Chapter 4
title('RD performance of Optimization, Foreman Sequence','FontSize',18);
xlabel('bit/pixel','FontSize',18);
ylabel('PSNR [dB]');
legend('Baseline2 (chapter5)','Baseline1 (chapter4)','FontSize',12);
