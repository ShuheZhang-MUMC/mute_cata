clc
clear

path = genpath('func_set');
addpath(path);


save_path = 'retinal image\';

type = 2; % 2: ycbcr 1: CIElab
mult = 2;
intermediate = 1;
type_retinex = 'Gaussian';

%% load image
[img_name, img_path]=uigetfile('*.*');
raw_img0 = imread([img_path,img_name]);
[mm,nn,~] = size(raw_img0);
if sqrt(mm*nn) > 875 
    raw_img = double(imresize(raw_img0,875/nn))/255;
    raw_img0 = double(raw_img0)/255;
else
    raw_img = double(raw_img0)/255;
end

%% prepossing, padding the image, for removing the boundary effect
tic
name_ = img_name(1:end-4);
[m0,n0,~] = size(raw_img);
gauss_size = fix(sqrt(m0*n0)/20);
[raw_img,mask] = prepossing(raw_img,20,gauss_size,mult);
toc
imshow(raw_img,[])
imwrite(raw_img,'padded_test4.png')
tic
if type == 1
    img_ycbcr = rgb2lab(raw_img);
    q_channel = img_ycbcr(:,:,1)/100;
else
    img_ycbcr = rgb2ycbcr_ITU(raw_img);
    q_channel = img_ycbcr(:,:,1);
end


[m,n,~] = size(q_channel);

w_ = fix(sqrt(m*n)/150);    % window size




%% perform intensity correction
temp = img_ycbcr(:,:,1);
q_channel_new = q_channel;
[q_channel_new,illu_est_crse] = fun_illu_corr_coarse(q_channel,gauss_size,mask,type_retinex);
[q_channel_new,illu_est_fine] = fun_illu_corr_fine(q_channel_new,w_,mask);  


%% perform dehazing
I0_aim = 0.52; % target intensity for dehazing 
level = 5;
old_q = q_channel_new;
[q_channel_new,tran_est_fine] = ...
     fun_haze_corr_multiscale_TV_NAdam(q_channel_new,...
                                       I0_aim, ...
                                       0, ...
                                       level, ...
                                       mask, ...
                                       name_);

if ~exist([save_path,img_name(1:end-4)],'dir')
    mkdir([save_path,img_name(1:end-4)])
end

save([save_path,img_name(1:end-4),'//', ...
           img_name(1:end-4),'_',num2str(I0_aim),'out.mat'], ...
           'q_channel_new','tran_est_fine','old_q','q_channel','mask')

%% perform denoising
NN = [1,-2, 1;
     -2, 4,-2;
      1,-2, 1];
SN = imfilter(old_q,NN,'conv','replicate');
delta = sqrt(pi/2)*sum(abs(SN(:)))/(6*(size(old_q,1)-2)*(size(old_q,2)-2));
q_channel_new = solving_o_adam(old_q,tran_est_fine,delta/30,65);
img_ycbcr(:,:,1) = q_channel_new;


%% resize the image to the input one
if sqrt(mm*nn) > 875
    kk = mult*gauss_size;
    q_channel_new = q_channel_new.*mask;
    q_channel_new = q_channel_new(kk+1:end-kk,kk+1:end-kk,:);
    TT = q_channel_new./(temp(kk+1:end-kk,kk+1:end-kk,:) + eps);
    TT = imresize(TT,[mm,nn]);
    img_ycbcr0 = rgb2ycbcr_ITU(raw_img0);
    img_ycbcr0(:,:,1) = img_ycbcr0(:,:,1).*TT;
    out_img_cut0 = ycbcr2rgb_ITU(img_ycbcr0);
else
    out_img0 = ycbcr2rgb_ITU(img_ycbcr);
    out_img0 = out_img0 .* mask;
    kk = mult*gauss_size;
    out_img_cut0 = out_img0(kk+1:end-kk,kk+1:end-kk,:);  
end




imwrite(out_img_cut0,[save_path,img_name(1:end-4),'//', ...
    img_name(1:end-4),'_',num2str(I0_aim),'out.png'])


toc
rmpath(path);









