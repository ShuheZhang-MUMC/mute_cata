function [illu_corr,raw_img_blur] = fun_illu_corr_coarse(input_img,gaussian_kernal_size,mask,type)


etp = 0.01;

%% low-pass filtering
if strcmp(type,'Gaussian')
    raw_img_blur = Gaussian_filter(input_img,gaussian_kernal_size,1,'fft');
elseif strcmp(type,'TV')
    raw_img_blur = TVL1denoise(input_img, mean(input_img(:))/50, 100);
else
    error("type should be 'Gaussian' or 'TV'")
end

log_illu_corr =  log(input_img + etp) - log(raw_img_blur + etp);
illu_corr = exp(log_illu_corr) - etp;

[~,~,c_chanel] = size(illu_corr);
%% grayvalue correction
for cc = 1:c_chanel
    temp = illu_corr(:,:,cc);
    max_ = max(max(temp(mask)));
    min_ = min(min(temp(mask)));
    
    temp = input_img(:,:,cc);
    max_0 = max(max(temp(mask)));
    min_0 = min(min(temp(mask)));
    illu_corr(:,:,cc) = (max_0 - min_0)*(illu_corr(:,:,cc) - min_) / (max_ - min_) + min_0;
end


end