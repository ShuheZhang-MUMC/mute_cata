function all_cor = color_adjustment(in1,in2,w,mask)


p = 0.3;

low_raw =  (Gaussian_filter(in2,w,1,'fft') .* mask);

low_new =  (Gaussian_filter(in1,w,1,'fft') .* mask);

hig_new = in1;
hig_new(mask) = in1(mask) ./ low_new(mask);

low_cor = low_new;

all_cor = hig_new .* low_cor;

apv1 = mean(mean(low_new(mask)));
apv2 = mean(mean(low_raw(mask)));

RPV = max(max(low_raw(mask))) - min(min(low_raw(mask)));

gamma1 = 1 - p;
gamma2 = 1 + p;
if abs(apv1-apv2)>1
    if apv1 < apv2
        gamma0 = gamma1;
    else
        gamma0 = gamma2;
    end
    low_new_gamma = low_new.^gamma0;
    
    temp1 =         low_new_gamma(mask)   - min(min(low_new_gamma(mask)));
    temp2 = max(max(low_new_gamma(mask))) - min(min(low_new_gamma(mask)));
    
    low_cor(mask) = temp1/temp2 * RPV + min(min(low_raw(mask)));
    apv1 = mean(mean(low_cor(mask)));
    for ccc = 1:20
        if abs(apv1-apv2)>1
            if gamma0 == gamma1
                if apv1<apv2
                    gamma0 = gamma0 - p;
                else
                    p = p/2;
                    gamma0 = gamma0 + p;
                end
            else
                if apv1<apv2
                    p = p/2;
                    gamma0 = gamma0 - p;
                else
                    gamma0 = gamma0 + p;
                end
            end
                low_new_gamma = low_new.^gamma0;
    
                temp1 =         low_new_gamma(mask)   - min(min(low_new_gamma(mask)));
                temp2 = max(max(low_new_gamma(mask))) - min(min(low_new_gamma(mask)));
    
                low_cor(mask) = temp1/temp2 * RPV + min(min(low_raw(mask)));
                apv1 = mean(mean(low_cor(mask)));
        else
            break;
        end
    end
    all_cor = hig_new .* low_cor;
else
    return;
end





end